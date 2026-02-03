// Geodesic Comparison Benchmark
// Compares Funnel Geodesics (GFR) vs FlipOut on various meshes
//
// Usage: geodesic_comparison <mesh.obj> [numPairs] [seed]

#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct ComparisonResult {
  size_t v0, v1;
  double gfrLength;
  double flipoutLength;
  size_t gfrIterations;
  size_t flipoutIterations;
  double gfrTimeMs;
  double flipoutTimeMs;
};

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <mesh.obj> [numPairs] [seed]" << std::endl;
    return 1;
  }

  std::string meshPath = argv[1];
  size_t numPairs = (argc >= 3) ? std::stoul(argv[2]) : 100;
  unsigned int seed = (argc >= 4) ? std::stoul(argv[3]) : 42;

  // Load mesh
  std::cout << "Loading mesh: " << meshPath << std::endl;
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);

  std::cout << "Mesh: " << mesh->nVertices() << " vertices, "
            << mesh->nFaces() << " faces" << std::endl;
  std::cout << "Testing " << numPairs << " random vertex pairs (seed=" << seed << ")" << std::endl;
  std::cout << std::endl;

  // Generate random vertex pairs
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> dist(0, mesh->nVertices() - 1);

  std::vector<std::pair<size_t, size_t>> vertexPairs;
  for (size_t i = 0; i < numPairs; i++) {
    size_t v0 = dist(rng);
    size_t v1 = dist(rng);
    if (v0 != v1) {
      vertexPairs.emplace_back(v0, v1);
    } else {
      i--; // retry
    }
  }

  // Run comparison
  std::vector<ComparisonResult> results;
  results.reserve(numPairs);

  double totalGfrTime = 0;
  double totalFlipoutTime = 0;
  double totalGfrLength = 0;
  double totalFlipoutLength = 0;
  size_t totalGfrIters = 0;
  size_t totalFlipoutIters = 0;
  size_t gfrWins = 0;
  size_t flipoutWins = 0;
  size_t ties = 0;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Progress: ";

  for (size_t i = 0; i < vertexPairs.size(); i++) {
    Vertex v0 = mesh->vertex(vertexPairs[i].first);
    Vertex v1 = mesh->vertex(vertexPairs[i].second);

    ComparisonResult r;
    r.v0 = vertexPairs[i].first;
    r.v1 = vertexPairs[i].second;

    // Run GFR
    auto gfrStart = std::chrono::high_resolution_clock::now();
    auto gfrPath = computeFunnelGeodesic(*mesh, *geometry, v0, v1);
    auto gfrEnd = std::chrono::high_resolution_clock::now();
    r.gfrLength = gfrPath->length();
    r.gfrIterations = gfrPath->iterationCount();
    r.gfrTimeMs = std::chrono::duration<double, std::milli>(gfrEnd - gfrStart).count();

    // Run FlipOut
    auto flipStart = std::chrono::high_resolution_clock::now();
    auto flipNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, v0, v1);
    flipNetwork->iterativeShorten();
    auto flipEnd = std::chrono::high_resolution_clock::now();
    r.flipoutLength = flipNetwork->length();
    r.flipoutIterations = flipNetwork->nShortenIters;
    r.flipoutTimeMs = std::chrono::duration<double, std::milli>(flipEnd - flipStart).count();

    results.push_back(r);

    // Accumulate totals
    totalGfrTime += r.gfrTimeMs;
    totalFlipoutTime += r.flipoutTimeMs;
    totalGfrLength += r.gfrLength;
    totalFlipoutLength += r.flipoutLength;
    totalGfrIters += r.gfrIterations;
    totalFlipoutIters += r.flipoutIterations;

    // Compare lengths (within tolerance)
    const double tol = 1e-8;
    if (r.gfrLength < r.flipoutLength - tol) {
      gfrWins++;
    } else if (r.flipoutLength < r.gfrLength - tol) {
      flipoutWins++;
    } else {
      ties++;
    }

    // Progress indicator
    if ((i + 1) % 10 == 0) {
      std::cout << "." << std::flush;
    }
  }
  std::cout << " done!" << std::endl << std::endl;

  // Print summary
  std::cout << "============ SUMMARY ============" << std::endl;
  std::cout << std::endl;
  std::cout << "Path Length Comparison:" << std::endl;
  std::cout << "  GFR shorter:      " << gfrWins << " (" << (100.0 * gfrWins / numPairs) << "%)" << std::endl;
  std::cout << "  FlipOut shorter:  " << flipoutWins << " (" << (100.0 * flipoutWins / numPairs) << "%)" << std::endl;
  std::cout << "  Equal (tie):      " << ties << " (" << (100.0 * ties / numPairs) << "%)" << std::endl;
  std::cout << std::endl;

  std::cout << "Average Path Length:" << std::endl;
  std::cout << "  GFR:      " << (totalGfrLength / numPairs) << std::endl;
  std::cout << "  FlipOut:  " << (totalFlipoutLength / numPairs) << std::endl;
  std::cout << std::endl;

  std::cout << "Average Iterations:" << std::endl;
  std::cout << "  GFR:      " << (totalGfrIters / (double)numPairs) << std::endl;
  std::cout << "  FlipOut:  " << (totalFlipoutIters / (double)numPairs) << std::endl;
  std::cout << std::endl;

  std::cout << "Total Time (ms):" << std::endl;
  std::cout << "  GFR:      " << totalGfrTime << std::endl;
  std::cout << "  FlipOut:  " << totalFlipoutTime << std::endl;
  std::cout << std::endl;

  std::cout << "Average Time per Path (ms):" << std::endl;
  std::cout << "  GFR:      " << (totalGfrTime / numPairs) << std::endl;
  std::cout << "  FlipOut:  " << (totalFlipoutTime / numPairs) << std::endl;
  std::cout << std::endl;

  double speedup = totalFlipoutTime / totalGfrTime;
  std::cout << "Speedup (FlipOut time / GFR time): " << speedup << "x" << std::endl;
  std::cout << std::endl;

  // Check for cases where GFR is longer (potential bugs)
  size_t gfrLonger = 0;
  double maxDiff = 0;
  for (const auto& r : results) {
    if (r.gfrLength > r.flipoutLength + 1e-6) {
      gfrLonger++;
      maxDiff = std::max(maxDiff, r.gfrLength - r.flipoutLength);
    }
  }

  if (gfrLonger > 0) {
    std::cout << "WARNING: GFR produced longer path in " << gfrLonger << " cases" << std::endl;
    std::cout << "  Max difference: " << maxDiff << std::endl;
  }

  // Detailed win margin analysis
  std::cout << std::endl;
  std::cout << "============ WIN MARGIN ANALYSIS ============" << std::endl;
  std::cout << std::endl;

  // Calculate percentage differences and bucket them
  // Negative = GFR shorter, Positive = FlipOut shorter
  std::vector<double> percentDiffs;
  for (const auto& r : results) {
    double avg = (r.gfrLength + r.flipoutLength) / 2.0;
    double pctDiff = 100.0 * (r.gfrLength - r.flipoutLength) / avg;
    percentDiffs.push_back(pctDiff);
  }

  // Sort for percentile analysis
  std::vector<double> sorted = percentDiffs;
  std::sort(sorted.begin(), sorted.end());

  // Histogram buckets
  size_t gfrBig = 0;     // GFR wins by > 1%
  size_t gfrMed = 0;     // GFR wins by 0.1-1%
  size_t gfrSmall = 0;   // GFR wins by 0.01-0.1%
  size_t gfrTiny = 0;    // GFR wins by 0.001-0.01%
  size_t tieZone = 0;    // Within ±0.001%
  size_t foTiny = 0;     // FlipOut wins by 0.001-0.01%
  size_t foSmall = 0;    // FlipOut wins by 0.01-0.1%
  size_t foMed = 0;      // FlipOut wins by 0.1-1%
  size_t foBig = 0;      // FlipOut wins by > 1%

  for (double d : percentDiffs) {
    if (d < -1.0) gfrBig++;
    else if (d < -0.1) gfrMed++;
    else if (d < -0.01) gfrSmall++;
    else if (d < -0.001) gfrTiny++;
    else if (d <= 0.001) tieZone++;
    else if (d < 0.01) foTiny++;
    else if (d < 0.1) foSmall++;
    else if (d < 1.0) foMed++;
    else foBig++;
  }

  std::cout << "Distribution (% difference, negative = GFR shorter):" << std::endl;
  std::cout << "  GFR wins by >1%:       " << std::setw(5) << gfrBig << " (" << std::setw(5) << std::fixed << std::setprecision(1) << (100.0 * gfrBig / numPairs) << "%)" << std::endl;
  std::cout << "  GFR wins by 0.1-1%:    " << std::setw(5) << gfrMed << " (" << std::setw(5) << (100.0 * gfrMed / numPairs) << "%)" << std::endl;
  std::cout << "  GFR wins by 0.01-0.1%: " << std::setw(5) << gfrSmall << " (" << std::setw(5) << (100.0 * gfrSmall / numPairs) << "%)" << std::endl;
  std::cout << "  GFR wins by <0.01%:    " << std::setw(5) << gfrTiny << " (" << std::setw(5) << (100.0 * gfrTiny / numPairs) << "%)" << std::endl;
  std::cout << "  Tie zone (±0.001%):    " << std::setw(5) << tieZone << " (" << std::setw(5) << (100.0 * tieZone / numPairs) << "%)" << std::endl;
  std::cout << "  FlipOut wins by <0.01%:" << std::setw(5) << foTiny << " (" << std::setw(5) << (100.0 * foTiny / numPairs) << "%)" << std::endl;
  std::cout << "  FlipOut wins by 0.01-0.1%:" << std::setw(4) << foSmall << " (" << std::setw(5) << (100.0 * foSmall / numPairs) << "%)" << std::endl;
  std::cout << "  FlipOut wins by 0.1-1%:" << std::setw(5) << foMed << " (" << std::setw(5) << (100.0 * foMed / numPairs) << "%)" << std::endl;
  std::cout << "  FlipOut wins by >1%:   " << std::setw(5) << foBig << " (" << std::setw(5) << (100.0 * foBig / numPairs) << "%)" << std::endl;
  std::cout << std::endl;

  // Percentiles
  std::cout << "Percentiles (% difference):" << std::endl;
  std::cout << "  Min (best GFR):  " << std::setw(8) << std::setprecision(4) << sorted.front() << "%" << std::endl;
  std::cout << "  5th percentile:  " << std::setw(8) << sorted[numPairs * 5 / 100] << "%" << std::endl;
  std::cout << "  25th percentile: " << std::setw(8) << sorted[numPairs * 25 / 100] << "%" << std::endl;
  std::cout << "  Median:          " << std::setw(8) << sorted[numPairs / 2] << "%" << std::endl;
  std::cout << "  75th percentile: " << std::setw(8) << sorted[numPairs * 75 / 100] << "%" << std::endl;
  std::cout << "  95th percentile: " << std::setw(8) << sorted[numPairs * 95 / 100] << "%" << std::endl;
  std::cout << "  Max (best FO):   " << std::setw(8) << sorted.back() << "%" << std::endl;
  std::cout << std::endl;

  // Mean difference
  double sumDiff = 0;
  for (double d : percentDiffs) sumDiff += d;
  double meanDiff = sumDiff / numPairs;
  std::cout << "Mean difference: " << std::setprecision(4) << meanDiff << "% (negative = GFR shorter overall)" << std::endl;
  std::cout << std::endl;

  // Iteration histogram
  std::cout << "GFR Iteration Distribution:" << std::endl;
  size_t iter0 = 0, iter1_5 = 0, iter6_20 = 0, iter21_50 = 0, iter50plus = 0;
  for (const auto& r : results) {
    if (r.gfrIterations == 0) iter0++;
    else if (r.gfrIterations <= 5) iter1_5++;
    else if (r.gfrIterations <= 20) iter6_20++;
    else if (r.gfrIterations <= 50) iter21_50++;
    else iter50plus++;
  }
  std::cout << "  0 iterations:     " << std::setw(5) << iter0 << " (" << std::setw(5) << (100.0 * iter0 / numPairs) << "%)" << std::endl;
  std::cout << "  1-5 iterations:   " << std::setw(5) << iter1_5 << " (" << std::setw(5) << (100.0 * iter1_5 / numPairs) << "%)" << std::endl;
  std::cout << "  6-20 iterations:  " << std::setw(5) << iter6_20 << " (" << std::setw(5) << (100.0 * iter6_20 / numPairs) << "%)" << std::endl;
  std::cout << "  21-50 iterations: " << std::setw(5) << iter21_50 << " (" << std::setw(5) << (100.0 * iter21_50 / numPairs) << "%)" << std::endl;
  std::cout << "  50+ iterations:   " << std::setw(5) << iter50plus << " (" << std::setw(5) << (100.0 * iter50plus / numPairs) << "%)" << std::endl;
  std::cout << std::endl;

  // Analyze 0-iteration paths specifically
  if (iter0 > 0) {
    std::cout << "Zero-Iteration Path Analysis:" << std::endl;
    size_t zeroIterGfrWins = 0, zeroIterFoWins = 0, zeroIterTies = 0;
    double zeroIterTotalDiff = 0;
    for (const auto& r : results) {
      if (r.gfrIterations == 0) {
        double avg = (r.gfrLength + r.flipoutLength) / 2.0;
        double pctDiff = 100.0 * (r.gfrLength - r.flipoutLength) / avg;
        zeroIterTotalDiff += pctDiff;
        if (r.gfrLength < r.flipoutLength - 1e-8) zeroIterGfrWins++;
        else if (r.flipoutLength < r.gfrLength - 1e-8) zeroIterFoWins++;
        else zeroIterTies++;
      }
    }
    std::cout << "  GFR wins:     " << zeroIterGfrWins << " (" << (100.0 * zeroIterGfrWins / iter0) << "%)" << std::endl;
    std::cout << "  FlipOut wins: " << zeroIterFoWins << " (" << (100.0 * zeroIterFoWins / iter0) << "%)" << std::endl;
    std::cout << "  Ties:         " << zeroIterTies << " (" << (100.0 * zeroIterTies / iter0) << "%)" << std::endl;
    std::cout << "  Mean diff:    " << std::setprecision(4) << (zeroIterTotalDiff / iter0) << "%" << std::endl;

    // Show the 0-iter paths where GFR wins big
    std::cout << "  Top GFR wins (0-iter):" << std::endl;
    std::vector<std::pair<double, size_t>> zeroIterDiffs;
    for (size_t i = 0; i < results.size(); i++) {
      if (results[i].gfrIterations == 0) {
        double avg = (results[i].gfrLength + results[i].flipoutLength) / 2.0;
        double pct = 100.0 * (results[i].gfrLength - results[i].flipoutLength) / avg;
        zeroIterDiffs.emplace_back(pct, i);
      }
    }
    std::sort(zeroIterDiffs.begin(), zeroIterDiffs.end());
    for (size_t i = 0; i < std::min(size_t(3), zeroIterDiffs.size()); i++) {
      const auto& r = results[zeroIterDiffs[i].second];
      std::cout << "    V" << r.v0 << " -> V" << r.v1
                << ": GFR=" << r.gfrLength << ", FO=" << r.flipoutLength
                << " (" << std::showpos << zeroIterDiffs[i].first << "%)" << std::noshowpos
                << std::endl;
    }
    std::cout << std::endl;
  }

  // Show worst 5 paths where FlipOut wins significantly
  std::cout << "Top 5 paths where FlipOut wins most:" << std::endl;
  std::vector<std::pair<double, size_t>> worstForGfr;
  for (size_t i = 0; i < results.size(); i++) {
    double avg = (results[i].gfrLength + results[i].flipoutLength) / 2.0;
    double pct = 100.0 * (results[i].gfrLength - results[i].flipoutLength) / avg;
    worstForGfr.emplace_back(pct, i);
  }
  std::sort(worstForGfr.rbegin(), worstForGfr.rend());
  for (size_t i = 0; i < std::min(size_t(5), worstForGfr.size()); i++) {
    const auto& r = results[worstForGfr[i].second];
    std::cout << "  V" << r.v0 << " -> V" << r.v1
              << ": GFR=" << std::setprecision(4) << r.gfrLength
              << ", FO=" << r.flipoutLength
              << " (+" << worstForGfr[i].first << "% longer)"
              << " GFR iters=" << r.gfrIterations << ", FO iters=" << r.flipoutIterations
              << std::endl;
  }

  return 0;
}
