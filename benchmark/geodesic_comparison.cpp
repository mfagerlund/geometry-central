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

  return 0;
}
