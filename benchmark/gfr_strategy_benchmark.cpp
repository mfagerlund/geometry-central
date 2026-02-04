// GFR Strategy Benchmark
// Compares GFR strategies (CornerGreedy, WedgeGreedy, CoherentMiniWedge) vs FlipOut
//
// Usage: gfr_strategy_benchmark <mesh.obj> [numPairs] [seed]

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

struct StrategySummary {
  std::string name;
  size_t gfrWins = 0;
  size_t flipoutWins = 0;
  size_t ties = 0;
  double meanDiff = 0;
  double totalRatio = 0;
};

static StrategySummary runStrategy(const std::string& name,
                                   ManifoldSurfaceMesh& mesh,
                                   VertexPositionGeometry& geometry,
                                   const std::vector<std::pair<size_t, size_t>>& pairs,
                                   const std::vector<double>& flipLengths,
                                   const FunnelGeodesicOptions& options) {
  StrategySummary summary;
  summary.name = name;

  double sumDiff = 0;
  double totalGfr = 0;
  double totalFo = 0;

  const double tol = 1e-8;
  for (size_t i = 0; i < pairs.size(); i++) {
    Vertex v0 = mesh.vertex(pairs[i].first);
    Vertex v1 = mesh.vertex(pairs[i].second);

    auto gfrPath = computeFunnelGeodesic(mesh, geometry, v0, v1, options);
    double gfrLen = gfrPath->length();
    double foLen = flipLengths[i];

    totalGfr += gfrLen;
    totalFo += foLen;

    double diff = (gfrLen - foLen) / foLen * 100.0;
    sumDiff += diff;

    if (gfrLen < foLen - tol) summary.gfrWins++;
    else if (foLen < gfrLen - tol) summary.flipoutWins++;
    else summary.ties++;
  }

  summary.meanDiff = sumDiff / pairs.size();
  summary.totalRatio = 100.0 * (totalGfr / totalFo - 1.0);
  return summary;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <mesh.obj> [numPairs] [seed]" << std::endl;
    return 1;
  }

  std::string meshPath = argv[1];
  size_t numPairs = (argc >= 3) ? std::stoul(argv[2]) : 1000;
  unsigned int seed = (argc >= 4) ? std::stoul(argv[3]) : 42;

  std::cout << "Loading mesh: " << meshPath << std::endl;
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);

  std::cout << "Mesh: " << mesh->nVertices() << " vertices, "
            << mesh->nFaces() << " faces" << std::endl;
  std::cout << "Testing " << numPairs << " random vertex pairs (seed=" << seed << ")" << std::endl;
  std::cout << std::endl;

  // Generate pairs
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> dist(0, mesh->nVertices() - 1);
  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.reserve(numPairs);
  for (size_t i = 0; i < numPairs; i++) {
    size_t v0 = dist(rng);
    size_t v1 = dist(rng);
    if (v0 == v1) { i--; continue; }
    pairs.emplace_back(v0, v1);
  }

  // Precompute FlipOut lengths
  std::vector<double> flipLengths;
  flipLengths.reserve(numPairs);
  std::cout << "Computing FlipOut baseline..." << std::endl;
  for (size_t i = 0; i < pairs.size(); i++) {
    Vertex v0 = mesh->vertex(pairs[i].first);
    Vertex v1 = mesh->vertex(pairs[i].second);
    auto flipNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, v0, v1);
    flipNetwork->iterativeShorten();
    flipLengths.push_back(flipNetwork->length());
  }
  std::cout << "FlipOut baseline ready." << std::endl << std::endl;

  // Run strategies
  std::vector<StrategySummary> summaries;

  FunnelGeodesicOptions cornerOpts;
  cornerOpts.strategy = FlipStrategy::CornerGreedy;
  summaries.push_back(runStrategy("CornerGreedy", *mesh, *geometry, pairs, flipLengths, cornerOpts));

  FunnelGeodesicOptions wedgeOpts;
  wedgeOpts.strategy = FlipStrategy::WedgeGreedy;
  summaries.push_back(runStrategy("WedgeGreedy", *mesh, *geometry, pairs, flipLengths, wedgeOpts));

  for (size_t depth : {1, 2, 3, 4}) {
    FunnelGeodesicOptions miniOpts;
    miniOpts.strategy = FlipStrategy::CoherentMiniWedge;
    miniOpts.coherentMiniWedgeDepth = depth;
    summaries.push_back(runStrategy("CoherentMiniWedge d=" + std::to_string(depth),
                                   *mesh, *geometry, pairs, flipLengths, miniOpts));
  }

  // Print summary
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "=== GFR Strategy Comparison vs FlipOut ===" << std::endl;
  for (const auto& s : summaries) {
    std::cout << s.name << ":" << std::endl;
    std::cout << "  GFR wins: " << s.gfrWins << " (" << (100.0 * s.gfrWins / pairs.size()) << "%)" << std::endl;
    std::cout << "  FlipOut wins: " << s.flipoutWins << " (" << (100.0 * s.flipoutWins / pairs.size()) << "%)" << std::endl;
    std::cout << "  Ties: " << s.ties << " (" << (100.0 * s.ties / pairs.size()) << "%)" << std::endl;
    std::cout << "  Mean diff: " << s.meanDiff << "%"<< std::endl;
    std::cout << "  Total distance ratio: " << s.totalRatio << "%" << std::endl;
    std::cout << std::endl;
  }

  return 0;
}
