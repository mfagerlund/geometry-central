// Epsilon vs FlipOut Test: Does accepting worsening flips help beat FlipOut on more paths?
//
// For each path, compares GFR (with various epsilon values) against FlipOut.
// Counts wins/ties/losses to see if relaxed epsilon helps on difficult paths.
//
// Usage: epsilon_vs_flipout_test [--mesh FILE] [--paths N]

#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/funnel_geodesics.h"
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

// Custom straightener with parameterized epsilon (same as epsilon_test.cpp)
double straightenWithEpsilon(
    std::vector<Face> sleeveFaces,
    Vertex startVert, Vertex endVert,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    double epsilon,
    bool allowAnyFlip
) {
  auto flatPos = funnel_internal::flattenSleeve(sleeveFaces, startVert, geom);
  auto portals = funnel_internal::buildPortals(sleeveFaces, flatPos);
  Vector2 entry2D = flatPos[startVert];
  Vector2 exit2D = flatPos[endVert];
  auto funnel = funnel_internal::runFunnel(portals, entry2D, exit2D);

  const size_t maxIters = 20000;
  std::vector<uint32_t> rejectedPhase(mesh.nVertices(), 0);
  const double relativeNegligible = 1e-5;

  uint32_t phase = 1;
  size_t iteration = 0;

  auto corners = funnel_internal::analyzeCorners(sleeveFaces, funnel, flatPos);

  auto selectCornerToFlip = [&]() -> funnel_internal::WaypointCorner* {
    funnel_internal::WaypointCorner* best = nullptr;
    double bestError = 0;
    for (auto& c : corners) {
      if (!c.wantsToFlip()) continue;
      if (rejectedPhase[c.vertex.getIndex()] == phase) continue;
      if (c.angleErrorDeg > bestError) {
        best = &c;
        bestError = c.angleErrorDeg;
      }
    }
    return best;
  };

  auto hasFlippableCorners = [&]() {
    for (const auto& c : corners) {
      if (c.wantsToFlip() && rejectedPhase[c.vertex.getIndex()] != phase) {
        return true;
      }
    }
    return false;
  };

  if (!hasFlippableCorners()) {
    return funnel.distance;
  }

  while (iteration < maxIters) {
    phase++;
    double phaseStartDistance = funnel.distance;

    while (iteration < maxIters) {
      iteration++;

      auto* cornerToFlip = selectCornerToFlip();
      if (cornerToFlip == nullptr) break;

      auto action = funnel_internal::computeFlipAction(sleeveFaces, *cornerToFlip, mesh);
      if (!action.canFlip) {
        rejectedPhase[cornerToFlip->vertex.getIndex()] = phase;
        continue;
      }

      size_t idx = cornerToFlip->vertex.getIndex();
      double distanceBefore = funnel.distance;

      auto newFaces = funnel_internal::applyFlip(sleeveFaces, action);
      auto newFlatPos = funnel_internal::flattenSleeve(newFaces, startVert, geom);
      auto newPortals = funnel_internal::buildPortals(newFaces, newFlatPos);
      auto newFunnel = funnel_internal::runFunnel(newPortals, entry2D, newFlatPos[endVert]);

      double actualImprovement = distanceBefore - newFunnel.distance;
      double relativeImprovement = actualImprovement / distanceBefore;

      bool acceptFlip;
      if (allowAnyFlip) {
        acceptFlip = true;
        rejectedPhase[idx] = phase;
      } else {
        acceptFlip = relativeImprovement > epsilon;
      }

      if (acceptFlip) {
        sleeveFaces = std::move(newFaces);
        flatPos = std::move(newFlatPos);
        portals = std::move(newPortals);
        funnel = std::move(newFunnel);
        exit2D = flatPos[endVert];

        corners = funnel_internal::analyzeCorners(sleeveFaces, funnel, flatPos);

        if (!allowAnyFlip && relativeImprovement < relativeNegligible) {
          rejectedPhase[idx] = phase;
        }
      } else {
        rejectedPhase[idx] = phase;
      }
    }

    double phaseEndDistance = funnel.distance;
    double phaseImprovement = (phaseStartDistance - phaseEndDistance) / phaseStartDistance;
    bool improved = phaseImprovement > relativeNegligible;

    corners = funnel_internal::analyzeCorners(sleeveFaces, funnel, flatPos);
    bool hasCorners = false;
    for (const auto& c : corners) {
      if (c.wantsToFlip()) {
        hasCorners = true;
        break;
      }
    }

    if (!improved || !hasCorners) break;
  }

  return funnel.distance;
}

std::string getStem(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  std::string filename = (pos == std::string::npos) ? path : path.substr(pos + 1);
  pos = filename.find_last_of('.');
  return (pos == std::string::npos) ? filename : filename.substr(0, pos);
}

int main(int argc, char** argv) {
  std::string meshPath = "C:/Dev/Colonel/Data/Meshes/stanford-bunny.obj";
  size_t numPaths = 200;
  unsigned int seed = 42;

  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--mesh" && i + 1 < argc) {
      meshPath = argv[++i];
    } else if (arg == "--paths" && i + 1 < argc) {
      numPaths = std::stoul(argv[++i]);
    } else if (arg == "--seed" && i + 1 < argc) {
      seed = std::stoul(argv[++i]);
    }
  }

  std::cerr << "Loading " << getStem(meshPath) << "..." << std::endl;

  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);

  std::cerr << "  " << mesh->nVertices() << " vertices" << std::endl;

  // Generate random vertex pairs
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> dist(0, mesh->nVertices() - 1);

  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.reserve(numPaths);
  while (pairs.size() < numPaths) {
    size_t v0 = dist(rng);
    size_t v1 = dist(rng);
    if (v0 != v1) pairs.emplace_back(v0, v1);
  }

  // Test configurations
  struct TestConfig {
    double epsilon;
    bool allowAny;
    std::string name;
  };

  std::vector<TestConfig> configs = {
    {0.0, false, "epsilon=0"},
    {-0.01, false, "epsilon=-1%"},
    {0.0, true, "allowAny"}
  };

  std::cout << "# Epsilon vs FlipOut: " << getStem(meshPath) << std::endl;
  std::cout << "Paths: " << numPaths << ", Seed: " << seed << std::endl;
  std::cout << std::endl;

  // For each config, track wins/ties/losses vs FlipOut
  struct ConfigResult {
    std::string name;
    size_t wins = 0;    // GFR shorter
    size_t ties = 0;    // Within 0.001%
    size_t losses = 0;  // FlipOut shorter
    double totalGfr = 0;
    double totalFlipout = 0;
    size_t flippedFromLoss = 0;  // Paths that were losses with baseline but wins with this config
  };

  std::vector<ConfigResult> results(configs.size());
  for (size_t i = 0; i < configs.size(); i++) {
    results[i].name = configs[i].name;
  }

  // Also track baseline results per path for comparison
  std::vector<double> baselineGfr(numPaths);
  std::vector<double> flipoutLengths(numPaths);
  std::vector<bool> baselineLost(numPaths);  // True if baseline lost to FlipOut

  std::cerr << "Computing paths..." << std::endl;

  size_t pathIdx = 0;
  for (const auto& pair : pairs) {
    Vertex v0 = mesh->vertex(pair.first);
    Vertex v1 = mesh->vertex(pair.second);

    // Get FlipOut length
    auto flipNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, v0, v1);
    flipNetwork->iterativeShorten();
    double flipoutLength = flipNetwork->length();
    flipoutLengths[pathIdx] = flipoutLength;

    // Get initial sleeve
    std::vector<Face> sleeveFaces = funnel_internal::buildFaceStripVeryDiscrete(*mesh, *geometry, v0, v1);
    if (sleeveFaces.empty()) {
      pathIdx++;
      continue;
    }

    // Test each epsilon config
    for (size_t i = 0; i < configs.size(); i++) {
      double gfrLength = straightenWithEpsilon(
          sleeveFaces, v0, v1, *mesh, *geometry, configs[i].epsilon, configs[i].allowAny);

      results[i].totalGfr += gfrLength;
      results[i].totalFlipout += flipoutLength;

      // Compare (tie threshold: 0.001% = 1e-5 relative)
      double diff = (gfrLength - flipoutLength) / flipoutLength;
      if (diff < -1e-5) {
        results[i].wins++;
      } else if (diff > 1e-5) {
        results[i].losses++;
      } else {
        results[i].ties++;
      }

      // Track baseline for comparison
      if (i == 0) {
        baselineGfr[pathIdx] = gfrLength;
        baselineLost[pathIdx] = (diff > 1e-5);
      } else if (baselineLost[pathIdx]) {
        // Check if this config converted a loss to a win
        double newDiff = (gfrLength - flipoutLength) / flipoutLength;
        if (newDiff < -1e-5) {
          results[i].flippedFromLoss++;
        }
      }
    }

    pathIdx++;
    if (pathIdx % 50 == 0) {
      std::cerr << "  " << pathIdx << "/" << numPaths << std::endl;
    }
  }

  // Print results
  std::cout << "| Config | Wins | Ties | Losses | Win% | Flipped→Win | Avg Diff |" << std::endl;
  std::cout << "|--------|------|------|--------|------|-------------|----------|" << std::endl;

  for (const auto& r : results) {
    size_t total = r.wins + r.ties + r.losses;
    double winPct = 100.0 * r.wins / total;
    double avgDiff = 100.0 * (r.totalGfr - r.totalFlipout) / r.totalFlipout;

    std::cout << "| " << r.name;
    std::cout << " | " << r.wins;
    std::cout << " | " << r.ties;
    std::cout << " | " << r.losses;
    std::cout << " | " << std::fixed << std::setprecision(1) << winPct << "%";
    std::cout << " | " << r.flippedFromLoss;
    std::cout << " | " << std::fixed << std::setprecision(3) << avgDiff << "%";
    std::cout << " |" << std::endl;
  }

  std::cout << std::endl;
  std::cout << "**Flipped→Win**: Paths where baseline (epsilon=0) lost to FlipOut, but this config won." << std::endl;

  return 0;
}
