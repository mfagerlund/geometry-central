// Epsilon Effect Test: What happens with different acceptance thresholds?
//
// Tests:
// 1. epsilon = 0 (current behavior - accept any positive improvement)
// 2. epsilon = -0.001 (accept 0.1% worsening)
// 3. epsilon = -0.01 (accept 1% worsening)
// 4. allowAny: accept ALL flips, but mark vertex done after flip to prevent oscillation
//
// Usage: epsilon_test [--mesh FILE] [--paths N]

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

struct EpsilonTestResult {
  std::string name;
  size_t numPaths;
  size_t totalIterations;
  size_t totalFlipsAccepted;
  size_t totalFlipsRejected;
  double totalInitialLength;
  double totalFinalLength;
  double totalTimeMs;
  size_t maxIterationsOnPath;
  size_t worseningFlipsAccepted;  // Flips that made path longer

  double avgIterations() const { return (double)totalIterations / numPaths; }
  double avgImprovement() const { return (totalInitialLength - totalFinalLength) / totalInitialLength * 100; }
  double avgTimeMs() const { return totalTimeMs / numPaths; }
};

// Custom straightener with parameterized epsilon
struct CustomStraightenResult {
  double finalDistance;
  size_t iterations;
  size_t flipsAccepted;
  size_t flipsRejected;
  size_t worseningFlipsAccepted;
};

CustomStraightenResult straightenWithEpsilon(
    std::vector<Face> sleeveFaces,  // Copy, we'll modify it
    Vertex startVert, Vertex endVert,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    double epsilon,      // Negative = accept relative worsening (e.g., -0.01 = accept 1% worse)
    bool allowAnyFlip    // If true, accept ALL flips but mark vertex done
) {
  CustomStraightenResult result = {0, 0, 0, 0, 0};

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
    result.finalDistance = funnel.distance;
    return result;
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
        // Accept ANY flip (even worsening), but mark vertex done for this phase
        acceptFlip = true;
        rejectedPhase[idx] = phase;  // Prevent re-flipping this vertex in same phase
      } else {
        // Accept based on epsilon threshold (relative)
        // epsilon=0: accept if improvement > 0
        // epsilon=-0.01: accept if relativeImprovement > -0.01 (up to 1% worsening OK)
        acceptFlip = relativeImprovement > epsilon;
      }

      if (acceptFlip) {
        if (actualImprovement < 0) {
          result.worseningFlipsAccepted++;
        }

        sleeveFaces = std::move(newFaces);
        flatPos = std::move(newFlatPos);
        portals = std::move(newPortals);
        funnel = std::move(newFunnel);
        exit2D = flatPos[endVert];
        result.iterations++;
        result.flipsAccepted++;

        corners = funnel_internal::analyzeCorners(sleeveFaces, funnel, flatPos);

        // For normal mode, mark tiny improvements as rejected to prevent infinite loops
        if (!allowAnyFlip && relativeImprovement < relativeNegligible) {
          rejectedPhase[idx] = phase;
        }
      } else {
        rejectedPhase[idx] = phase;
        result.flipsRejected++;
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

  result.finalDistance = funnel.distance;
  return result;
}

std::string getFilename(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

std::string getStem(const std::string& path) {
  std::string filename = getFilename(path);
  size_t pos = filename.find_last_of('.');
  return (pos == std::string::npos) ? filename : filename.substr(0, pos);
}

int main(int argc, char** argv) {
  std::string meshPath = "C:/Dev/Colonel/Data/Meshes/stanford-bunny.obj";
  size_t numPaths = 100;
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
    {0.0, false, "epsilon=0 (baseline: accept any improvement)"},
    {-0.001, false, "epsilon=-0.1% (accept up to 0.1% worsening)"},
    {-0.01, false, "epsilon=-1% (accept up to 1% worsening)"},
    {-0.1, false, "epsilon=-10% (accept up to 10% worsening)"},
    {0.0, true, "allowAny (accept ALL flips, 1 flip/vertex/phase)"}
  };

  std::cout << "# Epsilon Effect Test: " << getStem(meshPath) << std::endl;
  std::cout << "Paths: " << numPaths << ", Seed: " << seed << std::endl;
  std::cout << std::endl;

  std::vector<EpsilonTestResult> allResults;

  for (const auto& config : configs) {
    EpsilonTestResult result;
    result.name = config.name;
    result.numPaths = 0;
    result.totalIterations = 0;
    result.totalFlipsAccepted = 0;
    result.totalFlipsRejected = 0;
    result.totalInitialLength = 0;
    result.totalFinalLength = 0;
    result.totalTimeMs = 0;
    result.maxIterationsOnPath = 0;
    result.worseningFlipsAccepted = 0;

    std::cerr << "Testing: " << config.name << "..." << std::endl;

    for (const auto& pair : pairs) {
      Vertex v0 = mesh->vertex(pair.first);
      Vertex v1 = mesh->vertex(pair.second);

      // Get initial sleeve from VDG (same as normal GFR)
      std::vector<Face> sleeveFaces = funnel_internal::buildFaceStripVeryDiscrete(
          *mesh, *geometry, v0, v1);
      if (sleeveFaces.empty()) continue;

      result.numPaths++;

      // Get initial funnel length
      auto flatPos = funnel_internal::flattenSleeve(sleeveFaces, v0, *geometry);
      auto portals = funnel_internal::buildPortals(sleeveFaces, flatPos);
      auto funnel = funnel_internal::runFunnel(portals, flatPos[v0], flatPos[v1]);
      double initialLength = funnel.distance;
      result.totalInitialLength += initialLength;

      // Run straightening with this epsilon
      auto start = std::chrono::high_resolution_clock::now();
      auto straightenResult = straightenWithEpsilon(
        sleeveFaces, v0, v1, *mesh, *geometry, config.epsilon, config.allowAny);
      auto end = std::chrono::high_resolution_clock::now();

      result.totalFinalLength += straightenResult.finalDistance;
      result.totalIterations += straightenResult.iterations;
      result.totalFlipsAccepted += straightenResult.flipsAccepted;
      result.totalFlipsRejected += straightenResult.flipsRejected;
      result.worseningFlipsAccepted += straightenResult.worseningFlipsAccepted;
      result.totalTimeMs += std::chrono::duration<double, std::milli>(end - start).count();
      if (straightenResult.iterations > result.maxIterationsOnPath) {
        result.maxIterationsOnPath = straightenResult.iterations;
      }
    }

    allResults.push_back(result);
  }

  // Print results table
  std::cout << "| Configuration | Avg Iters | Max Iters | Accepted | Rejected | Worsening | Improvement | Avg Time |" << std::endl;
  std::cout << "|---------------|-----------|-----------|----------|----------|-----------|-------------|----------|" << std::endl;

  for (const auto& r : allResults) {
    std::cout << "| " << r.name;
    std::cout << " | " << std::fixed << std::setprecision(1) << r.avgIterations();
    std::cout << " | " << r.maxIterationsOnPath;
    std::cout << " | " << r.totalFlipsAccepted;
    std::cout << " | " << r.totalFlipsRejected;
    std::cout << " | " << r.worseningFlipsAccepted;
    std::cout << " | " << std::fixed << std::setprecision(3) << r.avgImprovement() << "%";
    std::cout << " | " << std::fixed << std::setprecision(2) << r.avgTimeMs() << "ms";
    std::cout << " |" << std::endl;
  }

  std::cout << std::endl;
  std::cout << "## Analysis" << std::endl;
  std::cout << std::endl;
  std::cout << "- **Worsening flips**: How many flips were accepted that made the path longer" << std::endl;
  std::cout << "- **Improvement**: Final path length vs initial funnel length (before ANY straightening)" << std::endl;
  std::cout << "- If worsening flips help escape local minima, we'd see higher improvement despite accepting worse flips" << std::endl;

  return 0;
}
