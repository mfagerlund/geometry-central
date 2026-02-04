// Initial Distance Invariant Test
// Verifies that GFR's initial funnel distance is NEVER longer than FlipOut's Dijkstra distance.
//
// The motivation: The sleeve should contain the shortest edge path as either
// its left or right boundary. The funnel algorithm finds the shortest path
// through the sleeve, which should be at most as long as any path along the
// sleeve boundaries (including the Dijkstra path).
//
// If this test fails, it indicates a bug in either:
// 1. Face strip construction (not including all necessary faces)
// 2. Sleeve flattening (distorting distances)
// 3. Portal/funnel algorithm (not finding true shortest path through sleeve)
//
// Usage: initial_distance_test [options]
//   --mesh FILE     Test specific mesh (can be repeated)
//   --paths N       Number of random paths per mesh (default: 1000)
//   --seed S        Random seed (default: 42)
//   --verbose       Show details for each failure
//   --all           Run on all available meshes

#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/very_discrete_geodesic.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _WIN32
#include <io.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif

using namespace geometrycentral;
using namespace geometrycentral::surface;

bool fileExists(const std::string& path) {
  return access(path.c_str(), F_OK) == 0;
}

std::string getStem(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  std::string filename = (pos == std::string::npos) ? path : path.substr(pos + 1);
  pos = filename.find_last_of('.');
  return (pos == std::string::npos) ? filename : filename.substr(0, pos);
}

struct FailureInfo {
  size_t v0, v1;
  double gfrInitial;
  double dijkstra;
  double difference;
  double percentDiff;
};

struct MeshTestResult {
  std::string name;
  size_t numPaths;
  size_t numFailures;
  double worstDiff;      // Worst case where GFR initial > Dijkstra (positive = bad)
  double avgDiff;        // Average (GFR initial - Dijkstra), negative = GFR shorter
  std::vector<FailureInfo> failures;
};

MeshTestResult testMesh(
    const std::string& meshPath,
    size_t numPaths,
    unsigned int seed,
    bool verbose) {

  MeshTestResult result;
  result.name = getStem(meshPath);
  result.numPaths = numPaths;
  result.numFailures = 0;
  result.worstDiff = 0;
  result.avgDiff = 0;

  // Load mesh
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  try {
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);
  } catch (const std::exception& e) {
    std::cerr << "Error loading " << meshPath << ": " << e.what() << std::endl;
    return result;
  }

  std::cout << "Testing " << result.name << " (" << mesh->nVertices() << " vertices, "
            << numPaths << " paths)... " << std::flush;

  // Generate random vertex pairs
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> dist(0, mesh->nVertices() - 1);

  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.reserve(numPaths);
  for (size_t i = 0; i < numPaths; i++) {
    size_t v0 = dist(rng);
    size_t v1 = dist(rng);
    if (v0 != v1) {
      pairs.emplace_back(v0, v1);
    } else {
      i--; // retry
    }
  }

  double sumDiff = 0;
  const double tol = 1e-9;  // Tolerance for floating point comparison

  for (const auto& pair : pairs) {
    Vertex v0 = mesh->vertex(pair.first);
    Vertex v1 = mesh->vertex(pair.second);

    // Compute GFR path (initial funnel distance is stored)
    auto gfrPath = computeFunnelGeodesic(*mesh, *geometry, v0, v1);
    double gfrInitial = gfrPath->initialFunnelLength();

    // Compute FlipOut Dijkstra path (length before iterativeShorten)
    auto flipNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, v0, v1);
    double dijkstra = flipNetwork->length();


    double diff = gfrInitial - dijkstra;
    sumDiff += diff;

    // Check invariant: GFR initial should be <= Dijkstra (within tolerance)
    if (gfrInitial > dijkstra + tol) {
      result.numFailures++;
      double pctDiff = 100.0 * diff / dijkstra;

      if (diff > result.worstDiff) {
        result.worstDiff = diff;
      }

      result.failures.push_back({
        pair.first, pair.second,
        gfrInitial, dijkstra,
        diff, pctDiff
      });
    }
  }

  result.avgDiff = sumDiff / numPaths;

  // Report
  if (result.numFailures == 0) {
    std::cout << "PASS (avg diff: " << std::fixed << std::setprecision(6)
              << result.avgDiff << ")" << std::endl;
  } else {
    std::cout << "FAIL - " << result.numFailures << " violations!" << std::endl;
    std::cout << "  Worst: GFR initial " << std::setprecision(4) << result.worstDiff
              << " longer than Dijkstra" << std::endl;

    if (verbose) {
      // Sort failures by severity
      std::sort(result.failures.begin(), result.failures.end(),
                [](const FailureInfo& a, const FailureInfo& b) {
                  return a.difference > b.difference;
                });

      std::cout << "  Top 10 violations:" << std::endl;
      for (size_t i = 0; i < std::min(size_t(10), result.failures.size()); i++) {
        const auto& f = result.failures[i];
        std::cout << "    V" << f.v0 << " -> V" << f.v1
                  << ": GFR=" << std::setprecision(6) << f.gfrInitial
                  << ", Dijkstra=" << f.dijkstra
                  << " (+" << std::setprecision(4) << f.percentDiff << "%)"
                  << std::endl;
      }
    }
  }

  return result;
}

int main(int argc, char** argv) {
  // Parse arguments
  std::vector<std::string> meshFiles;
  size_t numPaths = 1000;
  unsigned int seed = 42;
  bool verbose = false;
  bool useAll = false;
  std::string meshDir = "C:/Dev/Colonel/Data/Meshes";

  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--mesh" && i + 1 < argc) {
      meshFiles.push_back(argv[++i]);
    } else if (arg == "--paths" && i + 1 < argc) {
      numPaths = std::stoul(argv[++i]);
    } else if (arg == "--seed" && i + 1 < argc) {
      seed = std::stoul(argv[++i]);
    } else if (arg == "--verbose") {
      verbose = true;
    } else if (arg == "--all") {
      useAll = true;
    } else if (arg == "--meshdir" && i + 1 < argc) {
      meshDir = argv[++i];
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
      std::cout << "  --mesh FILE     Test specific mesh (can be repeated)" << std::endl;
      std::cout << "  --paths N       Number of random paths per mesh (default: 1000)" << std::endl;
      std::cout << "  --seed S        Random seed (default: 42)" << std::endl;
      std::cout << "  --verbose       Show details for each failure" << std::endl;
      std::cout << "  --all           Run on all available meshes" << std::endl;
      std::cout << "  --meshdir DIR   Directory for mesh files" << std::endl;
      return 0;
    }
  }

  // Default meshes if none specified
  if (meshFiles.empty() && !useAll) {
    meshFiles = {
      "bunny-small.obj",
      "spot.obj",
      "sphere320.obj",
      "pig.obj",
      "stanford-bunny.obj",
      "armadillo.obj",
      "max-planck.obj"
    };
  }

  // If --all, find all .obj files
  if (useAll) {
    meshFiles = {
      "bunny-small.obj", "spot.obj", "sphere320.obj", "pig.obj",
      "stanford-bunny.obj", "armadillo.obj", "max-planck.obj",
      "happy-buddha.obj", "cow.obj", "suzanne.obj", "teapot.obj",
      "torus.obj", "alligator.obj", "cheburashka.obj", "homer.obj",
      "fandisk.obj", "rocker-arm.obj"
    };
  }

  std::cout << "============================================" << std::endl;
  std::cout << "Initial Distance Invariant Test" << std::endl;
  std::cout << "============================================" << std::endl;
  std::cout << "Verifying: GFR initial funnel <= FlipOut Dijkstra" << std::endl;
  std::cout << "Paths per mesh: " << numPaths << ", seed: " << seed << std::endl;
  std::cout << std::endl;

  std::vector<MeshTestResult> results;
  size_t totalFailures = 0;

  for (const auto& meshFile : meshFiles) {
    std::string fullPath = meshFile;
    if (fullPath.find('/') == std::string::npos && fullPath.find('\\') == std::string::npos) {
      fullPath = meshDir + "/" + meshFile;
    }

    if (!fileExists(fullPath)) {
      std::cerr << "Skipping " << meshFile << " (not found)" << std::endl;
      continue;
    }

    auto result = testMesh(fullPath, numPaths, seed, verbose);
    results.push_back(result);
    totalFailures += result.numFailures;
  }

  // Summary
  std::cout << std::endl;
  std::cout << "============================================" << std::endl;
  std::cout << "SUMMARY" << std::endl;
  std::cout << "============================================" << std::endl;

  std::cout << std::left << std::setw(20) << "Mesh"
            << std::right << std::setw(8) << "Paths"
            << std::setw(10) << "Failures"
            << std::setw(12) << "Avg Diff"
            << std::setw(12) << "Worst" << std::endl;
  std::cout << std::string(62, '-') << std::endl;

  for (const auto& r : results) {
    std::cout << std::left << std::setw(20) << r.name
              << std::right << std::setw(8) << r.numPaths
              << std::setw(10) << r.numFailures
              << std::setw(12) << std::fixed << std::setprecision(6) << r.avgDiff
              << std::setw(12) << r.worstDiff << std::endl;
  }

  std::cout << std::endl;
  if (totalFailures == 0) {
    std::cout << "OVERALL: PASS - No invariant violations found!" << std::endl;
    std::cout << "GFR's initial funnel is always <= FlipOut's Dijkstra path." << std::endl;
    return 0;
  } else {
    std::cout << "OVERALL: FAIL - " << totalFailures << " invariant violations found!" << std::endl;
    std::cout << "This indicates a potential bug in face strip construction or funnel algorithm." << std::endl;
    return 1;
  }
}
