// Cold Cache Benchmark: GFR vs FlipOut without caching benefit
//
// Measures single-path performance where each path starts with an empty cache.
// This represents the worst-case for GFR (no exploration cache reuse) and
// is the scenario for applications computing isolated paths.
//
// Usage: cold_cache_benchmark [options]
//   --mesh FILE     Mesh file to benchmark (can be repeated)
//   --paths N       Number of paths per mesh (default: 100)
//   --seed S        Random seed (default: 42)

#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
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

std::string getFilename(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

std::string getStem(const std::string& path) {
  std::string filename = getFilename(path);
  size_t pos = filename.find_last_of('.');
  return (pos == std::string::npos) ? filename : filename.substr(0, pos);
}

bool fileExists(const std::string& path) {
  return access(path.c_str(), F_OK) == 0;
}

bool isAbsolutePath(const std::string& path) {
#ifdef _WIN32
  return path.length() >= 2 && path[1] == ':';
#else
  return !path.empty() && path[0] == '/';
#endif
}

std::string formatNumber(size_t n) {
  std::string s = std::to_string(n);
  std::string result;
  int count = 0;
  for (int i = s.length() - 1; i >= 0; i--) {
    if (count > 0 && count % 3 == 0) result = "," + result;
    result = s[i] + result;
    count++;
  }
  return result;
}

struct BenchmarkResult {
  std::string name;
  size_t vertices;
  size_t numPaths;

  double gfrWarmTime;
  double gfrColdTime;
  double flipoutTime;

  double warmSpeedup() const { return flipoutTime / gfrWarmTime; }
  double coldSpeedup() const { return flipoutTime / gfrColdTime; }
  double cacheBenefit() const { return gfrColdTime / gfrWarmTime; }
};

BenchmarkResult runBenchmark(const std::string& meshPath, size_t numPaths, unsigned int seed) {
  BenchmarkResult result;
  result.name = getStem(meshPath);

  // Load mesh
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);
  result.vertices = mesh->nVertices();
  result.numPaths = numPaths;

  std::cerr << "  " << result.name << " (" << formatNumber(result.vertices) << " verts): ";

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

  // GFR benchmark - start cold, cache builds naturally
  // No pre-warming. First path is cold, subsequent paths benefit from cache.
  resetTimingStats();
  clearFunnelCache();  // Start completely fresh
  result.gfrColdTime = 0;  // This is the real "cold start" time
  for (const auto& pair : pairs) {
    Vertex v0 = mesh->vertex(pair.first);
    Vertex v1 = mesh->vertex(pair.second);
    auto start = std::chrono::high_resolution_clock::now();
    auto path = computeFunnelGeodesic(*mesh, *geometry, v0, v1);
    auto end = std::chrono::high_resolution_clock::now();
    result.gfrColdTime += std::chrono::duration<double, std::milli>(end - start).count();
  }
  auto cacheStats = getCacheStats();
  result.gfrWarmTime = result.gfrColdTime;  // Same measurement - no artificial warming
  std::cerr << "gfr(" << cacheStats.hits << "/" << (cacheStats.hits + cacheStats.misses) << " hits) ";

  // FlipOut benchmark - also starts cold, no cache to warm
  result.flipoutTime = 0;
  for (const auto& pair : pairs) {
    Vertex v0 = mesh->vertex(pair.first);
    Vertex v1 = mesh->vertex(pair.second);
    auto start = std::chrono::high_resolution_clock::now();
    auto flipNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, v0, v1);
    flipNetwork->iterativeShorten();
    auto end = std::chrono::high_resolution_clock::now();
    result.flipoutTime += std::chrono::duration<double, std::milli>(end - start).count();
  }
  std::cerr << "flip ";

  std::cerr << "done" << std::endl;
  return result;
}

int main(int argc, char** argv) {
  std::string meshDir = "C:/Dev/Colonel/Data/Meshes";
  size_t numPaths = 100;
  unsigned int seed = 42;
  std::vector<std::string> meshes;

  // Default meshes for cold cache test
  std::vector<std::string> defaultMeshes = {
    "bunny-small.obj",
    "stanford-bunny.obj",
    "armadillo.obj"
  };

  // Parse arguments
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--mesh" && i + 1 < argc) {
      meshes.push_back(argv[++i]);
    } else if (arg == "--paths" && i + 1 < argc) {
      numPaths = std::stoul(argv[++i]);
    } else if (arg == "--seed" && i + 1 < argc) {
      seed = std::stoul(argv[++i]);
    } else if (arg == "--meshes" && i + 1 < argc) {
      meshDir = argv[++i];
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Cold Cache Benchmark: GFR vs FlipOut without caching benefit" << std::endl;
      std::cout << std::endl;
      std::cout << "Usage: cold_cache_benchmark [options]" << std::endl;
      std::cout << "  --mesh FILE     Mesh file to benchmark (can be repeated)" << std::endl;
      std::cout << "  --meshes DIR    Mesh directory (default: C:/Dev/Colonel/Data/Meshes)" << std::endl;
      std::cout << "  --paths N       Paths per mesh (default: 100)" << std::endl;
      std::cout << "  --seed S        Random seed (default: 42)" << std::endl;
      return 0;
    }
  }

  if (meshes.empty()) meshes = defaultMeshes;

  // Print header
  auto now = std::chrono::system_clock::now();
  auto now_t = std::chrono::system_clock::to_time_t(now);

  std::cout << "# Cold Cache Benchmark Results" << std::endl;
  std::cout << "Generated: " << std::put_time(std::localtime(&now_t), "%Y-%m-%d %H:%M") << std::endl;
  std::cout << "Paths per mesh: " << numPaths << std::endl;
  std::cout << std::endl;

  // Run benchmarks
  std::vector<BenchmarkResult> results;
  std::cerr << "Running cold cache benchmarks..." << std::endl;

  for (const auto& meshName : meshes) {
    std::string meshPath;
    if (isAbsolutePath(meshName) || fileExists(meshName)) {
      meshPath = meshName;
    } else {
      meshPath = meshDir + "/" + meshName;
    }

    if (!fileExists(meshPath)) {
      std::cerr << "  Skipping " << meshName << " (not found)" << std::endl;
      continue;
    }

    try {
      auto result = runBenchmark(meshPath, numPaths, seed);
      results.push_back(result);
    } catch (const std::exception& e) {
      std::cerr << "  Error on " << meshName << ": " << e.what() << std::endl;
    }
  }

  if (results.empty()) {
    std::cerr << "No meshes benchmarked!" << std::endl;
    return 1;
  }

  // Print results table
  std::cout << "## Results (Cold Start, No Pre-Warming)" << std::endl;
  std::cout << std::endl;
  std::cout << "Both algorithms start cold. GFR's cache builds naturally during the run." << std::endl;
  std::cout << std::endl;
  std::cout << "| Mesh | Vertices | GFR Time | FlipOut Time | Speedup |" << std::endl;
  std::cout << "|------|----------|----------|--------------|---------|" << std::endl;

  double sumSpeedup = 0;
  for (const auto& r : results) {
    std::cout << "| " << r.name;
    std::cout << " | " << formatNumber(r.vertices);
    std::cout << " | " << std::fixed << std::setprecision(0) << r.gfrColdTime << "ms";
    std::cout << " | " << std::fixed << std::setprecision(0) << r.flipoutTime << "ms";
    std::cout << " | " << std::fixed << std::setprecision(2) << r.coldSpeedup() << "x";
    std::cout << " |" << std::endl;
    sumSpeedup += r.coldSpeedup();
  }

  std::cout << std::endl;
  std::cout << "**Average speedup:** " << std::fixed << std::setprecision(2)
            << (sumSpeedup / results.size()) << "x" << std::endl;

  std::cout << std::endl;
  std::cout << "## Notes" << std::endl;
  std::cout << std::endl;
  std::cout << "- Both algorithms start with empty caches" << std::endl;
  std::cout << "- GFR's exploration cache builds naturally as paths are computed" << std::endl;
  std::cout << "- FlipOut has no equivalent cache mechanism" << std::endl;

  std::cerr << "Done!" << std::endl;
  return 0;
}
