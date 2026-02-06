// Multi-Mesh GFR vs FlipOut Benchmark
// Produces paper-ready markdown tables comparing GFR and FlipOut across multiple meshes.
//
// Usage: multi_mesh_benchmark [options]
//   --meshes DIR    Directory containing mesh files (default: C:/Dev/Colonel/Data/Meshes)
//   --paths N       Number of random vertex pairs per mesh (default: 1000)
//   --seed S        Random seed for reproducibility (default: 42)
//   --output FILE   Output file (default: stdout)
//   --mesh FILE     Run on specific mesh only (can be repeated)
//   --verbose       Show per-mesh detailed statistics

#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <fstream>
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

// Simple path utilities (C++14 compatible)
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

// Mesh result for a single mesh
struct MeshResult {
  std::string name;
  size_t vertices;
  size_t faces;
  size_t numPaths;

  double gfrTotalTime;
  double flipoutTotalTime;
  double gfrTotalLength;
  double flipoutTotalLength;

  size_t gfrWins;
  size_t flipoutWins;
  size_t ties;

  size_t gfrTotalIters;
  size_t flipoutTotalIters;

  // Investigation stats
  size_t gfrTotalFaces;           // Total sleeve face count
  double gfrTotalInitialLength;   // Total initial funnel length (before straightening)
  size_t cacheHits;
  size_t cacheMisses;
  double aStarMs;
  double flattenMs;
  double straightenMs;

  double speedup() const { return flipoutTotalTime / gfrTotalTime; }
  double vsFlipout() const { return 100.0 * (gfrTotalLength - flipoutTotalLength) / flipoutTotalLength; }
  double avgGfrTime() const { return gfrTotalTime / numPaths; }
  double avgFlipoutTime() const { return flipoutTotalTime / numPaths; }
  double avgFaces() const { return (double)gfrTotalFaces / numPaths; }
  double avgIters() const { return (double)gfrTotalIters / numPaths; }
  double itersPerFace() const { return (double)gfrTotalIters / gfrTotalFaces; }
  double straighteningImprovement() const {
    return 100.0 * (gfrTotalInitialLength - gfrTotalLength) / gfrTotalInitialLength;
  }
  double cacheHitRate() const {
    return 100.0 * cacheHits / (cacheHits + cacheMisses);
  }
};

// Default meshes to benchmark (matching paper-draft.md)
const std::vector<std::string> DEFAULT_MESHES = {
  "bunny-small.obj",  // ~2.5K vertices
  "spot.obj",         // ~2.9K vertices
  "pig.obj",          // ~1.8K vertices
  "stanford-bunny.obj", // ~35K vertices
  "armadillo.obj",    // ~50K vertices
  "max-planck.obj",   // ~50K vertices
  "happy-buddha.obj"  // ~49K vertices
};

std::string formatTime(double ms) {
  std::ostringstream ss;
  if (ms < 1000) {
    ss << std::fixed << std::setprecision(0) << ms << "ms";
  } else {
    ss << std::fixed << std::setprecision(1) << (ms / 1000.0) << "s";
  }
  return ss.str();
}

std::string formatNumber(size_t n) {
  // Format with comma thousands separator
  std::string s = std::to_string(n);
  std::string result;
  int count = 0;
  for (int i = s.length() - 1; i >= 0; i--) {
    if (count > 0 && count % 3 == 0) {
      result = "," + result;
    }
    result = s[i] + result;
    count++;
  }
  return result;
}

MeshResult benchmarkMesh(const std::string& meshPath, size_t numPaths, unsigned int seed, bool verbose) {
  MeshResult result;

  // Extract mesh name from path
  result.name = getStem(meshPath);

  // Load mesh
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);

  result.vertices = mesh->nVertices();
  result.faces = mesh->nFaces();
  result.numPaths = numPaths;

  if (verbose) {
    std::cerr << "  " << result.name << ": " << result.vertices << " vertices, "
              << result.faces << " faces" << std::endl;
  }

  // Generate random vertex pairs
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> dist(0, mesh->nVertices() - 1);

  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.reserve(numPaths);
  while (pairs.size() < numPaths) {
    size_t v0 = dist(rng);
    size_t v1 = dist(rng);
    if (v0 != v1) {
      pairs.emplace_back(v0, v1);
    }
  }

  // Reset timing stats and clear cache - start cold, no pre-warming
  resetTimingStats();
  clearFunnelCache();

  // Initialize accumulators
  result.gfrTotalTime = 0;
  result.flipoutTotalTime = 0;
  result.gfrTotalLength = 0;
  result.flipoutTotalLength = 0;
  result.gfrWins = 0;
  result.flipoutWins = 0;
  result.ties = 0;
  result.gfrTotalIters = 0;
  result.flipoutTotalIters = 0;
  result.gfrTotalFaces = 0;
  result.gfrTotalInitialLength = 0;

  // Progress indicator
  size_t progressStep = std::max(numPaths / 20, (size_t)1);

  for (size_t i = 0; i < pairs.size(); i++) {
    Vertex v0 = mesh->vertex(pairs[i].first);
    Vertex v1 = mesh->vertex(pairs[i].second);

    // Run GFR
    auto gfrStart = std::chrono::high_resolution_clock::now();
    auto gfrPath = computeFunnelGeodesic(*mesh, *geometry, v0, v1);
    auto gfrEnd = std::chrono::high_resolution_clock::now();
    double gfrTime = std::chrono::duration<double, std::milli>(gfrEnd - gfrStart).count();
    double gfrLen = gfrPath->length();
    size_t gfrIters = gfrPath->iterationCount();
    size_t gfrFaces = gfrPath->faceCount();
    double gfrInitialLen = gfrPath->initialFunnelLength();

    // Run FlipOut
    auto flipStart = std::chrono::high_resolution_clock::now();
    auto flipNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, v0, v1);
    flipNetwork->iterativeShorten();
    auto flipEnd = std::chrono::high_resolution_clock::now();
    double flipTime = std::chrono::duration<double, std::milli>(flipEnd - flipStart).count();
    double flipLen = flipNetwork->length();
    size_t flipIters = flipNetwork->nShortenIters;

    // Accumulate
    result.gfrTotalTime += gfrTime;
    result.flipoutTotalTime += flipTime;
    result.gfrTotalLength += gfrLen;
    result.flipoutTotalLength += flipLen;
    result.gfrTotalIters += gfrIters;
    result.flipoutTotalIters += flipIters;
    result.gfrTotalFaces += gfrFaces;
    result.gfrTotalInitialLength += gfrInitialLen;

    // Compare lengths (within tolerance)
    const double tol = 1e-8;
    if (gfrLen < flipLen - tol) {
      result.gfrWins++;
    } else if (flipLen < gfrLen - tol) {
      result.flipoutWins++;
    } else {
      result.ties++;
    }

    // Progress
    if (verbose && (i + 1) % progressStep == 0) {
      std::cerr << "." << std::flush;
    }
  }

  if (verbose) {
    std::cerr << " done" << std::endl;
  }

  // Capture cache and timing stats
  auto cacheStats = getCacheStats();
  auto timingStats = getTimingStats();
  result.cacheHits = cacheStats.hits;
  result.cacheMisses = cacheStats.misses;
  result.aStarMs = timingStats.aStarMs;
  result.flattenMs = timingStats.flattenMs;
  result.straightenMs = timingStats.straightenMs;

  return result;
}

void printMarkdownTable(std::ostream& out, const std::vector<MeshResult>& results) {
  // Header
  out << "| Mesh | Vertices | Paths | GFR Time | FlipOut Time | Speedup | vs FlipOut |" << std::endl;
  out << "|------|----------|-------|----------|--------------|---------|------------|" << std::endl;

  // Rows
  for (const auto& r : results) {
    out << "| " << r.name;
    out << " | " << formatNumber(r.vertices);
    out << " | " << formatNumber(r.numPaths);
    out << " | " << formatTime(r.gfrTotalTime);
    out << " | " << formatTime(r.flipoutTotalTime);
    out << " | **" << std::fixed << std::setprecision(1) << r.speedup() << "x**";
    out << " | **" << std::fixed << std::setprecision(2) << r.vsFlipout() << "%**";
    out << " |" << std::endl;
  }
}

void printAggregate(std::ostream& out, const std::vector<MeshResult>& results) {
  size_t totalPaths = 0;
  double sumSpeedup = 0;
  double sumVsFlipout = 0;

  for (const auto& r : results) {
    totalPaths += r.numPaths;
    sumSpeedup += r.speedup();
    sumVsFlipout += r.vsFlipout();
  }

  // Simple average across meshes (each mesh weighted equally)
  double avgSpeedup = sumSpeedup / results.size();
  double avgVsFlipout = sumVsFlipout / results.size();

  out << std::endl;
  out << "**Average across " << results.size() << " meshes (" << formatNumber(totalPaths) << " paths):** ";
  out << "GFR is " << std::fixed << std::setprecision(2) << avgSpeedup << "x faster ";
  out << "and produces " << std::fixed << std::setprecision(2) << std::abs(avgVsFlipout) << "% ";
  out << (avgVsFlipout < 0 ? "shorter" : "longer") << " paths." << std::endl;
}

void printDiagnosticTable(std::ostream& out, const std::vector<MeshResult>& results) {
  out << std::endl;
  out << "## Investigation: Speedup vs Mesh Characteristics" << std::endl;
  out << std::endl;
  out << "| Mesh | Vertices | Speedup | AvgFaces | AvgIters | Iters/Face | Cache% | Straighten% |" << std::endl;
  out << "|------|----------|---------|----------|----------|------------|--------|-------------|" << std::endl;

  for (const auto& r : results) {
    out << "| " << r.name;
    out << " | " << formatNumber(r.vertices);
    out << " | " << std::fixed << std::setprecision(2) << r.speedup() << "x";
    out << " | " << std::fixed << std::setprecision(1) << r.avgFaces();
    out << " | " << std::fixed << std::setprecision(1) << r.avgIters();
    out << " | " << std::fixed << std::setprecision(3) << r.itersPerFace();
    out << " | " << std::fixed << std::setprecision(1) << r.cacheHitRate() << "%";
    out << " | " << std::fixed << std::setprecision(2) << r.straighteningImprovement() << "%";
    out << " |" << std::endl;
  }

  // Analysis summary
  out << std::endl;
  out << "**Analysis:**" << std::endl;

  // Find lowest and highest speedup meshes
  auto minIt = std::min_element(results.begin(), results.end(),
      [](const MeshResult& a, const MeshResult& b) { return a.speedup() < b.speedup(); });
  auto maxIt = std::max_element(results.begin(), results.end(),
      [](const MeshResult& a, const MeshResult& b) { return a.speedup() < b.speedup(); });

  out << "- **Lowest speedup:** " << minIt->name << " (" << std::fixed << std::setprecision(2) << minIt->speedup() << "x)" << std::endl;
  out << "  - Avg faces: " << std::setprecision(1) << minIt->avgFaces()
      << ", Iters/face: " << std::setprecision(3) << minIt->itersPerFace()
      << ", Cache: " << std::setprecision(1) << minIt->cacheHitRate() << "%" << std::endl;
  out << "- **Highest speedup:** " << maxIt->name << " (" << std::fixed << std::setprecision(2) << maxIt->speedup() << "x)" << std::endl;
  out << "  - Avg faces: " << std::setprecision(1) << maxIt->avgFaces()
      << ", Iters/face: " << std::setprecision(3) << maxIt->itersPerFace()
      << ", Cache: " << std::setprecision(1) << maxIt->cacheHitRate() << "%" << std::endl;
}

void printVerboseStats(std::ostream& out, const MeshResult& r) {
  out << std::endl;
  out << "### " << r.name << " (" << formatNumber(r.vertices) << " vertices, "
      << formatNumber(r.numPaths) << " paths)" << std::endl;
  out << "- GFR wins: " << r.gfrWins << " (" << std::fixed << std::setprecision(1)
      << (100.0 * r.gfrWins / r.numPaths) << "%)" << std::endl;
  out << "- FlipOut wins: " << r.flipoutWins << " (" << std::fixed << std::setprecision(1)
      << (100.0 * r.flipoutWins / r.numPaths) << "%)" << std::endl;
  out << "- Ties: " << r.ties << " (" << std::fixed << std::setprecision(1)
      << (100.0 * r.ties / r.numPaths) << "%)" << std::endl;
  out << "- GFR avg iterations: " << std::fixed << std::setprecision(1)
      << r.avgIters() << std::endl;
  out << "- FlipOut avg iterations: " << std::fixed << std::setprecision(1)
      << (double)r.flipoutTotalIters / r.numPaths << std::endl;
  out << std::endl;
  out << "**Investigation Stats:**" << std::endl;
  out << "- Avg sleeve faces: " << std::fixed << std::setprecision(1) << r.avgFaces() << std::endl;
  out << "- Iters/face ratio: " << std::fixed << std::setprecision(3) << r.itersPerFace() << std::endl;
  out << "- Straightening improvement: " << std::fixed << std::setprecision(2)
      << r.straighteningImprovement() << "%" << std::endl;
  out << "- Cache hit rate: " << std::fixed << std::setprecision(1) << r.cacheHitRate()
      << "% (" << formatNumber(r.cacheHits) << " hits, " << formatNumber(r.cacheMisses) << " misses)" << std::endl;
  out << "- Timing: A*=" << std::fixed << std::setprecision(1) << r.aStarMs << "ms, "
      << "Flatten=" << r.flattenMs << "ms, "
      << "Straighten=" << r.straightenMs << "ms" << std::endl;
}

int main(int argc, char** argv) {
  // Default arguments
  std::string meshDir = "C:/Dev/Colonel/Data/Meshes";
  size_t numPaths = 1000;
  unsigned int seed = 42;
  std::string outputFile;
  std::vector<std::string> specificMeshes;
  bool verbose = false;
  bool diagnostic = false;

  // Parse arguments
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--meshes" && i + 1 < argc) {
      meshDir = argv[++i];
    } else if (arg == "--paths" && i + 1 < argc) {
      numPaths = std::stoul(argv[++i]);
    } else if (arg == "--seed" && i + 1 < argc) {
      seed = std::stoul(argv[++i]);
    } else if (arg == "--output" && i + 1 < argc) {
      outputFile = argv[++i];
    } else if (arg == "--mesh" && i + 1 < argc) {
      specificMeshes.push_back(argv[++i]);
    } else if (arg == "--verbose") {
      verbose = true;
    } else if (arg == "--diagnostic") {
      diagnostic = true;
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: multi_mesh_benchmark [options]" << std::endl;
      std::cout << "  --meshes DIR    Mesh directory (default: C:/Dev/Colonel/Data/Meshes)" << std::endl;
      std::cout << "  --paths N       Paths per mesh (default: 1000)" << std::endl;
      std::cout << "  --seed S        Random seed (default: 42)" << std::endl;
      std::cout << "  --output FILE   Output file (default: stdout)" << std::endl;
      std::cout << "  --mesh FILE     Specific mesh (can repeat)" << std::endl;
      std::cout << "  --verbose       Show detailed stats per mesh" << std::endl;
      std::cout << "  --diagnostic    Show speedup investigation table" << std::endl;
      return 0;
    }
  }

  // Determine which meshes to run
  std::vector<std::string> meshesToRun;
  if (!specificMeshes.empty()) {
    meshesToRun = specificMeshes;
  } else {
    meshesToRun = DEFAULT_MESHES;
  }

  // Open output stream
  std::ofstream fileOut;
  std::ostream* out = &std::cout;
  if (!outputFile.empty()) {
    fileOut.open(outputFile);
    out = &fileOut;
  }

  // Print header
  auto now = std::chrono::system_clock::now();
  auto now_t = std::chrono::system_clock::to_time_t(now);

  *out << "# GFR vs FlipOut Benchmark Results" << std::endl;
  *out << "Generated: " << std::put_time(std::localtime(&now_t), "%Y-%m-%d %H:%M") << std::endl;
  *out << "Seed: " << seed << std::endl;
  *out << std::endl;

  // Run benchmarks
  std::vector<MeshResult> results;

  std::cerr << "Running benchmarks..." << std::endl;
  for (const auto& meshName : meshesToRun) {
    std::string meshPath;

    // Check if it's an absolute path or just a filename
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
      auto result = benchmarkMesh(meshPath, numPaths, seed, verbose);
      results.push_back(result);
    } catch (const std::exception& e) {
      std::cerr << "  Error on " << meshName << ": " << e.what() << std::endl;
    }
  }

  if (results.empty()) {
    std::cerr << "No meshes were benchmarked!" << std::endl;
    return 1;
  }

  // Print results
  *out << "## Executive Summary" << std::endl;
  *out << std::endl;
  printMarkdownTable(*out, results);
  printAggregate(*out, results);

  if (diagnostic) {
    printDiagnosticTable(*out, results);
  }

  if (verbose) {
    *out << std::endl;
    *out << "## Detailed Statistics" << std::endl;
    for (const auto& r : results) {
      printVerboseStats(*out, r);
    }
  }

  std::cerr << "Done!" << std::endl;
  return 0;
}
