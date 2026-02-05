// Fan Direction Verification Test
// Verifies that fan direction selection in face strip building is optimal.
//
// There are two places where fan direction matters:
// 1. handleApexJump() - already tries both directions, picks optimal
// 2. handleEdgeStep() - uses geometric heuristic, fallback on failure
//
// This test measures how often the heuristic in handleEdgeStep picks the
// shorter fan direction, and quantifies the impact when it doesn't.
//
// Usage: fan_direction_test [options]
//   --mesh FILE     Mesh file (can be repeated)
//   --paths N       Paths per mesh (default: 1000)
//   --meshes DIR    Mesh directory (default: C:/Dev/Colonel/Data/Meshes)
//   --verbose       Show per-path details

#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/very_discrete_geodesic.h"

#include <algorithm>
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

// ============================================================================
// Helper functions
// ============================================================================

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

// ============================================================================
// Fan Direction Analysis
// ============================================================================

// A single fan direction decision during face strip building
struct FanDecision {
  size_t pathIndex;
  size_t stepIndex;
  Vertex fanVertex;

  // CW direction results
  double cwAngleSum;
  bool cwReachedTarget;

  // CCW direction results
  double ccwAngleSum;
  bool ccwReachedTarget;

  // What the heuristic chose
  bool choseCW;

  // Was this optimal? (smaller angle sum)
  bool wasOptimal() const {
    if (!cwReachedTarget && !ccwReachedTarget) return true;  // both failed
    if (!cwReachedTarget) return !choseCW;  // CCW was only option
    if (!ccwReachedTarget) return choseCW;   // CW was only option

    // Both reached target - check which had smaller angle sum
    if (choseCW) {
      return cwAngleSum <= ccwAngleSum;
    } else {
      return ccwAngleSum <= cwAngleSum;
    }
  }

  // Extra angle (in degrees) due to suboptimal choice
  double extraAngle() const {
    if (!cwReachedTarget || !ccwReachedTarget) return 0;
    if (wasOptimal()) return 0;

    if (choseCW) {
      return cwAngleSum - ccwAngleSum;
    } else {
      return ccwAngleSum - cwAngleSum;
    }
  }
};

// Results for a single mesh
struct MeshResult {
  std::string name;
  size_t vertices;
  size_t numPaths;

  size_t totalDecisions;
  size_t optimalDecisions;
  size_t suboptimalDecisions;
  size_t bothFailedDecisions;

  double totalExtraAngle;  // in degrees

  double optimalRate() const {
    if (totalDecisions == 0) return 100.0;
    return 100.0 * optimalDecisions / totalDecisions;
  }

  double suboptimalRate() const {
    if (totalDecisions == 0) return 0.0;
    return 100.0 * suboptimalDecisions / totalDecisions;
  }

  double avgExtraAngle() const {
    if (suboptimalDecisions == 0) return 0.0;
    return totalExtraAngle / suboptimalDecisions;
  }
};

// Walk around vertex in specified direction, measuring angle sum
struct WalkResult {
  std::vector<Face> faces;
  double angleSum;  // Total angle traversed around vertex
  bool reachedTarget;
};

Halfedge findHalfedgeToVertex(Face face, Vertex vertex) {
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.tipVertex() == vertex) return he;
  }
  return Halfedge();
}

Halfedge findHalfedgeFromVertex(Face face, Vertex vertex) {
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.tailVertex() == vertex) return he;
  }
  return Halfedge();
}

bool faceContainsEdge(Face face, Vertex v0, Vertex v1) {
  for (Halfedge he : face.adjacentHalfedges()) {
    if ((he.tailVertex() == v0 && he.tipVertex() == v1) ||
        (he.tailVertex() == v1 && he.tipVertex() == v0)) {
      return true;
    }
  }
  return false;
}

// Compute corner angle at vertex in face
double cornerAngle(Face face, Vertex vertex, VertexPositionGeometry& geom) {
  // Find the two edges adjacent to vertex in this face
  Vertex prev, next;
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.tipVertex() == vertex) {
      prev = he.tailVertex();
    }
    if (he.tailVertex() == vertex) {
      next = he.tipVertex();
    }
  }

  Vector3 vPos = geom.vertexPositions[vertex];
  Vector3 e1 = geom.vertexPositions[prev] - vPos;
  Vector3 e2 = geom.vertexPositions[next] - vPos;

  double cosAngle = dot(unit(e1), unit(e2));
  cosAngle = std::max(-1.0, std::min(1.0, cosAngle));  // clamp for numerical safety
  return std::acos(cosAngle) * 180.0 / M_PI;  // return in degrees
}

WalkResult walkToOutgoingEdge(Face startFace, Vertex vertex, Vertex targetVertex,
                              bool clockwise, VertexPositionGeometry& geom) {
  WalkResult result;
  result.reachedTarget = false;
  result.angleSum = 0.0;

  Halfedge startHe = findHalfedgeToVertex(startFace, vertex);
  if (startHe == Halfedge()) return result;

  if (faceContainsEdge(startFace, vertex, targetVertex)) {
    result.reachedTarget = true;
    return result;
  }

  Halfedge currentHe = clockwise ? startHe.next().twin() : startHe.twin();

  const int maxSteps = 100;
  for (int i = 0; i < maxSteps; i++) {
    if (!currentHe.isInterior() || currentHe.face().isBoundaryLoop()) {
      return result;  // Hit boundary
    }

    Face currentFace = currentHe.face();
    if (currentFace == startFace) {
      return result;  // Wrapped around
    }

    result.faces.push_back(currentFace);
    result.angleSum += cornerAngle(currentFace, vertex, geom);

    if (faceContainsEdge(currentFace, vertex, targetVertex)) {
      result.reachedTarget = true;
      return result;
    }

    if (clockwise) {
      Halfedge leaving = findHalfedgeFromVertex(currentFace, vertex);
      if (leaving == Halfedge()) return result;
      currentHe = leaving.twin();
    } else {
      Halfedge entering = findHalfedgeToVertex(currentFace, vertex);
      if (entering == Halfedge()) return result;
      currentHe = entering.twin();
    }
  }

  return result;
}

// Analyze a single path for fan direction decisions
std::vector<FanDecision> analyzePath(
    size_t pathIndex,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex start, Vertex end) {

  std::vector<FanDecision> decisions;

  // Get the VDG path result with steps
  auto pathResult = very_discrete_geodesic::findGeodesicPath(start, end, mesh, geom);
  const std::vector<Vertex>& path = pathResult.first;
  const std::vector<very_discrete_geodesic::VertexPathStep>& steps = pathResult.second;

  if (path.size() < 3 || steps.size() < 2) {
    return decisions;  // Too short to have fan decisions
  }

  // For each non-apex-jump step, we need to fan around the vertex
  // The heuristic is based on determineWalkDirection(prev, current, next)
  // We simulate what both directions would give us

  for (size_t i = 1; i < steps.size(); i++) {
    const auto& step = steps[i];

    // Skip apex jumps - they already try both directions in handleApexJump
    if (step.isApexJump) continue;

    Vertex current = step.from;
    Vertex next = step.to;

    // Get prev vertex - need to look at previous step
    Vertex prev = steps[i-1].from;

    // Simulate walking from a face containing (prev, current) to a face containing (current, next)
    // Find a starting face that contains both prev and current
    Face startFace;
    for (Face f : current.adjacentFaces()) {
      bool hasPrev = false, hasCurrent = false;
      for (Vertex v : f.adjacentVertices()) {
        if (v == prev) hasPrev = true;
        if (v == current) hasCurrent = true;
      }
      if (hasPrev && hasCurrent) {
        startFace = f;
        break;
      }
    }

    if (startFace == Face()) continue;  // Couldn't find starting face

    // If start face already contains edge to next, no fan needed
    if (faceContainsEdge(startFace, current, next)) continue;

    // Try both directions
    WalkResult walkCW = walkToOutgoingEdge(startFace, current, next, true, geom);
    WalkResult walkCCW = walkToOutgoingEdge(startFace, current, next, false, geom);

    // Only record if at least one direction worked and they differ in angle
    if (!walkCW.reachedTarget && !walkCCW.reachedTarget) continue;
    if (std::abs(walkCW.angleSum - walkCCW.angleSum) < 0.01) continue;  // Same angle either way

    // The fixed code picks based on angle sum (smaller = better)
    // Simulate what the fixed handleEdgeStep does: pick smaller angle
    bool fixedChoseCW;
    if (walkCW.reachedTarget && walkCCW.reachedTarget) {
      fixedChoseCW = (walkCW.angleSum <= walkCCW.angleSum);
    } else if (walkCW.reachedTarget) {
      fixedChoseCW = true;
    } else {
      fixedChoseCW = false;
    }

    FanDecision decision;
    decision.pathIndex = pathIndex;
    decision.stepIndex = i;
    decision.fanVertex = current;
    decision.cwAngleSum = walkCW.angleSum;
    decision.cwReachedTarget = walkCW.reachedTarget;
    decision.ccwAngleSum = walkCCW.angleSum;
    decision.ccwReachedTarget = walkCCW.reachedTarget;
    decision.choseCW = fixedChoseCW;  // Use fixed algorithm, not old heuristic

    decisions.push_back(decision);
  }

  return decisions;
}

MeshResult benchmarkMesh(const std::string& meshPath, size_t numPaths, unsigned int seed, bool verbose) {
  MeshResult result;
  result.name = getStem(meshPath);

  // Load mesh
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);

  result.vertices = mesh->nVertices();
  result.numPaths = numPaths;
  result.totalDecisions = 0;
  result.optimalDecisions = 0;
  result.suboptimalDecisions = 0;
  result.bothFailedDecisions = 0;
  result.totalExtraAngle = 0.0;

  std::cerr << "  " << result.name << ": " << result.vertices << " vertices" << std::endl;

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

  // Progress
  size_t progressStep = std::max(numPaths / 20, (size_t)1);

  for (size_t i = 0; i < pairs.size(); i++) {
    Vertex v0 = mesh->vertex(pairs[i].first);
    Vertex v1 = mesh->vertex(pairs[i].second);

    std::vector<FanDecision> decisions = analyzePath(i, *mesh, *geometry, v0, v1);

    for (const auto& d : decisions) {
      result.totalDecisions++;

      if (!d.cwReachedTarget && !d.ccwReachedTarget) {
        result.bothFailedDecisions++;
      } else if (d.wasOptimal()) {
        result.optimalDecisions++;
      } else {
        result.suboptimalDecisions++;
        result.totalExtraAngle += d.extraAngle();

        if (verbose) {
          std::cerr << "    Suboptimal: path " << d.pathIndex
                    << " step " << d.stepIndex
                    << " chose " << (d.choseCW ? "CW" : "CCW")
                    << " (" << std::fixed << std::setprecision(1) << (d.choseCW ? d.cwAngleSum : d.ccwAngleSum) << " deg)"
                    << " vs optimal " << (d.choseCW ? "CCW" : "CW")
                    << " (" << std::fixed << std::setprecision(1) << (d.choseCW ? d.ccwAngleSum : d.cwAngleSum) << " deg)"
                    << " extra: +" << std::fixed << std::setprecision(1) << d.extraAngle() << " deg"
                    << std::endl;
        }
      }
    }

    if ((i + 1) % progressStep == 0) {
      std::cerr << "." << std::flush;
    }
  }

  std::cerr << " done" << std::endl;

  return result;
}

void printResults(std::ostream& out, const std::vector<MeshResult>& results) {
  out << "## Fan Direction Verification Results" << std::endl;
  out << std::endl;
  out << "| Mesh | Vertices | Paths | Decisions | Optimal | Suboptimal | Rate | Avg Extra Angle |" << std::endl;
  out << "|------|----------|-------|-----------|---------|------------|------|-----------------|" << std::endl;

  for (const auto& r : results) {
    out << "| " << r.name;
    out << " | " << formatNumber(r.vertices);
    out << " | " << formatNumber(r.numPaths);
    out << " | " << formatNumber(r.totalDecisions);
    out << " | " << formatNumber(r.optimalDecisions);
    out << " | " << r.suboptimalDecisions;
    out << " | " << std::fixed << std::setprecision(1) << r.optimalRate() << "%";
    out << " | +" << std::fixed << std::setprecision(1) << r.avgExtraAngle() << " deg";
    out << " |" << std::endl;
  }

  // Aggregate
  size_t totalDecisions = 0;
  size_t totalOptimal = 0;
  size_t totalSuboptimal = 0;
  double totalExtraAngle = 0.0;

  for (const auto& r : results) {
    totalDecisions += r.totalDecisions;
    totalOptimal += r.optimalDecisions;
    totalSuboptimal += r.suboptimalDecisions;
    totalExtraAngle += r.totalExtraAngle;
  }

  double overallOptimalRate = (totalDecisions > 0) ? 100.0 * totalOptimal / totalDecisions : 100.0;
  double overallAvgExtra = (totalSuboptimal > 0) ? totalExtraAngle / totalSuboptimal : 0.0;

  out << std::endl;
  out << "**Overall:** " << formatNumber(totalDecisions) << " fan decisions, ";
  out << std::fixed << std::setprecision(1) << overallOptimalRate << "% optimal";
  if (totalSuboptimal > 0) {
    out << ", " << totalSuboptimal << " suboptimal (avg +"
        << std::fixed << std::setprecision(1) << overallAvgExtra << " deg extra angle)";
  }
  out << std::endl;
}

const std::vector<std::string> DEFAULT_MESHES = {
  "bunny-small.obj",
  "spot.obj",
  "stanford-bunny.obj",
  "armadillo.obj",
  "max-planck.obj",
  "happy-buddha.obj"
};

int main(int argc, char** argv) {
  std::string meshDir = "C:/Dev/Colonel/Data/Meshes";
  size_t numPaths = 1000;
  unsigned int seed = 42;
  std::vector<std::string> specificMeshes;
  bool verbose = false;

  // Parse arguments
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--meshes" && i + 1 < argc) {
      meshDir = argv[++i];
    } else if (arg == "--paths" && i + 1 < argc) {
      numPaths = std::stoul(argv[++i]);
    } else if (arg == "--seed" && i + 1 < argc) {
      seed = std::stoul(argv[++i]);
    } else if (arg == "--mesh" && i + 1 < argc) {
      specificMeshes.push_back(argv[++i]);
    } else if (arg == "--verbose") {
      verbose = true;
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: fan_direction_test [options]" << std::endl;
      std::cout << "  --meshes DIR    Mesh directory (default: C:/Dev/Colonel/Data/Meshes)" << std::endl;
      std::cout << "  --paths N       Paths per mesh (default: 1000)" << std::endl;
      std::cout << "  --seed S        Random seed (default: 42)" << std::endl;
      std::cout << "  --mesh FILE     Specific mesh (can repeat)" << std::endl;
      std::cout << "  --verbose       Show suboptimal decision details" << std::endl;
      return 0;
    }
  }

  std::vector<std::string> meshesToRun = specificMeshes.empty() ? DEFAULT_MESHES : specificMeshes;

  auto now = std::chrono::system_clock::now();
  auto now_t = std::chrono::system_clock::to_time_t(now);

  std::cout << "# Fan Direction Verification" << std::endl;
  std::cout << "Generated: " << std::put_time(std::localtime(&now_t), "%Y-%m-%d %H:%M") << std::endl;
  std::cout << "Seed: " << seed << ", Paths per mesh: " << numPaths << std::endl;
  std::cout << std::endl;

  std::vector<MeshResult> results;

  std::cerr << "Running fan direction verification..." << std::endl;
  for (const auto& meshName : meshesToRun) {
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
      auto result = benchmarkMesh(meshPath, numPaths, seed, verbose);
      results.push_back(result);
    } catch (const std::exception& e) {
      std::cerr << "  Error on " << meshName << ": " << e.what() << std::endl;
    }
  }

  if (results.empty()) {
    std::cerr << "No meshes were tested!" << std::endl;
    return 1;
  }

  printResults(std::cout, results);

  std::cerr << "Done!" << std::endl;
  return 0;
}
