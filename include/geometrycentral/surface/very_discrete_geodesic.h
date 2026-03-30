#pragma once

// VeryDiscreteGeodesic - A* pathfinder with L5 multi-face exploration
// Produces better initial face strips than edge-only Dijkstra by unfolding
// up to 5 faces deep to find better discrete path approximations.

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <array>

namespace geometrycentral {
namespace surface {
namespace very_discrete_geodesic {

// Candidate vertices reachable via unfolding (L5 max depth)
enum class CandidateName : uint8_t {
  L1 = 0,
  L2L = 1, L3L = 2, L4L = 3, L5L = 4, L5LM = 5,
  L2R = 6, L3R = 7, L4R = 8, L5R = 9, L5RM = 10,
  None = 11
};
constexpr size_t NUM_CANDIDATES = 11;

struct CandidateVertex {
  CandidateName name = CandidateName::None;
  Vertex vertex;
  Vector2 flatPosition = {std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};
  double distance = std::numeric_limits<double>::quiet_NaN();
  bool isReachable = false;

  bool hasVertex() const { return vertex != Vertex(); }
};

struct ExplorationResult {
  Corner corner;
  std::array<CandidateVertex, NUM_CANDIDATES> candidates;

  CandidateVertex& operator[](CandidateName name) { return candidates[static_cast<size_t>(name)]; }
  const CandidateVertex& operator[](CandidateName name) const { return candidates[static_cast<size_t>(name)]; }
};

struct VertexPathStep {
  Vertex from;
  Vertex to;
  bool isApexJump = false;
  std::vector<Face> crossedFaces;
};

struct PathStep {
  Vertex from;
  Vertex to;
  bool isExplorerJump = false;
  CandidateName candidateName = CandidateName::None;
  double distance = 0.0;
  Corner sourceCorner;
};

struct PathResult {
  std::vector<Vertex> path;
  std::vector<PathStep> steps;
  bool isComplete = false;
  bool isFallback = false;
};

enum class WalkDirection {
  Clockwise,
  CounterClockwise
};

struct WalkResult {
  std::vector<Face> faces;
  Face finalFace;
  bool reachedTarget = false;
  std::string error = "";
  bool hasFinalFace() const { return finalFace != Face(); }
};

constexpr double HEURISTIC_WEIGHT = 1.0;

// ============================================================================
// Explorer
// ============================================================================

ExplorationResult explore(Corner corner, VertexPositionGeometry& geom);

// ============================================================================
// Geometry helpers
// ============================================================================

Vector2 flattenVertex(Vertex v, Halfedge portal, Vector2 flatA, Vector2 flatB,
                      VertexPositionGeometry& geom);

Vector2 getFlatPosition(Vertex v, Halfedge portal, Vector2 flatP1, Vector2 flatP2, Vector2 flatApex);

bool segmentCrossesPortal(Vector2 v0, Vector2 target, Vector2 a, Vector2 b, double eps);

Vector2 computeTriangleApex(Vector2 a, Vector2 b, double distA, double distB, bool pickPositiveY);

// ============================================================================
// Pathfinding (standalone, primarily for testing)
// ============================================================================

PathResult findPath(Vertex from, Vertex to,
                    ManifoldSurfaceMesh& mesh,
                    VertexPositionGeometry& geom);

std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> findGeodesicPath(
    Vertex from, Vertex to,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom);

std::pair<std::vector<Face>, std::vector<Vertex>> findFaceStripWithPath(
    Vertex from, Vertex to,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom);

double computePathDistance(const std::vector<PathStep>& steps);

std::vector<Face> getCrossedFaces(Corner corner, CandidateName candidateName);

// ============================================================================
// Face strip construction
// ============================================================================

WalkDirection determineWalkDirection(Vertex prev, Vertex current, Vertex next,
                                     VertexPositionGeometry& geom);

WalkResult walkToOutgoingEdge(Face startFace, Vertex vertex, Vertex targetVertex,
                               WalkDirection direction);

WalkResult walkToFace(Face startFace, Face targetFace, Vertex vertex, WalkDirection direction);

Face selectFirstFace(Vertex v0, Vertex v1, Vertex v2,
                     VertexPositionGeometry& geom,
                     Face targetFace = Face());

std::vector<Face> convertPath(const std::vector<Vertex>& path,
                               const std::vector<VertexPathStep>& steps,
                               Vertex entryVertex, Vertex exitVertex,
                               VertexPositionGeometry& geom);

// ============================================================================
// Face strip utilities
// ============================================================================

Face getSharedFace(Vertex v1, Vertex v2);
std::vector<Face> getAllSharedFaces(Vertex v1, Vertex v2);
bool areFacesAdjacent(Face f1, Face f2);
bool faceContainsEdge(Face face, Vertex v1, Vertex v2);
bool faceContainsVertex(Face face, Vertex vertex);
Halfedge findHalfedgeToVertex(Face face, Vertex vertex);
Halfedge findHalfedgeFromVertex(Face face, Vertex vertex);

// ============================================================================
// CachedVeryDiscreteGeodesicPathfinder
//
// Production pathfinder with aggressive caching:
// 1. Per-corner exploration results (L5 unfolding)
// 2. Per-vertex boundary edge counts
// 3. Per-vertex adjacency arrays
// 4. Reusable A* containers
// ============================================================================

class CachedVeryDiscreteGeodesicPathfinder {
public:
  CachedVeryDiscreteGeodesicPathfinder(ManifoldSurfaceMesh& mesh,
                                        VertexPositionGeometry& geom);

  std::pair<std::vector<Face>, std::vector<Vertex>> findFaceStripWithPath(
      Vertex from, Vertex to);

  std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> findGeodesicPath(
      Vertex from, Vertex to);

  PathResult findPath(Vertex from, Vertex to);

  void clearCache();
  void precomputeVertexData();

  size_t getCacheHits() const { return cacheHits; }
  size_t getCacheMisses() const { return cacheMisses; }
  size_t getExplorationCacheSize() const { return explorationCache.size(); }

private:
  ManifoldSurfaceMesh& mesh;
  VertexPositionGeometry& geom;

  std::unordered_map<size_t, ExplorationResult> explorationCache;

  VertexData<int> boundaryEdgeCount;
  VertexData<std::vector<Vertex>> adjacentVertices;
  VertexData<std::vector<Corner>> adjacentCorners;
  bool vertexDataComputed = false;

  std::unordered_map<size_t, double> gScore;
  std::unordered_map<size_t, std::pair<Vertex, PathStep>> cameFrom;
  std::unordered_set<size_t> closed;
  std::unordered_set<size_t> seenVertices;
  std::vector<std::tuple<Vertex, double, PathStep>> neighbors;

  size_t cacheHits = 0;
  size_t cacheMisses = 0;

  void ensureVertexData();
  const ExplorationResult& getOrComputeExploration(Corner corner);
  const std::vector<std::tuple<Vertex, double, PathStep>>& getNeighbors(
      Vertex current, Vertex goal, double currentG);
  std::tuple<std::vector<Vertex>, std::vector<PathStep>, bool> reconstructPath(
      Vertex from, Vertex to);
  std::vector<Face> getCrossedFaces(Corner corner, CandidateName name);
};

} // namespace very_discrete_geodesic
} // namespace surface
} // namespace geometrycentral
