#pragma once

// VeryDiscreteGeodesicPathfinder - Line-by-line port from C#
//
// Original C# files:
//   - VeryDiscreteGeodesicPathfinder.cs
//   - VeryDiscreteGeodesicExplorer.cs
//   - VeryDiscreteGeodesicExplorerHelper.cs
//   - FaceStripWalker.cs
//   - VertexPathToFaceStripConverter.cs
//   - GeodesicPathJoiner.cs
//   - FaceStripUtils.cs
//
// This is an A* pathfinder using multi-face jumps via unfolding.
// Produces better initial approximations than edge-only Dijkstra.

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

// GFR uses fixed L5 exploration depth (no configuration needed)

// Candidate vertices reachable via unfolding (L5 max depth)
enum class CandidateName {
  L1,
  L2L, L3L, L4L, L5L,
  L2R, L3R, L4R, L5R,
  L5LM, L5RM,
  None
};

// ============================================================================
// C#: public enum BlockedReason
// ============================================================================
enum class BlockedReason {
  None_,  // Can't use "None" because it conflicts
  NotComputed,
  F1Boundary,
  V0FlatteningFailed,
  ParentFaceBoundary,
  FlatteningFailed,
  TargetFaceBoundary,
  PortalBlocked
};

// ============================================================================
// C#: public record CandidateVertex(...)
// ============================================================================
struct CandidateVertex {
  // C#: CandidateName Name,
  CandidateName name = CandidateName::None;

  // C#: Hedge.Vertex? Vertex,
  Vertex vertex;  // Invalid vertex means null

  // C#: Vector2 FlatPosition,
  Vector2 flatPosition = {std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};

  // C#: float Distance,
  double distance = std::numeric_limits<double>::quiet_NaN();

  // C#: bool IsReachable,
  bool isReachable = false;

  // C#: BlockedReason Blocked,
  BlockedReason blocked = BlockedReason::NotComputed;

  // C#: string? BlockedPortal = null
  std::string blockedPortal = "";

  // Helper to check if vertex is valid (not null in C# terms)
  bool hasVertex() const { return vertex != Vertex(); }
};

// Result of exploring from a corner (L5 depth = 11 candidates)
struct ExplorationResult {
  Corner corner;

  // L5 candidates: L1 + 4 left + 4 right + 2 middle
  CandidateVertex L1;
  CandidateVertex L2L, L3L, L4L, L5L, L5LM;
  CandidateVertex L2R, L3R, L4R, L5R, L5RM;

  std::vector<CandidateVertex> getReachableCandidates() const {
    std::vector<CandidateVertex> result;
    if (L1.isReachable) result.push_back(L1);
    if (L2L.isReachable) result.push_back(L2L);
    if (L3L.isReachable) result.push_back(L3L);
    if (L4L.isReachable) result.push_back(L4L);
    if (L5L.isReachable) result.push_back(L5L);
    if (L5LM.isReachable) result.push_back(L5LM);
    if (L2R.isReachable) result.push_back(L2R);
    if (L3R.isReachable) result.push_back(L3R);
    if (L4R.isReachable) result.push_back(L4R);
    if (L5R.isReachable) result.push_back(L5R);
    if (L5RM.isReachable) result.push_back(L5RM);
    return result;
  }
};

// ============================================================================
// C#: public record VertexPathStep(Hedge.Vertex From, Hedge.Vertex To, bool IsApexJump, Hedge.Face[]? CrossedFaces)
// ============================================================================
struct VertexPathStep {
  // C#: Hedge.Vertex From,
  Vertex from;

  // C#: Hedge.Vertex To,
  Vertex to;

  // C#: bool IsApexJump,
  bool isApexJump = false;

  // C#: Hedge.Face[]? CrossedFaces
  std::vector<Face> crossedFaces;  // Empty = null in C#
};

// ============================================================================
// C#: public readonly struct PathStep(...)
// ============================================================================
struct PathStep {
  // C#: Hedge.Vertex From,
  Vertex from;

  // C#: Hedge.Vertex To,
  Vertex to;

  // C#: bool IsExplorerJump,
  bool isExplorerJump = false;

  // C#: CandidateName CandidateName,
  CandidateName candidateName = CandidateName::None;

  // C#: float Distance,
  double distance = 0.0;

  // C#: Hedge.Face[]? CrossedFaces = null
  std::vector<Face> crossedFaces;  // Empty = null in C#
};

// ============================================================================
// C#: public record PathResult(...)
// ============================================================================
struct PathResult {
  // C#: List<Hedge.Vertex> Path,
  std::vector<Vertex> path;

  // C#: List<PathStep> Steps,
  std::vector<PathStep> steps;

  // C#: bool IsComplete,
  bool isComplete = false;

  // C#: bool IsFallback
  bool isFallback = false;
};

// ============================================================================
// C#: public enum WalkDirection { Clockwise, CounterClockwise }
// ============================================================================
enum class WalkDirection {
  Clockwise,
  CounterClockwise
};

// ============================================================================
// C#: public record WalkResult(...)
// ============================================================================
struct WalkResult {
  // C#: List<Hedge.Face> Faces,
  std::vector<Face> faces;

  // C#: Hedge.Face? FinalFace,
  Face finalFace;  // Invalid face = null

  // C#: bool ReachedTarget,
  bool reachedTarget = false;

  // C#: string? Error
  std::string error = "";

  // Helper to check if finalFace is valid
  bool hasFinalFace() const { return finalFace != Face(); }
};

// A* heuristic weight (1.0 = standard A*)
constexpr double HEURISTIC_WEIGHT = 1.0;

// ============================================================================
// VeryDiscreteGeodesicExplorer functions
// ============================================================================

// Explore reachable vertices from a corner via L5 unfolding
ExplorationResult explore(Corner corner, VertexPositionGeometry& geom);

// ============================================================================
// VeryDiscreteGeodesicExplorerHelper functions
// ============================================================================

// C#: internal static Vector2 FlattenVertex(Hedge.Vertex v, Hedge.Halfedge portal, Vector2 flatA, Vector2 flatB)
Vector2 flattenVertex(Vertex v, Halfedge portal, Vector2 flatA, Vector2 flatB,
                      VertexPositionGeometry& geom);

// C#: internal static Vector2 GetFlatPosition(Hedge.Vertex v, Hedge.Halfedge portal, Vector2 flatP1, Vector2 flatP2, Vector2 flatApex)
Vector2 getFlatPosition(Vertex v, Halfedge portal, Vector2 flatP1, Vector2 flatP2, Vector2 flatApex);

// C#: private static bool SegmentCrossesPortal(Vector2 v0, Vector2 target, Vector2 a, Vector2 b, float eps)
bool segmentCrossesPortal(Vector2 v0, Vector2 target, Vector2 a, Vector2 b, double eps);

// C#: CheckReachability1/2/3/4/5/6/7
CandidateVertex checkReachability1(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, double scaleTolerance);
CandidateVertex checkReachability2(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, Vector2 p2a, Vector2 p2b, double scaleTolerance);
CandidateVertex checkReachability3(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, Vector2 p2a, Vector2 p2b, Vector2 p3a, Vector2 p3b, double scaleTolerance);
CandidateVertex checkReachability4(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, Vector2 p2a, Vector2 p2b, Vector2 p3a, Vector2 p3b,
                                   Vector2 p4a, Vector2 p4b, double scaleTolerance);
CandidateVertex checkReachability5(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, Vector2 p2a, Vector2 p2b, Vector2 p3a, Vector2 p3b,
                                   Vector2 p4a, Vector2 p4b, Vector2 p5a, Vector2 p5b, double scaleTolerance);

// ============================================================================
// VeryDiscreteGeodesicPathfinder functions
// ============================================================================

// C#: public static PathResult FindPath(Hedge.Vertex from, Hedge.Vertex to)
PathResult findPath(Vertex from, Vertex to,
                    ManifoldSurfaceMesh& mesh,
                    VertexPositionGeometry& geom);

// C#: public static (List<Hedge.Vertex> path, List<VertexPathStep> steps) FindGeodesicPath(...)
std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> findGeodesicPath(
    Vertex from, Vertex to,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom);

// C#: public static (List<Hedge.Face> faces, List<Hedge.Vertex> vertexPath) FindFaceStripWithPath(...)
std::pair<std::vector<Face>, std::vector<Vertex>> findFaceStripWithPath(
    Vertex from, Vertex to,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom);

// C#: public static float ComputePathDistance(List<PathStep> steps)
double computePathDistance(const std::vector<PathStep>& steps);

// ============================================================================
// FaceStripWalker functions
// ============================================================================

// C#: public static WalkDirection DetermineWalkDirection(...)
WalkDirection determineWalkDirection(Vertex prev, Vertex current, Vertex next,
                                     VertexPositionGeometry& geom);

// C#: public static WalkResult WalkToOutgoingEdge(...)
WalkResult walkToOutgoingEdge(Face startFace, Vertex vertex, Vertex targetVertex,
                               WalkDirection direction);

// C#: private static WalkResult WalkToFace(...)
WalkResult walkToFace(Face startFace, Face targetFace, Vertex vertex, WalkDirection direction);

// C#: public static Hedge.Face? SelectFirstFace(...)
Face selectFirstFace(Vertex v0, Vertex v1, Vertex v2,
                     VertexPositionGeometry& geom,
                     Face targetFace = Face());

// C#: public static List<Hedge.Face> ConvertPath(...)
std::vector<Face> convertPath(const std::vector<Vertex>& path,
                               const std::vector<VertexPathStep>& steps,
                               Vertex entryVertex, Vertex exitVertex,
                               VertexPositionGeometry& geom);

// ============================================================================
// FaceStripUtils functions
// ============================================================================

// C#: public static Hedge.Face? GetSharedFace(Hedge.Vertex v1, Hedge.Vertex v2)
Face getSharedFace(Vertex v1, Vertex v2);

// C#: public static List<Hedge.Face> GetAllSharedFaces(...)
std::vector<Face> getAllSharedFaces(Vertex v1, Vertex v2);

// C#: public static bool AreFacesAdjacent(Hedge.Face f1, Hedge.Face f2)
bool areFacesAdjacent(Face f1, Face f2);

// C#: public static bool FaceContainsEdge(...)
bool faceContainsEdge(Face face, Vertex v1, Vertex v2);

// C#: public static bool FaceContainsVertex(...)
bool faceContainsVertex(Face face, Vertex vertex);

// C#: public static Hedge.Halfedge? FindHalfedgeToVertex(...)
Halfedge findHalfedgeToVertex(Face face, Vertex vertex);

// C#: public static Hedge.Halfedge? FindHalfedgeFromVertex(...)
Halfedge findHalfedgeFromVertex(Face face, Vertex vertex);

// ============================================================================
// GeodesicPathJoiner functions
// ============================================================================

// C#: public static (List<Hedge.Vertex> path, List<VertexPathStep> steps) JoinPaths(...)
std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> joinPaths(
    const std::vector<std::pair<std::vector<Vertex>, std::vector<VertexPathStep>>>& segments,
    Vertex entryVertex, Vertex exitVertex,
    bool removeLoops = true);

// ============================================================================
// Helper: V2.ComputeTriangleApex from C#
// ============================================================================

// C#: V2.ComputeTriangleApex(Vector2 a, Vector2 b, float distA, float distB, bool pickPositiveY)
Vector2 computeTriangleApex(Vector2 a, Vector2 b, double distA, double distB, bool pickPositiveY);

// ============================================================================
// CachedVeryDiscreteGeodesicPathfinder
//
// Instance-based A* pathfinder with aggressive caching:
// 1. Per-corner exploration results (L5 unfolding, cached per corner)
// 2. Per-vertex boundary edge counts
// 3. Per-vertex adjacency arrays
// 4. Reusable A* containers
// ============================================================================

class CachedVeryDiscreteGeodesicPathfinder {
public:
  CachedVeryDiscreteGeodesicPathfinder(ManifoldSurfaceMesh& mesh,
                                        VertexPositionGeometry& geom);

  // Main API
  std::pair<std::vector<Face>, std::vector<Vertex>> findFaceStripWithPath(
      Vertex from, Vertex to);

  std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> findGeodesicPath(
      Vertex from, Vertex to);

  PathResult findPath(Vertex from, Vertex to);

  // Cache management
  void clearCache();
  void precomputeVertexData();

  // Statistics
  size_t getCacheHits() const { return cacheHits; }
  size_t getCacheMisses() const { return cacheMisses; }
  size_t getExplorationCacheSize() const { return explorationCache.size(); }

private:
  ManifoldSurfaceMesh& mesh;
  VertexPositionGeometry& geom;

  // Per-corner exploration cache (L5 unfolding results)
  std::unordered_map<size_t, ExplorationResult> explorationCache;

  // Per-vertex caches
  VertexData<int> boundaryEdgeCount;
  VertexData<std::vector<Vertex>> adjacentVertices;
  VertexData<std::vector<Corner>> adjacentCorners;
  bool vertexDataComputed = false;

  // Reusable A* containers
  std::unordered_map<size_t, double> gScore;
  std::unordered_map<size_t, std::pair<Vertex, PathStep>> cameFrom;
  std::unordered_set<size_t> closed;
  std::unordered_set<size_t> seenVertices;
  std::vector<std::tuple<Vertex, double, PathStep>> neighbors;

  // Statistics
  size_t cacheHits = 0;
  size_t cacheMisses = 0;

  // Internal helpers
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
