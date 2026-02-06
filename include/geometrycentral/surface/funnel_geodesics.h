#pragma once

// Funnel Geodesics (GFR) - Greedy Funnel Refinement
//
// An alternative to FlipOut for computing geodesic paths on triangle meshes.
// Uses sleeve refinement and quality-gated corner flipping instead of
// intrinsic edge flipping.
//
// Reference: See docs/GFR_PORT_PLAN.md for algorithm details.
//
// Ported from C# implementation in Colonel.Meshing.GreedyFunnelRefinement

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"
// C# counterparts:
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/CachedGreedyFunnelRefinementPathfinder.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/FlipOutComparisonTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/GeodesicComparisonTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/GeodesicTestBase.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/GreedyFunnelRefinementPathfinder.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/IterativeStraightener/CornerAnalyzerRenderer.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/IterativeStraightener/InterativeStraightener.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/IterativeStraightener/InterativeStraightenerTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/IterativeStraightener/PartialFunnelRecomputationTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/IterativeStraightener/WaypointCorner.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/IterativeStraightener/WaypointCornerAnalyzer.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/IterativeStraightener/WaypointCornerAnalyzerTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/IterativeStraightener/WaypointCornerFlipAction.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/RandomPathStressTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/DepthComparisonTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/FlatVertex.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/FlipActionVisualizationTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/FlipRoundTripTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/GeodesicPathJoiner.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/Portal.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/PortalCrossing.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/ShortestEdgePathfinder.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/Sleeve.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/SleeveBadPathsFinderTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/SleeveFlipVisualizationTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/SleevePerformanceTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/SleeveRegressionTests.cs
// - C:/Dev/Colonel/Colonel.Meshing/GreedyFunnelRefinement/SleeveBuilding/SleeveTests.cs


#include <memory>
#include <vector>

namespace geometrycentral {
namespace surface {

// Forward declarations
class FunnelGeodesicPath;

// ============================================================================
// Main API
// ============================================================================

/// Strategy for selecting which corner to flip during straightening.
enum class FlipStrategy {
  /// Always pick the globally most acute corner (quality-gated).
  CornerGreedy,
  /// Flip most acute corner, then flip newly exposed vertices once (ungated).
  WedgeGreedy,
  /// Quality-gated corner greedy with a short coherent follow-up window.
  CoherentMiniWedge
};

/// Options for GFR straightening behavior.
struct FunnelGeodesicOptions {
  FlipStrategy strategy = FlipStrategy::CornerGreedy;
  size_t maxIterations = 20000;
  size_t coherentMiniWedgeDepth = 1;
  size_t coherentMiniWedgeMaxCandidates = 64;
};

/// Compute a geodesic path between two vertices using Funnel Refinement.
/// This is the main entry point - mirrors FlipEdgeNetwork::constructFromDijkstraPath()
std::unique_ptr<FunnelGeodesicPath> computeFunnelGeodesic(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex startVert,
    Vertex endVert);

/// Compute a geodesic path with custom options (strategy, limits).
std::unique_ptr<FunnelGeodesicPath> computeFunnelGeodesic(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex startVert,
    Vertex endVert,
    const FunnelGeodesicOptions& options);

/// Get cache statistics for the VeryDiscreteGeodesic pathfinder (for debugging)
struct CacheStats {
  size_t hits = 0;
  size_t misses = 0;
  size_t cacheSize = 0;
};
CacheStats getCacheStats();

/// Clear the internal cache (for cold-cache benchmarking)
void clearFunnelCache();

/// Get timing breakdown for profiling
struct TimingStats {
  double aStarMs = 0;
  double flattenMs = 0;
  double straightenMs = 0;
  size_t pathCount = 0;
  // HYPOTHESIS 2: Flip discovery stats
  size_t flipAttempts = 0;
  size_t flipFailures = 0;
};
TimingStats getTimingStats();
void resetTimingStats();

// ============================================================================
// Result class
// ============================================================================

class FunnelGeodesicPath {
public:
  FunnelGeodesicPath(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom);

  // Query results
  double length() const;
  const std::vector<SurfacePoint>& getPath() const;

  // Statistics
  size_t iterationCount() const;
  size_t faceCount() const;  // faces in final sleeve
  double initialFunnelLength() const;  // length before any straightening

  // For visualization (optional)
  VertexPositionGeometry* posGeom = nullptr;
  std::vector<Vector3> getPathPolyline3D() const;

  // For debugging
  const std::vector<Face>& getSleeveFaces() const { return sleeveFaces; }

private:
  ManifoldSurfaceMesh& mesh;
  VertexPositionGeometry& geom;

  // Results
  std::vector<SurfacePoint> pathPoints;
  double pathLength = 0.0;
  double initialLength = 0.0;  // funnel distance before any straightening
  size_t nIterations = 0;
  size_t nFaces = 0;

  // Internal: Sleeve data
  std::vector<Face> sleeveFaces;

  friend std::unique_ptr<FunnelGeodesicPath> computeFunnelGeodesic(
      ManifoldSurfaceMesh&, VertexPositionGeometry&, Vertex, Vertex);
  friend std::unique_ptr<FunnelGeodesicPath> computeFunnelGeodesic(
      ManifoldSurfaceMesh&, VertexPositionGeometry&, Vertex, Vertex, const FunnelGeodesicOptions&);
};

// ============================================================================
// Internal structures (exposed for testing)
// ============================================================================

namespace funnel_internal {

/// 2D position of a vertex in the unfolded sleeve
struct FlatVertex {
  Vertex vertex;
  Vector2 pos;
};

/// Portal between two adjacent faces in the sleeve
struct Portal {
  Vector2 left, right;
  Vertex leftVert, rightVert;
};

/// Result of funnel algorithm
struct FunnelResult {
  std::vector<Vector2> waypoints2D;
  std::vector<Vertex> waypointVertices;
  double distance;
};

/// Corner that may want to flip
struct WaypointCorner {
  Vertex vertex;
  size_t faceIndex;
  double angleErrorDeg;
  // C#: public bool WantToFlip() => AngleErrorDeg > 0 && CanStraighten;
  // Note: CanStraighten is checked in computeFlipAction (canFlip)
  bool wantsToFlip() const { return angleErrorDeg > 0; }
};

/// Action to flip a corner (add/remove faces)
struct FlipAction {
  std::vector<Face> removeFaces;
  std::vector<Face> addFaces;
  size_t spliceAfterIndex;  // Insert addFaces after this face index
  size_t spliceBeforeIndex; // And before this face index
  bool canFlip;
};

// Phase 1: Flatten face strip to 2D
VertexData<Vector2> flattenSleeve(
    const std::vector<Face>& faces,
    Vertex entry,
    VertexPositionGeometry& geom);

// Phase 2: Build portals from flattened sleeve
std::vector<Portal> buildPortals(
    const std::vector<Face>& faces,
    const VertexData<Vector2>& flatPos);

// Phase 3: Run funnel algorithm
FunnelResult runFunnel(
    const std::vector<Portal>& portals,
    Vector2 entry,
    Vector2 exit);

// Phase 4: Build face strip from Dijkstra path
std::vector<Face> buildFaceStrip(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex start,
    Vertex end);

// Phase 4b: Build face strip via VeryDiscreteGeodesic (better initial corridors)
std::vector<Face> buildFaceStripVeryDiscrete(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex start,
    Vertex end);

// Phase 5: Analyze corners
std::vector<WaypointCorner> analyzeCorners(
    const std::vector<Face>& faces,
    const FunnelResult& funnel,
    const VertexData<Vector2>& flatPos);

// Phase 6: Compute flip action
FlipAction computeFlipAction(
    const std::vector<Face>& faces,
    const WaypointCorner& corner,
    ManifoldSurfaceMesh& mesh);

// Phase 7: Apply flip (returns new face list)
std::vector<Face> applyFlip(
    const std::vector<Face>& faces,
    const FlipAction& action);

} // namespace funnel_internal

} // namespace surface
} // namespace geometrycentral
