#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"

namespace geometrycentral {
namespace surface {

// ============================================================================
// Main API
// ============================================================================

std::unique_ptr<FunnelGeodesicPath> computeFunnelGeodesic(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex startVert,
    Vertex endVert) {

  auto result = std::make_unique<FunnelGeodesicPath>(mesh, geom);

  // Phase 4: Build initial face strip via Dijkstra
  result->sleeveFaces = funnel_internal::buildFaceStrip(mesh, geom, startVert, endVert);
  result->nFaces = result->sleeveFaces.size();

  // Phase 1: Flatten to 2D
  auto flatPos = funnel_internal::flattenSleeve(result->sleeveFaces, startVert, geom);

  // Phase 2: Build portals
  auto portals = funnel_internal::buildPortals(result->sleeveFaces, flatPos);

  // Phase 3: Initial funnel
  Vector2 entry2D = flatPos[startVert];
  Vector2 exit2D = flatPos[endVert];
  auto funnel = funnel_internal::runFunnel(portals, entry2D, exit2D);

  // Phase 8: Iterative straightening
  std::set<size_t> rejected;
  const size_t maxIters = 20000;
  const double minImprovement = 0.00001;

  for (size_t iter = 0; iter < maxIters; iter++) {
    // Phase 5: Analyze corners
    auto corners = funnel_internal::analyzeCorners(result->sleeveFaces, funnel, flatPos);

    // Select best corner to flip (most acute, not rejected)
    funnel_internal::WaypointCorner* best = nullptr;
    double bestError = 0;
    for (auto& c : corners) {
      if (!c.wantsToFlip()) continue;
      if (rejected.count(c.vertex.getIndex())) continue;
      if (c.angleErrorDeg > bestError) {
        best = &c;
        bestError = c.angleErrorDeg;
      }
    }

    if (best == nullptr) break;  // Converged

    // Phase 6: Compute flip action
    auto action = funnel_internal::computeFlipAction(result->sleeveFaces, *best, mesh);
    if (!action.canFlip) {
      rejected.insert(best->vertex.getIndex());
      continue;
    }

    double oldDist = funnel.distance;

    // Phase 7: Apply flip speculatively
    auto newFaces = funnel_internal::applyFlip(result->sleeveFaces, action);
    auto newFlatPos = funnel_internal::flattenSleeve(newFaces, startVert, geom);
    auto newPortals = funnel_internal::buildPortals(newFaces, newFlatPos);
    auto newFunnel = funnel_internal::runFunnel(newPortals, entry2D, newFlatPos[endVert]);

    if (newFunnel.distance < oldDist - minImprovement) {
      // Good flip - keep it
      result->sleeveFaces = std::move(newFaces);
      flatPos = std::move(newFlatPos);
      portals = std::move(newPortals);
      funnel = std::move(newFunnel);
      exit2D = flatPos[endVert];
      result->nIterations++;
    } else {
      // Bad flip - reject
      rejected.insert(best->vertex.getIndex());
    }
  }

  result->pathLength = funnel.distance;
  result->nFaces = result->sleeveFaces.size();

  // Convert waypoints to SurfacePoints
  // TODO: Implement proper SurfacePoint conversion

  return result;
}

// ============================================================================
// FunnelGeodesicPath implementation
// ============================================================================

FunnelGeodesicPath::FunnelGeodesicPath(ManifoldSurfaceMesh& mesh_, VertexPositionGeometry& geom_)
    : mesh(mesh_), geom(geom_) {}

double FunnelGeodesicPath::length() const { return pathLength; }

const std::vector<SurfacePoint>& FunnelGeodesicPath::getPath() const { return pathPoints; }

size_t FunnelGeodesicPath::iterationCount() const { return nIterations; }

size_t FunnelGeodesicPath::faceCount() const { return nFaces; }

std::vector<Vector3> FunnelGeodesicPath::getPathPolyline3D() const {
  std::vector<Vector3> result;
  for (const auto& sp : pathPoints) {
    result.push_back(sp.interpolate(geom.vertexPositions));
  }
  return result;
}

// ============================================================================
// Internal implementations
// ============================================================================

namespace funnel_internal {

// ----------------------------------------------------------------------------
// Phase 1: Flatten sleeve to 2D
// ----------------------------------------------------------------------------
VertexData<Vector2> flattenSleeve(
    const std::vector<Face>& faces,
    Vertex entry,
    VertexPositionGeometry& geom) {

  // TODO: Implement - Port from FunnelAlgorithm.FlattenFaceStripD()
  //
  // Algorithm:
  // 1. Place entry vertex at origin (0,0)
  // 2. Place first edge along +X axis
  // 3. For each subsequent face, unfold using edge lengths and angles
  // 4. Use double precision throughout

  VertexData<Vector2> flatPos(geom.mesh, Vector2::zero());

  // STUB: Return empty for now
  return flatPos;
}

// ----------------------------------------------------------------------------
// Phase 2: Build portals
// ----------------------------------------------------------------------------
std::vector<Portal> buildPortals(
    const std::vector<Face>& faces,
    const VertexData<Vector2>& flatPos) {

  // TODO: Implement - Port from Sleeve.BuildPortalsList()
  //
  // Algorithm:
  // 1. For each pair of adjacent faces, find shared edge
  // 2. Determine left/right based on sleeve traversal direction
  // 3. Store 2D positions from flatPos

  std::vector<Portal> portals;

  // STUB: Return empty for now
  return portals;
}

// ----------------------------------------------------------------------------
// Phase 3: Funnel algorithm
// ----------------------------------------------------------------------------
FunnelResult runFunnel(
    const std::vector<Portal>& portals,
    Vector2 entry,
    Vector2 exit) {

  // TODO: Implement - Port from FunnelAlgorithm.Populate()
  //
  // Algorithm (Lee-Preparata):
  // 1. Initialize funnel apex at entry
  // 2. For each portal:
  //    - Try to narrow left edge
  //    - Try to narrow right edge
  //    - If funnel inverts, emit waypoint and restart
  // 3. Connect to exit

  FunnelResult result;
  result.distance = 0.0;

  // STUB: Direct distance for now
  result.distance = (exit - entry).norm();

  return result;
}

// ----------------------------------------------------------------------------
// Phase 4: Build face strip from Dijkstra path
// ----------------------------------------------------------------------------
std::vector<Face> buildFaceStrip(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex start,
    Vertex end) {

  // TODO: Implement - Port from FaceStripBuilder.FindSleeve()
  //
  // Algorithm:
  // 1. Run Dijkstra to get vertex path (use shortestEdgePath from geometry-central)
  // 2. Collect faces adjacent to path edges
  // 3. Order and deduplicate

  std::vector<Face> faces;

  // Use geometry-central's existing Dijkstra
  std::vector<Halfedge> edgePath = shortestEdgePath(geom, start, end);

  if (edgePath.empty()) {
    return faces;  // No path found
  }

  // STUB: Collect faces touching the path
  std::set<Face> faceSet;
  for (Halfedge he : edgePath) {
    faceSet.insert(he.face());
    if (!he.twin().isInterior()) continue;
    faceSet.insert(he.twin().face());
  }

  faces.assign(faceSet.begin(), faceSet.end());
  // TODO: Order faces properly along the path

  return faces;
}

// ----------------------------------------------------------------------------
// Phase 5: Analyze corners
// ----------------------------------------------------------------------------
std::vector<WaypointCorner> analyzeCorners(
    const std::vector<Face>& faces,
    const FunnelResult& funnel,
    const VertexData<Vector2>& flatPos) {

  // TODO: Implement - Port from WaypointCornerAnalyzer.Populate()
  //
  // Algorithm:
  // 1. For each interior waypoint vertex
  // 2. Compute angle in 2D (from prev waypoint, through this, to next)
  // 3. angleError = 180 - angle (deviation from straight)
  // 4. wantsToFlip if angleError > threshold

  std::vector<WaypointCorner> corners;

  // STUB: Return empty for now
  return corners;
}

// ----------------------------------------------------------------------------
// Phase 6: Compute flip action
// ----------------------------------------------------------------------------
FlipAction computeFlipAction(
    const std::vector<Face>& faces,
    const WaypointCorner& corner,
    ManifoldSurfaceMesh& mesh) {

  // TODO: Implement - Port from WaypointCornerFlipAction.Compute()
  //
  // Algorithm:
  // 1. At corner vertex, identify the "wedge" on the short side of the path
  // 2. removeFaces = faces in sleeve that are in the wedge
  // 3. addFaces = faces on opposite side of wedge
  // 4. canFlip = false if we hit mesh boundary

  FlipAction action;
  action.canFlip = false;

  // STUB: Return non-flippable for now
  return action;
}

// ----------------------------------------------------------------------------
// Phase 7: Apply flip
// ----------------------------------------------------------------------------
std::vector<Face> applyFlip(
    const std::vector<Face>& faces,
    const FlipAction& action) {

  // TODO: Implement - Port from Sleeve.ApplyFlip() (simple rebuild path)
  //
  // Algorithm:
  // 1. Find position of removeFaces in the face list
  // 2. Remove them
  // 3. Insert addFaces at same position
  // 4. Return new list

  std::vector<Face> newFaces = faces;

  // STUB: Return unchanged for now
  return newFaces;
}

} // namespace funnel_internal

} // namespace surface
} // namespace geometrycentral
