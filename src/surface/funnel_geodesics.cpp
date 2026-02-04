// Funnel Geodesics - GFR Algorithm Implementation
// Computes geodesic paths via sleeve unfolding and iterative corner straightening.

#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/very_discrete_geodesic.h"

#include <chrono>
#include <cmath>
#include <set>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace geometrycentral {
namespace surface {

// Timing accumulators for profiling
namespace {
  double totalAStarTimeMs = 0;
  double totalFlattenTimeMs = 0;
  double totalStraightenTimeMs = 0;
  size_t totalPathsComputed = 0;
  size_t totalFlipAttempts = 0;
  size_t totalFlipFailures = 0;
}

TimingStats getTimingStats() {
  TimingStats stats;
  stats.aStarMs = totalAStarTimeMs;
  stats.flattenMs = totalFlattenTimeMs;
  stats.straightenMs = totalStraightenTimeMs;
  stats.pathCount = totalPathsComputed;
  stats.flipAttempts = totalFlipAttempts;
  stats.flipFailures = totalFlipFailures;
  return stats;
}

void resetTimingStats() {
  totalAStarTimeMs = 0;
  totalFlattenTimeMs = 0;
  totalStraightenTimeMs = 0;
  totalPathsComputed = 0;
  totalFlipAttempts = 0;
  totalFlipFailures = 0;
}

// ============================================================================
// Main API
// ============================================================================

std::unique_ptr<FunnelGeodesicPath> computeFunnelGeodesic(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex startVert,
    Vertex endVert) {

  auto result = std::make_unique<FunnelGeodesicPath>(mesh, geom);
  totalPathsComputed++;

  // Phase 4: Build initial face strip via VeryDiscreteGeodesic (TIMED)
  auto aStarStart = std::chrono::high_resolution_clock::now();
  result->sleeveFaces = funnel_internal::buildFaceStripVeryDiscrete(mesh, geom, startVert, endVert);
  result->nFaces = result->sleeveFaces.size();
  auto aStarEnd = std::chrono::high_resolution_clock::now();
  totalAStarTimeMs += std::chrono::duration<double, std::milli>(aStarEnd - aStarStart).count();

  // Handle degenerate case: empty sleeve
  // This can happen if vertices share an edge (valid), or if pathfinding failed (bug)
  if (result->sleeveFaces.empty()) {
    // Check if vertices share an edge directly
    bool sharesEdge = false;
    double edgeLength = 0;
    for (Halfedge he : startVert.outgoingHalfedges()) {
      if (he.tipVertex() == endVert) {
        sharesEdge = true;
        edgeLength = geom.edgeLengths[he.edge()];
        break;
      }
    }

    result->pathPoints.push_back(SurfacePoint(startVert));
    if (startVert != endVert) {
      result->pathPoints.push_back(SurfacePoint(endVert));
      if (sharesEdge) {
        // Adjacent vertices - return edge length
        result->pathLength = edgeLength;
      } else {
        // Pathfinding failed - this shouldn't happen on connected meshes
        // Fall back to Euclidean as a last resort
        result->pathLength = (geom.vertexPositions[startVert] - geom.vertexPositions[endVert]).norm();
      }
    } else {
      result->pathLength = 0;
    }
    result->nIterations = 0;
    result->nFaces = 0;
    return result;
  }

  // Phase 1-3: Flatten, portals, funnel (TIMED as "flatten")
  auto flattenStart = std::chrono::high_resolution_clock::now();
  auto flatPos = funnel_internal::flattenSleeve(result->sleeveFaces, startVert, geom);
  auto portals = funnel_internal::buildPortals(result->sleeveFaces, flatPos);
  Vector2 entry2D = flatPos[startVert];
  Vector2 exit2D = flatPos[endVert];
  auto funnel = funnel_internal::runFunnel(portals, entry2D, exit2D);
  result->initialLength = funnel.distance;  // Store initial funnel length before straightening
  auto flattenEnd = std::chrono::high_resolution_clock::now();
  totalFlattenTimeMs += std::chrono::duration<double, std::milli>(flattenEnd - flattenStart).count();

  // Straightening start (TIMED)
  auto straightenStart = std::chrono::high_resolution_clock::now();

  // Phase 8: Iterative straightening with phase-based optimization
  // Port of C# InterativeStraightener.StraightenWithPhases()
  //
  // Each phase tries all corners; phases repeat until no improvement.
  // This guarantees local optimality: if any single flip could improve the path,
  // we would have found it.
  //
  // Use phase-stamped array instead of std::set for O(1) lookups
  // Start phase at 1 so initial check (rejectedPhase[idx] != phase) works correctly
  std::vector<uint32_t> rejectedPhase(mesh.nVertices(), 0);
  const size_t maxIters = 20000;

  // Use relative thresholds based on path length to handle paths of all scales
  // Accept any improvement > 0 (floating point epsilon for numerical stability)
  // Mark corner as "done" if relative improvement is negligible (< 0.001% of path length)
  const double relativeNegligible = 1e-5;  // 0.001% relative improvement considered negligible

  size_t iteration = 0;
  uint32_t phase = 1;

  // Initial analysis
  auto corners = funnel_internal::analyzeCorners(result->sleeveFaces, funnel, flatPos);

  // Helper lambda to check if there are flippable corners
  auto hasFlippableCorners = [&]() {
    for (const auto& c : corners) {
      if (c.wantsToFlip() && rejectedPhase[c.vertex.getIndex()] != phase) {
        return true;
      }
    }
    return false;
  };

  // Helper lambda to select best corner
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

  // Check if already optimal (no corners want to flip)
  if (!hasFlippableCorners()) {
    // Already optimal - nothing to do
  } else {
    // Main optimization loop
    while (iteration < maxIters) {
      phase++;  // Incrementing phase implicitly "clears" rejectedPhase for this phase
      double phaseStartDistance = funnel.distance;

      // Run one phase - try all corners until none want to flip (or are all rejected)
      while (iteration < maxIters) {
        iteration++;

        auto* cornerToFlip = selectCornerToFlip();
        if (cornerToFlip == nullptr) {
          break;  // Phase complete - no more corners to try
        }

        auto action = funnel_internal::computeFlipAction(result->sleeveFaces, *cornerToFlip, mesh);
        totalFlipAttempts++;
        if (!action.canFlip) {
          // Corner can't be flipped (boundary constraint) - reject for this phase
          totalFlipFailures++;
          rejectedPhase[cornerToFlip->vertex.getIndex()] = phase;
          continue;
        }

        size_t idx = cornerToFlip->vertex.getIndex();
        double distanceBefore = funnel.distance;

        // Apply flip speculatively
        auto newFaces = funnel_internal::applyFlip(result->sleeveFaces, action);
        auto newFlatPos = funnel_internal::flattenSleeve(newFaces, startVert, geom);
        auto newPortals = funnel_internal::buildPortals(newFaces, newFlatPos);
        auto newFunnel = funnel_internal::runFunnel(newPortals, entry2D, newFlatPos[endVert]);

        // Accept any improvement (flip is speculative, so no harm if we reject)
        double actualImprovement = distanceBefore - newFunnel.distance;

        if (actualImprovement <= 0) {
          // No improvement or worsening - reject for this phase
          rejectedPhase[idx] = phase;
          continue;
        }

        // Good flip - keep it
        result->sleeveFaces = std::move(newFaces);
        flatPos = std::move(newFlatPos);
        portals = std::move(newPortals);
        funnel = std::move(newFunnel);
        exit2D = flatPos[endVert];
        result->nIterations++;

        // Re-analyze corners after flip
        corners = funnel_internal::analyzeCorners(result->sleeveFaces, funnel, flatPos);

        // If relative improvement is negligible, mark corner as rejected to prevent infinite loops
        // This catches cases where floating point drift causes tiny "improvements" that go nowhere
        double relativeImprovement = actualImprovement / distanceBefore;
        if (relativeImprovement < relativeNegligible) {
          rejectedPhase[idx] = phase;
        }
      }

      // Phase complete - check for improvement (relative threshold)
      double phaseEndDistance = funnel.distance;
      double phaseImprovement = (phaseStartDistance - phaseEndDistance) / phaseStartDistance;
      bool improved = phaseImprovement > relativeNegligible;

      // Re-analyze corners to check if any want to flip
      corners = funnel_internal::analyzeCorners(result->sleeveFaces, funnel, flatPos);
      bool hasCorners = false;
      for (const auto& c : corners) {
        if (c.wantsToFlip()) {
          hasCorners = true;
          break;
        }
      }

      if (!improved || !hasCorners) {
        // No improvement this phase, or no more corners - we're locally optimal
        break;
      }

      // Phase improved - continue to next phase
      // (rejectedThisPhase.clear() happens at top of outer loop)
    }
  }

  // End straightening timing
  auto straightenEnd = std::chrono::high_resolution_clock::now();
  totalStraightenTimeMs += std::chrono::duration<double, std::milli>(straightenEnd - straightenStart).count();

  result->pathLength = funnel.distance;
  result->nFaces = result->sleeveFaces.size();

  // Convert waypoints to SurfacePoints
  // Path: startVert -> waypointVertices -> endVert
  result->pathPoints.push_back(SurfacePoint(startVert));
  for (Vertex v : funnel.waypointVertices) {
    result->pathPoints.push_back(SurfacePoint(v));
  }
  result->pathPoints.push_back(SurfacePoint(endVert));

  result->posGeom = &geom;

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

double FunnelGeodesicPath::initialFunnelLength() const { return initialLength; }

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
// Helper: Circle intersection for 2D unfolding
// Given two circles centered at c1 (radius r1) and c2 (radius r2),
// return the intersection point on the opposite side from refPoint.
// ----------------------------------------------------------------------------
static Vector2 circleIntersect(Vector2 c1, double r1, Vector2 c2, double r2, Vector2 refPoint) {
  double dx = c2.x - c1.x;
  double dy = c2.y - c1.y;
  double d = std::sqrt(dx * dx + dy * dy);

  if (d < 1e-14) {
    return Vector2{c1.x + r1, c1.y};
  }

  // Handle non-intersecting circles (degenerate case)
  if (d > r1 + r2 + 1e-10 || d < std::abs(r1 - r2) - 1e-10) {
    return Vector2{c1.x + dx / d * r1, c1.y + dy / d * r1};
  }

  // Standard circle intersection formula
  double a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
  double hSq = r1 * r1 - a * a;
  if (hSq < 0) hSq = 0;
  double h = std::sqrt(hSq);

  double invD = 1.0 / d;
  double dirX = dx * invD;
  double dirY = dy * invD;
  double midX = c1.x + dirX * a;
  double midY = c1.y + dirY * a;

  if (h < 1e-10) {
    return Vector2{midX, midY};
  }

  double perpX = -dirY * h;
  double perpY = dirX * h;

  double p1X = midX + perpX;
  double p1Y = midY + perpY;

  // Pick the side opposite to refPoint
  double refSide = dx * (refPoint.y - c1.y) - dy * (refPoint.x - c1.x);
  double p1Side = dx * (p1Y - c1.y) - dy * (p1X - c1.x);

  if (p1Side * refSide < 0) {
    return Vector2{p1X, p1Y};
  } else {
    return Vector2{midX - perpX, midY - perpY};
  }
}

// Helper: Circle intersection picking positive Y side (for first face)
static Vector2 circleIntersectPositiveY(Vector2 c1, double r1, Vector2 c2, double r2) {
  double dx = c2.x - c1.x;
  double dy = c2.y - c1.y;
  double d = std::sqrt(dx * dx + dy * dy);

  if (d < 1e-14) {
    return Vector2{c1.x + r1, c1.y};
  }

  if (d > r1 + r2 + 1e-10 || d < std::abs(r1 - r2) - 1e-10) {
    return Vector2{c1.x + dx / d * r1, c1.y + dy / d * r1};
  }

  double a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
  double hSq = r1 * r1 - a * a;
  if (hSq < 0) hSq = 0;
  double h = std::sqrt(hSq);

  double invD = 1.0 / d;
  double dirX = dx * invD;
  double dirY = dy * invD;
  double midX = c1.x + dirX * a;
  double midY = c1.y + dirY * a;

  if (h < 1e-10) {
    return Vector2{midX, midY};
  }

  double perpX = -dirY * h;
  double perpY = dirX * h;

  // Pick positive Y
  return Vector2{midX + perpX, midY + perpY};
}

// ----------------------------------------------------------------------------
// Phase 1: Flatten sleeve to 2D
// Port of FunnelAlgorithm.FlattenFaceStripD() from C#
// ----------------------------------------------------------------------------
VertexData<Vector2> flattenSleeve(
    const std::vector<Face>& faces,
    Vertex entry,
    VertexPositionGeometry& geom) {

  VertexData<Vector2> flatPos(geom.mesh, Vector2::zero());
  VertexData<char> known(geom.mesh, 0);  // Track which vertices have been positioned

  if (faces.empty()) {
    return flatPos;
  }

  // Get 3D vertex positions
  const VertexData<Vector3>& positions = geom.vertexPositions;

  // Flatten first face with entry vertex at origin
  Face firstFace = faces[0];
  Halfedge he = firstFace.halfedge();
  Vertex v0, v1, v2;

  // Find the entry vertex and order the other two
  if (he.vertex() == entry) {
    v0 = he.vertex();
    v1 = he.next().vertex();
    v2 = he.next().next().vertex();
  } else if (he.next().vertex() == entry) {
    v0 = he.next().vertex();
    v1 = he.next().next().vertex();
    v2 = he.vertex();
  } else {
    v0 = he.next().next().vertex();
    v1 = he.vertex();
    v2 = he.next().vertex();
  }

  // Place v0 at origin, v1 along +X axis
  Vector2 p0 = Vector2::zero();
  double d01 = (positions[v0] - positions[v1]).norm();
  Vector2 p1 = Vector2{d01, 0};

  // Place v2 using circle intersection (positive Y side)
  double r02 = (positions[v0] - positions[v2]).norm();
  double r12 = (positions[v1] - positions[v2]).norm();
  Vector2 p2 = circleIntersectPositiveY(p0, r02, p1, r12);

  flatPos[v0] = p0;
  flatPos[v1] = p1;
  flatPos[v2] = p2;
  known[v0] = 1;
  known[v1] = 1;
  known[v2] = 1;

  if (faces.size() == 1) {
    return flatPos;
  }

  // Determine initial reference point (the vertex that "falls off" from face 0 to face 1)
  Face face1 = faces[1];
  Halfedge f1he = face1.halfedge();
  Vertex f1v0 = f1he.vertex();
  Vertex f1v1 = f1he.next().vertex();
  Vertex f1v2 = f1he.next().next().vertex();

  bool v0InF1 = (v0 == f1v0 || v0 == f1v1 || v0 == f1v2);
  bool v1InF1 = (v1 == f1v0 || v1 == f1v1 || v1 == f1v2);

  Vector2 refPos;
  if (!v0InF1) refPos = p0;
  else if (!v1InF1) refPos = p1;
  else refPos = p2;

  // Process remaining faces
  for (size_t i = 1; i < faces.size(); i++) {
    Face face = faces[i];
    Halfedge fhe = face.halfedge();
    Vertex fv0 = fhe.vertex();
    Vertex fv1 = fhe.next().vertex();
    Vertex fv2 = fhe.next().next().vertex();

    // Find which vertex is new (not yet positioned) - O(1) lookup via known mask
    bool fv0Known = known[fv0];
    bool fv1Known = known[fv1];
    bool fv2Known = known[fv2];

    Vertex newVert;
    Vertex hinge1, hinge2;

    if (!fv0Known) {
      newVert = fv0;
      hinge1 = fv1;
      hinge2 = fv2;
    } else if (!fv1Known) {
      newVert = fv1;
      hinge1 = fv0;
      hinge2 = fv2;
    } else if (!fv2Known) {
      newVert = fv2;
      hinge1 = fv0;
      hinge2 = fv1;
    } else {
      // All three vertices already known - skip
      continue;
    }

    // Compute new vertex position using circle intersection
    Vector2 h1Pos = flatPos[hinge1];
    Vector2 h2Pos = flatPos[hinge2];
    double r1 = (positions[hinge1] - positions[newVert]).norm();
    double r2 = (positions[hinge2] - positions[newVert]).norm();

    Vector2 newPos = circleIntersect(h1Pos, r1, h2Pos, r2, refPos);
    flatPos[newVert] = newPos;
    known[newVert] = 1;

    // Update refPos for next iteration: find vertex that "falls off"
    if (i + 1 < faces.size()) {
      Face nextFace = faces[i + 1];
      Halfedge nhe = nextFace.halfedge();
      Vertex nv0 = nhe.vertex();
      Vertex nv1 = nhe.next().vertex();
      Vertex nv2 = nhe.next().next().vertex();

      bool fv0InNext = (fv0 == nv0 || fv0 == nv1 || fv0 == nv2);
      bool fv1InNext = (fv1 == nv0 || fv1 == nv1 || fv1 == nv2);
      bool fv2InNext = (fv2 == nv0 || fv2 == nv1 || fv2 == nv2);

      if (!fv0InNext) refPos = flatPos[fv0];
      else if (!fv1InNext) refPos = flatPos[fv1];
      else if (!fv2InNext) refPos = flatPos[fv2];
    }
  }

  return flatPos;
}

// ----------------------------------------------------------------------------
// Phase 2: Build portals
// Port of Sleeve.BuildPortalsFullList() from C#
// ----------------------------------------------------------------------------
std::vector<Portal> buildPortals(
    const std::vector<Face>& faces,
    const VertexData<Vector2>& flatPos) {

  std::vector<Portal> portals;

  if (faces.size() < 2) {
    return portals;
  }

  portals.reserve(faces.size() - 1);

  // Find entry vertex (the one at origin)
  Vertex entryVert;
  for (Vertex v : faces[0].adjacentVertices()) {
    if (flatPos[v].norm() < 1e-10) {
      entryVert = v;
      break;
    }
  }
  Vector2 prevPoint = flatPos[entryVert];

  for (size_t i = 0; i < faces.size() - 1; i++) {
    Face currentFace = faces[i];
    Face nextFace = faces[i + 1];

    // Find shared halfedge
    Halfedge sharedHe;
    bool found = false;
    for (Halfedge he : currentFace.adjacentHalfedges()) {
      if (he.twin().face() == nextFace) {
        sharedHe = he;
        found = true;
        break;
      }
    }

    if (!found) {
      // Faces not adjacent - this shouldn't happen in a valid sleeve
      continue;
    }

    Vertex v1 = sharedHe.vertex();
    Vertex v2 = sharedHe.twin().vertex();

    Vector2 p1 = flatPos[v1];
    Vector2 p2 = flatPos[v2];

    // Determine left/right based on travel direction
    Vector2 mid = (p1 + p2) * 0.5;
    Vector2 dir = mid - prevPoint;
    double dirLen = dir.norm();
    if (dirLen > 1e-10) {
      dir = dir / dirLen;
    }

    // Cross product to determine which side each vertex is on
    double cross1 = dir.x * (p1.y - prevPoint.y) - dir.y * (p1.x - prevPoint.x);
    double cross2 = dir.x * (p2.y - prevPoint.y) - dir.y * (p2.x - prevPoint.x);

    Portal portal;
    if (cross1 > cross2) {
      portal.left = p1;
      portal.right = p2;
      portal.leftVert = v1;
      portal.rightVert = v2;
    } else {
      portal.left = p2;
      portal.right = p1;
      portal.leftVert = v2;
      portal.rightVert = v1;
    }
    portals.push_back(portal);

    prevPoint = mid;
  }

  return portals;
}

// ----------------------------------------------------------------------------
// Helper: Signed triangle area (2x area, used for orientation tests)
// Positive = counter-clockwise, negative = clockwise
// ----------------------------------------------------------------------------
static inline double triArea2D(Vector2 a, Vector2 b, Vector2 c) {
  return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}

// ----------------------------------------------------------------------------
// Phase 3: Funnel algorithm
// Port of FunnelAlgorithm.PopulateEager() from C# (simplified version)
// Lee-Preparata algorithm for shortest path through a simple polygon
// ----------------------------------------------------------------------------
FunnelResult runFunnel(
    const std::vector<Portal>& portals,
    Vector2 entry,
    Vector2 exit) {

  const double CROSS_EPS = 1e-6;

  FunnelResult result;
  result.distance = 0.0;

  // Trivial case: no portals
  if (portals.empty()) {
    result.waypoints2D.push_back(entry);
    result.waypoints2D.push_back(exit);
    result.distance = (exit - entry).norm();
    return result;
  }

  std::vector<Vector2> path2D;
  path2D.push_back(entry);

  Vector2 apex = entry;

  // Initialize funnel from first portal
  Vector2 leftPos = portals[0].left;
  Vector2 rightPos = portals[0].right;
  Vertex leftVertex = portals[0].leftVert;
  Vertex rightVertex = portals[0].rightVert;
  size_t leftIndex = 0;
  size_t rightIndex = 0;

  for (size_t i = 1; i <= portals.size(); i++) {
    Vector2 newLeft, newRight;
    Vertex newLeftVertex, newRightVertex;

    if (i < portals.size()) {
      newLeft = portals[i].left;
      newRight = portals[i].right;
      newLeftVertex = portals[i].leftVert;
      newRightVertex = portals[i].rightVert;
    } else {
      // Final "portal" is the exit point
      newLeft = exit;
      newRight = exit;
      newLeftVertex = Vertex();
      newRightVertex = Vertex();
    }

    // Check if funnel is pinched (left == right)
    bool isPinched = (leftPos.x == rightPos.x && leftPos.y == rightPos.y);

    // Check right side narrowing
    double rightNarrowCheck = triArea2D(apex, rightPos, newRight);
    if (rightNarrowCheck >= -CROSS_EPS) {
      bool rightNotNarrowing = (newRight.x == rightPos.x && newRight.y == rightPos.y);
      double leftCross = triArea2D(apex, leftPos, newRight);

      if (isPinched || rightNotNarrowing ||
          (apex.x == rightPos.x && apex.y == rightPos.y) ||
          leftCross <= CROSS_EPS) {
        // Just update right edge
        rightPos = newRight;
        rightVertex = newRightVertex;
        rightIndex = i;
      } else {
        // Funnel closes on left side - add waypoint at left vertex
        result.waypoints2D.push_back(leftPos);
        result.waypointVertices.push_back(leftVertex);
        path2D.push_back(leftPos);

        apex = leftPos;
        i = leftIndex;

        if (i + 1 < portals.size()) {
          leftPos = portals[i + 1].left;
          rightPos = portals[i + 1].right;
          leftVertex = portals[i + 1].leftVert;
          rightVertex = portals[i + 1].rightVert;
        }
        leftIndex = rightIndex = i + 1;
        continue;
      }
    }

    // Check left side narrowing
    double leftNarrow = triArea2D(apex, leftPos, newLeft);
    if (leftNarrow <= CROSS_EPS) {
      double rightCross = triArea2D(apex, rightPos, newLeft);
      bool leftNotNarrowing = (newLeft.x == leftPos.x && newLeft.y == leftPos.y);

      if (isPinched || leftNotNarrowing ||
          (apex.x == leftPos.x && apex.y == leftPos.y) ||
          rightCross >= -CROSS_EPS) {
        // Just update left edge
        leftPos = newLeft;
        leftVertex = newLeftVertex;
        leftIndex = i;
      } else {
        // Funnel closes on right side - add waypoint at right vertex
        result.waypoints2D.push_back(rightPos);
        result.waypointVertices.push_back(rightVertex);
        path2D.push_back(rightPos);

        apex = rightPos;
        i = rightIndex;

        if (i + 1 < portals.size()) {
          leftPos = portals[i + 1].left;
          rightPos = portals[i + 1].right;
          leftVertex = portals[i + 1].leftVert;
          rightVertex = portals[i + 1].rightVert;
        }
        leftIndex = rightIndex = i + 1;
        continue;
      }
    }
  }

  path2D.push_back(exit);

  // Compute total distance
  double distance = 0.0;
  for (size_t i = 0; i < path2D.size() - 1; i++) {
    distance += (path2D[i + 1] - path2D[i]).norm();
  }

  result.distance = distance;
  return result;
}

// ----------------------------------------------------------------------------
// Phase 4: Build face strip from Dijkstra path (Walk-based approach)
// Port of FaceStripWalker.cs - uses halfedge topology to walk around vertices
// ----------------------------------------------------------------------------

enum class WalkDirection { Clockwise, CounterClockwise };

// Determine walk direction based on turn at a vertex (signed angle test)
static WalkDirection determineWalkDirection(
    Vector3 prevPos, Vector3 currPos, Vector3 nextPos, Vector3 normal) {
  Vector3 incoming = prevPos - currPos;
  Vector3 outgoing = nextPos - currPos;

  // Project onto the tangent plane and compute signed angle
  incoming = incoming - normal * dot(incoming, normal);
  outgoing = outgoing - normal * dot(outgoing, normal);

  // Cross product gives signed area - positive if CCW turn, negative if CW
  Vector3 crossProd = cross(incoming, outgoing);
  double signedAngle = dot(crossProd, normal);

  // If left turn (CCW), walk CCW to find the exit edge
  return (signedAngle > 0) ? WalkDirection::CounterClockwise : WalkDirection::Clockwise;
}

// Use shared utilities from very_discrete_geodesic namespace
using very_discrete_geodesic::findHalfedgeToVertex;
using very_discrete_geodesic::findHalfedgeFromVertex;
using very_discrete_geodesic::faceContainsEdge;
using very_discrete_geodesic::getAllSharedFaces;

// Walk around a vertex from startFace until finding a face containing edge to targetVertex
// Returns the faces traversed (not including startFace, not including final face)
static std::vector<Face> walkToOutgoingEdge(
    Face startFace, Vertex vertex, Vertex targetVertex, WalkDirection direction) {

  std::vector<Face> walked;

  // Check if startFace already contains the target edge
  if (faceContainsEdge(startFace, vertex, targetVertex)) {
    return walked;
  }

  // Find starting halfedge
  Halfedge startHe = findHalfedgeToVertex(startFace, vertex);
  if (startHe == Halfedge()) {
    return walked;
  }

  // Move to next face based on direction
  Halfedge currentHe;
  if (direction == WalkDirection::Clockwise) {
    currentHe = startHe.next().twin();
  } else {
    currentHe = startHe.twin();
  }

  const int maxSteps = 100;
  for (int step = 0; step < maxSteps; step++) {
    if (!currentHe.isInterior()) {
      // Hit boundary
      return walked;
    }

    Face currentFace = currentHe.face();

    if (currentFace == startFace) {
      // Wrapped around without finding target edge
      return walked;
    }

    walked.push_back(currentFace);

    // Check if this face contains the target edge
    if (faceContainsEdge(currentFace, vertex, targetVertex)) {
      return walked;
    }

    // Move to next face
    if (direction == WalkDirection::Clockwise) {
      Halfedge leaving = findHalfedgeFromVertex(currentFace, vertex);
      if (leaving == Halfedge()) return walked;
      currentHe = leaving.twin();
    } else {
      Halfedge entering = findHalfedgeToVertex(currentFace, vertex);
      if (entering == Halfedge()) return walked;
      currentHe = entering.twin();
    }
  }

  return walked;
}

std::vector<Face> buildFaceStrip(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex start,
    Vertex end) {

  std::vector<Face> faces;

  if (start == end) {
    for (Face f : start.adjacentFaces()) {
      faces.push_back(f);
      return faces;
    }
    return faces;
  }

  // Use geometry-central's existing Dijkstra to get edge path
  std::vector<Halfedge> edgePath = shortestEdgePath(geom, start, end);

  if (edgePath.empty()) {
    std::cerr << "DEBUG: shortestEdgePath returned empty for V" << start.getIndex()
              << " -> V" << end.getIndex() << std::endl;
    // Start and end on same face
    for (Face f : start.adjacentFaces()) {
      for (Vertex v : f.adjacentVertices()) {
        if (v == end) {
          std::cerr << "DEBUG: Found shared face F" << f.getIndex() << std::endl;
          faces.push_back(f);
          return faces;
        }
      }
    }
    std::cerr << "DEBUG: No shared face found - returning empty" << std::endl;
    return faces;
  }

  // Convert edge path to vertex path
  std::vector<Vertex> vertexPath;
  vertexPath.push_back(start);
  for (Halfedge he : edgePath) {
    vertexPath.push_back(he.tipVertex());
  }

  const VertexData<Vector3>& positions = geom.vertexPositions;

  // Select first face - prefer the one that leads toward the second edge
  Vertex v0 = vertexPath[0];
  Vertex v1 = vertexPath[1];
  auto firstEdgeFaces = getAllSharedFaces(v0, v1);

  if (firstEdgeFaces.empty()) {
    return faces;
  }

  Face currentFace;
  if (vertexPath.size() >= 3 && firstEdgeFaces.size() > 1) {
    // Pick the face that minimizes walk distance to the next edge
    Vertex v2 = vertexPath[2];
    int bestWalkLen = INT_MAX;

    for (Face candFace : firstEdgeFaces) {
      // Check if this face directly contains the next edge
      if (faceContainsEdge(candFace, v1, v2)) {
        currentFace = candFace;
        bestWalkLen = 0;
        break;
      }

      // Try both walk directions and pick shorter
      auto walkCW = walkToOutgoingEdge(candFace, v1, v2, WalkDirection::Clockwise);
      auto walkCCW = walkToOutgoingEdge(candFace, v1, v2, WalkDirection::CounterClockwise);

      int cwLen = faceContainsEdge(candFace, v1, v2) ? 0 :
                  (walkCW.empty() ? INT_MAX : static_cast<int>(walkCW.size()));
      int ccwLen = faceContainsEdge(candFace, v1, v2) ? 0 :
                   (walkCCW.empty() ? INT_MAX : static_cast<int>(walkCCW.size()));

      // Check if walk actually reaches target
      if (!walkCW.empty() && !faceContainsEdge(walkCW.back(), v1, v2)) cwLen = INT_MAX;
      if (!walkCCW.empty() && !faceContainsEdge(walkCCW.back(), v1, v2)) ccwLen = INT_MAX;

      int minLen = std::min(cwLen, ccwLen);
      if (minLen < bestWalkLen) {
        bestWalkLen = minLen;
        currentFace = candFace;
      }
    }

    if (currentFace == Face()) {
      currentFace = firstEdgeFaces[0];
    }
  } else {
    currentFace = firstEdgeFaces[0];
  }

  faces.push_back(currentFace);

  // Process each vertex transition using walk-based approach
  for (size_t i = 1; i + 1 < vertexPath.size(); i++) {
    Vertex prev = vertexPath[i - 1];
    Vertex curr = vertexPath[i];
    Vertex next = vertexPath[i + 1];

    // Check if current face already contains the next edge
    if (faceContainsEdge(currentFace, curr, next)) {
      continue;
    }

    // Determine walk direction based on turn
    Vector3 normal{0, 0, 1};  // Default, will be overridden
    for (Face f : curr.adjacentFaces()) {
      // Use any face's normal as approximation
      Halfedge he = f.halfedge();
      Vector3 p0 = positions[he.vertex()];
      Vector3 p1 = positions[he.next().vertex()];
      Vector3 p2 = positions[he.next().next().vertex()];
      normal = cross(p1 - p0, p2 - p0);
      double len = norm(normal);
      if (len > 1e-10) {
        normal = normal / len;
        break;
      }
    }

    WalkDirection direction = determineWalkDirection(
        positions[prev], positions[curr], positions[next], normal);

    // Walk from current face to find one containing the next edge
    auto walk = walkToOutgoingEdge(currentFace, curr, next, direction);

    // If primary direction failed, try opposite
    if (walk.empty() || !faceContainsEdge(walk.back(), curr, next)) {
      WalkDirection opposite = (direction == WalkDirection::Clockwise)
          ? WalkDirection::CounterClockwise
          : WalkDirection::Clockwise;
      auto walkOpp = walkToOutgoingEdge(currentFace, curr, next, opposite);

      if (!walkOpp.empty() && faceContainsEdge(walkOpp.back(), curr, next)) {
        walk = walkOpp;
      }
    }

    // Add walked faces to the strip
    for (Face f : walk) {
      // Avoid duplicates
      if (faces.empty() || faces.back() != f) {
        faces.push_back(f);
        currentFace = f;
      }
    }
  }

  // Ensure last face contains end vertex
  bool lastHasEnd = false;
  for (Vertex v : faces.back().adjacentVertices()) {
    if (v == end) {
      lastHasEnd = true;
      break;
    }
  }

  if (!lastHasEnd) {
    // Walk to a face containing the end vertex
    Vertex prev = vertexPath[vertexPath.size() - 2];
    for (Face f : end.adjacentFaces()) {
      // Check if adjacent to current face
      for (Halfedge he : faces.back().adjacentHalfedges()) {
        if (he.twin().isInterior() && he.twin().face() == f) {
          faces.push_back(f);
          goto found_end_face;
        }
      }
    }
    found_end_face:;
  }

  return faces;
}

// ----------------------------------------------------------------------------
// Phase 4b: Build face strip via VeryDiscreteGeodesic
// Uses A* with multi-face jumps for better initial corridors
// Uses cached pathfinder for performance (caches explorer results per corner)
// ----------------------------------------------------------------------------

// Static cached pathfinder - reused across calls on the same mesh
static ManifoldSurfaceMesh* cachedMesh = nullptr;
static VertexPositionGeometry* cachedGeom = nullptr;
static std::unique_ptr<very_discrete_geodesic::CachedVeryDiscreteGeodesicPathfinder> cachedPathfinder;

std::vector<Face> buildFaceStripVeryDiscrete(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex start,
    Vertex end) {

  if (start == end) {
    std::vector<Face> faces;
    for (Face f : start.adjacentFaces()) {
      faces.push_back(f);
      return faces;
    }
    return faces;
  }

  // Create or reuse cached pathfinder
  if (cachedMesh != &mesh || cachedGeom != &geom) {
    cachedPathfinder = std::unique_ptr<very_discrete_geodesic::CachedVeryDiscreteGeodesicPathfinder>(
        new very_discrete_geodesic::CachedVeryDiscreteGeodesicPathfinder(mesh, geom));
    cachedMesh = &mesh;
    cachedGeom = &geom;
  }

  // Use cached VeryDiscreteGeodesic to get face strip and path
  std::pair<std::vector<Face>, std::vector<Vertex>> result =
      cachedPathfinder->findFaceStripWithPath(start, end);

  std::vector<Face>& faces = result.first;

  // If VeryDiscreteGeodesic failed, fall back to Dijkstra
  if (faces.empty()) {
    return buildFaceStrip(mesh, geom, start, end);
  }

  // Validate first face contains start vertex
  bool firstHasStart = false;
  for (Vertex v : faces.front().adjacentVertices()) {
    if (v == start) {
      firstHasStart = true;
      break;
    }
  }

  // Validate last face contains end vertex
  bool lastHasEnd = false;
  for (Vertex v : faces.back().adjacentVertices()) {
    if (v == end) {
      lastHasEnd = true;
      break;
    }
  }

  // Fall back to Dijkstra if face strip is invalid
  if (!firstHasStart || !lastHasEnd) {
    return buildFaceStrip(mesh, geom, start, end);
  }

  // Validate connectivity: each consecutive pair must share an edge
  for (size_t i = 1; i < faces.size(); i++) {
    Face f1 = faces[i - 1];
    Face f2 = faces[i];

    bool sharesEdge = false;
    for (Halfedge he : f1.adjacentHalfedges()) {
      if (he.twin().face() == f2) {
        sharesEdge = true;
        break;
      }
    }

    if (!sharesEdge) {
      // Disconnected face strip - fall back to Dijkstra
      return buildFaceStrip(mesh, geom, start, end);
    }
  }

  return faces;
}

// ----------------------------------------------------------------------------
// Phase 5: Analyze corners
// Port of WaypointCornerAnalyzer.Populate() from C#
// ----------------------------------------------------------------------------
std::vector<WaypointCorner> analyzeCorners(
    const std::vector<Face>& faces,
    const FunnelResult& funnel,
    const VertexData<Vector2>& flatPos) {

  std::vector<WaypointCorner> corners;

  if (funnel.waypointVertices.empty()) {
    return corners;
  }

  const double RAD_TO_DEG = 180.0 / M_PI;

  // Build path2D: [entry, wp0, wp1, ..., wpN, exit]
  // We need this to compute angles
  std::vector<Vector2> path2D;

  // Get entry position (the vertex at origin)
  for (Vertex v : faces[0].adjacentVertices()) {
    if (flatPos[v].norm() < 1e-10) {
      path2D.push_back(flatPos[v]);
      break;
    }
  }

  // Add waypoint positions
  for (size_t i = 0; i < funnel.waypoints2D.size(); i++) {
    path2D.push_back(funnel.waypoints2D[i]);
  }

  // Get exit position (vertex in last face that has largest distance from entry)
  Vector2 entry2D = path2D[0];
  Vertex exitVert;
  double maxDist = -1;
  for (Vertex v : faces.back().adjacentVertices()) {
    double d = (flatPos[v] - entry2D).norm();
    if (d > maxDist) {
      maxDist = d;
      exitVert = v;
    }
  }
  path2D.push_back(flatPos[exitVert]);

  // Analyze each waypoint
  for (size_t i = 0; i < funnel.waypointVertices.size(); i++) {
    Vertex vertex = funnel.waypointVertices[i];

    // Get positions in path: prev -> waypoint -> next
    // Path2D is: [entry, wp0, wp1, ..., wpN, exit]
    Vector2 prev = path2D[i];
    Vector2 curr = funnel.waypoints2D[i];
    Vector2 next = path2D[i + 2];

    // Compute turn angle error
    Vector2 incoming = curr - prev;
    Vector2 outgoing = next - curr;
    double cross = incoming.x * outgoing.y - incoming.y * outgoing.x;
    double dot = incoming.x * outgoing.x + incoming.y * outgoing.y;
    double signedAngleRad = std::atan2(cross, dot);
    double angleErrorDeg = std::abs(signedAngleRad * RAD_TO_DEG);

    WaypointCorner corner;
    corner.vertex = vertex;
    corner.faceIndex = 0;  // TODO: Find actual face index
    corner.angleErrorDeg = angleErrorDeg;

    // Find which face in the strip contains this vertex
    for (size_t j = 0; j < faces.size(); j++) {
      for (Vertex v : faces[j].adjacentVertices()) {
        if (v == vertex) {
          corner.faceIndex = j;
          break;
        }
      }
    }

    corners.push_back(corner);
  }

  return corners;
}

// ----------------------------------------------------------------------------
// Helper: Walk halfedge fan around a vertex until reaching exitFace
// Walks in one direction around the vertex, collecting faces until hitting exit or boundary
// ----------------------------------------------------------------------------
static std::vector<Face> walkFanToFace(Halfedge startHe, Face exitFace, bool walkClockwise) {
  std::vector<Face> result;
  Halfedge current = startHe;

  for (int step = 0; step < 50; step++) {
    if (!current.isInterior()) break;

    Face face = current.face();
    if (face == exitFace) {
      // Reached the exit face - success
      return result;
    }

    result.push_back(face);

    // Move to next face around the vertex
    // For CW: go to prev halfedge's twin (prev = next.next for triangles)
    // For CCW: go to next halfedge's twin
    Halfedge neighbor = walkClockwise ? current.next().next() : current.next();
    Halfedge twin = neighbor.twin();
    if (!twin.isInterior()) break;
    current = twin;
  }

  // C# returns partial result (what was collected) even if didn't reach exitFace
  return result;
}

// ----------------------------------------------------------------------------
// Phase 6: Compute flip action
// Port of WaypointCornerFlipAction.Compute() from C#
//
// A corner flip replaces faces around a waypoint vertex with alternate faces
// on the "other side" of the vertex. This allows the path to take a different
// route around the corner.
// ----------------------------------------------------------------------------
FlipAction computeFlipAction(
    const std::vector<Face>& faces,
    const WaypointCorner& corner,
    ManifoldSurfaceMesh& mesh) {

  FlipAction action;
  action.canFlip = false;
  action.spliceAfterIndex = 0;
  action.spliceBeforeIndex = 0;

  Vertex vertex = corner.vertex;

  // Boundary vertices cannot be flipped
  if (vertex.isBoundary()) {
    return action;
  }

  // Find the first and last face in the sleeve that contains this vertex
  size_t firstFaceWithVertex = faces.size();
  size_t lastFaceWithVertex = 0;
  for (size_t i = 0; i < faces.size(); i++) {
    for (Vertex v : faces[i].adjacentVertices()) {
      if (v == vertex) {
        if (i < firstFaceWithVertex) firstFaceWithVertex = i;
        if (i > lastFaceWithVertex) lastFaceWithVertex = i;
        break;
      }
    }
  }

  if (firstFaceWithVertex >= faces.size()) {
    return action;
  }

  // Entry and exit are the first and last faces containing the vertex
  Face entryFace = faces[firstFaceWithVertex];
  Face exitFace = faces[lastFaceWithVertex];

  // Store splice indices for applyFlip
  action.spliceAfterIndex = firstFaceWithVertex;
  action.spliceBeforeIndex = lastFaceWithVertex;

  // If vertex only appears in one face, no flip possible
  if (lastFaceWithVertex <= firstFaceWithVertex) {
    return action;
  }

  // RemoveFaces = faces STRICTLY BETWEEN entry and exit (not including entry or exit)
  for (size_t i = firstFaceWithVertex + 1; i < lastFaceWithVertex; i++) {
    action.removeFaces.push_back(faces[i]);
  }

  // Find the halfedge in entry face that points TO the vertex
  Halfedge heToVertex;
  bool foundHe = false;
  for (Halfedge he : entryFace.adjacentHalfedges()) {
    if (he.tipVertex() == vertex) {
      heToVertex = he;
      foundHe = true;
      break;
    }
  }

  if (!foundHe) {
    return action;
  }

  // Walk both directions of the fan from entry toward exit
  // walkA: start from heToVertex.twin(), walk CW
  // walkB: start from the other halfedge's twin, walk CCW
  Halfedge startA = heToVertex.twin();

  // Find the halfedge from vertex in entryFace (the outgoing one)
  Halfedge heFromVertex;
  for (Halfedge he : entryFace.adjacentHalfedges()) {
    if (he.tailVertex() == vertex) {
      heFromVertex = he;
      break;
    }
  }
  Halfedge startB = heFromVertex.twin();

  std::vector<Face> walkA = walkFanToFace(startA, exitFace, true);
  std::vector<Face> walkB = walkFanToFace(startB, exitFace, false);

  // C# checks overlap with removeFaces only, NOT with entry/exit faces
  // Build removeFaces set (faces strictly between entry and exit)
  std::set<Face> removeFaceSet;
  for (const Face& f : action.removeFaces) {
    removeFaceSet.insert(f);
  }

  bool aValid = !walkA.empty() || (startA.isInterior() && startA.face() == exitFace);
  bool bValid = !walkB.empty() || (startB.isInterior() && startB.face() == exitFace);

  // walkA/walkB are empty if they reached exitFace immediately
  if (startA.isInterior() && startA.face() == exitFace) {
    aValid = true;
    walkA.clear();
  }
  if (startB.isInterior() && startB.face() == exitFace) {
    bValid = true;
    walkB.clear();
  }

  // C# checks: aOverlapsRemove = walkA.Any(f => removeFaceSet.Contains(f))
  // Check for overlap with removeFaces (not the full sleeve range)
  if (aValid) {
    for (Face f : walkA) {
      if (removeFaceSet.count(f)) { aValid = false; break; }
    }
  }
  if (bValid) {
    for (Face f : walkB) {
      if (removeFaceSet.count(f)) { bValid = false; break; }
    }
  }

  // Pick the shorter valid walk
  if (aValid && bValid) {
    action.addFaces = (walkA.size() <= walkB.size()) ? walkA : walkB;
    action.canFlip = true;
  } else if (aValid) {
    action.addFaces = walkA;
    action.canFlip = true;
  } else if (bValid) {
    action.addFaces = walkB;
    action.canFlip = true;
  }

  return action;
}

// ----------------------------------------------------------------------------
// Phase 7: Apply flip
// Port of Sleeve.ComputeNewFaces() from C#
//
// Splices addFaces into the sleeve between spliceAfterIndex and spliceBeforeIndex,
// removing any intermediate faces.
// ----------------------------------------------------------------------------
std::vector<Face> applyFlip(
    const std::vector<Face>& faces,
    const FlipAction& action) {

  if (!action.canFlip) {
    return faces;
  }

  // Use the splice indices from computeFlipAction
  size_t spliceAfter = action.spliceAfterIndex;
  size_t spliceBefore = action.spliceBeforeIndex;

  if (spliceAfter >= faces.size() || spliceBefore >= faces.size() || spliceBefore <= spliceAfter) {
    return faces;
  }

  std::vector<Face> newFaces;
  newFaces.reserve(faces.size() - (spliceBefore - spliceAfter - 1) + action.addFaces.size());

  // Prefix: faces up to and including spliceAfter (entry face)
  for (size_t i = 0; i <= spliceAfter; i++) {
    newFaces.push_back(faces[i]);
  }

  // Add replacement faces
  for (Face f : action.addFaces) {
    newFaces.push_back(f);
  }

  // Suffix: faces from spliceBefore onwards (exit face and beyond)
  for (size_t i = spliceBefore; i < faces.size(); i++) {
    newFaces.push_back(faces[i]);
  }

  return newFaces;
}

} // namespace funnel_internal

CacheStats getCacheStats() {
  CacheStats stats;
  if (funnel_internal::cachedPathfinder) {
    stats.hits = funnel_internal::cachedPathfinder->getCacheHits();
    stats.misses = funnel_internal::cachedPathfinder->getCacheMisses();
    stats.cacheSize = funnel_internal::cachedPathfinder->getExplorationCacheSize();
  }
  return stats;
}

} // namespace surface
} // namespace geometrycentral
