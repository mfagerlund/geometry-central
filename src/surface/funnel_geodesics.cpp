#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"

#include <cmath>
#include <set>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

    // Find which vertex is new (not yet positioned)
    bool fv0Known = (flatPos[fv0] != Vector2::zero() || fv0 == entry);
    bool fv1Known = (flatPos[fv1] != Vector2::zero() || fv1 == entry);
    bool fv2Known = (flatPos[fv2] != Vector2::zero() || fv2 == entry);

    // More robust check: see if vertex appears in any previous face
    for (size_t j = 0; j < i && (!fv0Known || !fv1Known || !fv2Known); j++) {
      Face prevFace = faces[j];
      for (Vertex v : prevFace.adjacentVertices()) {
        if (v == fv0) fv0Known = true;
        if (v == fv1) fv1Known = true;
        if (v == fv2) fv2Known = true;
      }
    }

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
// Phase 4: Build face strip from Dijkstra path
// Converts a halfedge path into an ordered face strip suitable for funnel algorithm
//
// The face strip connects start to end vertex. For each edge in the Dijkstra path,
// we add one of its two adjacent faces. We ensure consecutive faces share an edge.
// ----------------------------------------------------------------------------
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

  // Use geometry-central's existing Dijkstra
  std::vector<Halfedge> edgePath = shortestEdgePath(geom, start, end);

  if (edgePath.empty()) {
    // Start and end on same face
    for (Face f : start.adjacentFaces()) {
      for (Vertex v : f.adjacentVertices()) {
        if (v == end) {
          faces.push_back(f);
          return faces;
        }
      }
    }
    return faces;
  }

  // For each edge in the path, we have two faces. We need to pick a sequence
  // of faces where consecutive faces share an edge.
  //
  // Strategy: For each edge, if the current face touches this edge, stay.
  // Otherwise, cross to the other face.

  // Initialize with a face from the first edge
  Halfedge firstHe = edgePath[0];
  Face currentFace;
  if (firstHe.isInterior()) {
    currentFace = firstHe.face();
  } else {
    currentFace = firstHe.twin().face();
  }
  faces.push_back(currentFace);

  // Process each edge - decide whether to cross to the other face
  for (size_t i = 0; i < edgePath.size(); i++) {
    Halfedge he = edgePath[i];
    Edge e = he.edge();

    // Check if current face touches this edge
    bool currentTouchesEdge = false;
    for (Halfedge fhe : currentFace.adjacentHalfedges()) {
      if (fhe.edge() == e) {
        currentTouchesEdge = true;
        break;
      }
    }

    if (!currentTouchesEdge) {
      // We need to cross. Find a face adjacent to currentFace that touches this edge.
      for (Face adjFace : currentFace.adjacentFaces()) {
        for (Halfedge afhe : adjFace.adjacentHalfedges()) {
          if (afhe.edge() == e) {
            faces.push_back(adjFace);
            currentFace = adjFace;
            currentTouchesEdge = true;
            break;
          }
        }
        if (currentTouchesEdge) break;
      }
    }

    // After processing edge i, we should be on a face touching edge i.
    // For the next iteration, we may need to cross to touch edge i+1.
    // The crossing happens at the shared vertex between edges i and i+1.

    if (i + 1 < edgePath.size()) {
      Halfedge nextHe = edgePath[i + 1];
      Edge nextE = nextHe.edge();

      // Check if current face touches the next edge
      bool touchesNext = false;
      for (Halfedge fhe : currentFace.adjacentHalfedges()) {
        if (fhe.edge() == nextE) {
          touchesNext = true;
          break;
        }
      }

      if (!touchesNext) {
        // Cross to a face that touches both current edge and next edge,
        // or just next edge if needed
        Vertex sharedVertex = he.twin().vertex();  // The vertex connecting edges i and i+1

        // Walk around the shared vertex to find a face touching the next edge
        for (Halfedge vhe : sharedVertex.outgoingHalfedges()) {
          if (!vhe.isInterior()) continue;
          Face candFace = vhe.face();

          for (Halfedge cfhe : candFace.adjacentHalfedges()) {
            if (cfhe.edge() == nextE) {
              // This face touches the next edge
              // Check if it's adjacent to current face
              bool adjToCurrent = false;
              for (Halfedge che : currentFace.adjacentHalfedges()) {
                if (che.twin().isInterior() && che.twin().face() == candFace) {
                  adjToCurrent = true;
                  break;
                }
              }

              if (adjToCurrent && candFace != currentFace) {
                faces.push_back(candFace);
                currentFace = candFace;
              } else if (!adjToCurrent && candFace != currentFace) {
                // Need an intermediate face
                for (Face midFace : currentFace.adjacentFaces()) {
                  // Check if midFace is adjacent to candFace
                  for (Halfedge mhe : midFace.adjacentHalfedges()) {
                    if (mhe.twin().isInterior() && mhe.twin().face() == candFace) {
                      faces.push_back(midFace);
                      faces.push_back(candFace);
                      currentFace = candFace;
                      goto done_crossing;
                    }
                  }
                }
                // Couldn't find intermediate - just add candFace
                faces.push_back(candFace);
                currentFace = candFace;
              }
              goto done_crossing;
            }
          }
        }
        done_crossing:;
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
    for (Face f : end.adjacentFaces()) {
      bool adjacent = false;
      for (Halfedge he : faces.back().adjacentHalfedges()) {
        if (he.twin().isInterior() && he.twin().face() == f) {
          adjacent = true;
          break;
        }
      }
      if (adjacent) {
        faces.push_back(f);
        break;
      }
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
// ----------------------------------------------------------------------------
static std::vector<Face> walkFan(Halfedge startHe, Face exitFace, bool usePrev = true) {
  std::vector<Face> result;
  Halfedge current = startHe;

  for (int step = 0; step < 50; step++) {
    if (!current.isInterior()) break;

    Face face = current.face();
    if (face == exitFace) {
      return result;
    }

    result.push_back(face);

    Halfedge neighbor = usePrev ? current.next().next() : current.next();
    Halfedge twin = neighbor.twin();
    if (!twin.isInterior()) break;
    current = twin;
  }

  return result;
}

// ----------------------------------------------------------------------------
// Phase 6: Compute flip action
// Port of WaypointCornerFlipAction.Compute() from C#
// ----------------------------------------------------------------------------
FlipAction computeFlipAction(
    const std::vector<Face>& faces,
    const WaypointCorner& corner,
    ManifoldSurfaceMesh& mesh) {

  FlipAction action;
  action.canFlip = false;

  Vertex vertex = corner.vertex;

  // Boundary vertices cannot be flipped
  if (vertex.isBoundary()) {
    return action;
  }

  // Find entry/exit faces - the faces before and after the corner vertex in the sleeve
  size_t faceIndex = corner.faceIndex;

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

  if (firstFaceWithVertex >= faces.size() || lastFaceWithVertex < firstFaceWithVertex) {
    return action;
  }

  Face entryFace = faces[firstFaceWithVertex];
  Face exitFace = faces[std::min(lastFaceWithVertex + 1, faces.size() - 1)];

  // If there are no faces strictly between entry and exit, no flip possible
  if (lastFaceWithVertex <= firstFaceWithVertex) {
    return action;
  }

  // RemoveFaces = faces strictly between entry and exit that contain this vertex
  for (size_t i = firstFaceWithVertex + 1; i <= lastFaceWithVertex; i++) {
    action.removeFaces.push_back(faces[i]);
  }

  if (action.removeFaces.empty()) {
    return action;
  }

  // Find the halfedge in entry face that points TO the vertex
  Halfedge heToVertex;
  bool foundHe = false;
  for (Halfedge he : entryFace.adjacentHalfedges()) {
    if (he.next().vertex() == vertex) {
      heToVertex = he;
      foundHe = true;
      break;
    }
  }

  if (!foundHe) {
    return action;
  }

  // Walk both directions of the fan from entry toward exit
  std::vector<Face> walkA = walkFan(heToVertex.twin(), exitFace, true);
  std::vector<Face> walkB = walkFan(heToVertex.next().twin(), exitFace, false);

  // AddFaces = the fan walk that does NOT overlap with removeFaces
  std::set<Face> removeFaceSet(action.removeFaces.begin(), action.removeFaces.end());

  bool aOverlapsRemove = false;
  for (Face f : walkA) {
    if (removeFaceSet.count(f)) {
      aOverlapsRemove = true;
      break;
    }
  }

  bool bOverlapsRemove = false;
  for (Face f : walkB) {
    if (removeFaceSet.count(f)) {
      bOverlapsRemove = true;
      break;
    }
  }

  if (!aOverlapsRemove && !walkA.empty()) {
    action.addFaces = walkA;
    action.canFlip = true;
  } else if (!bOverlapsRemove && !walkB.empty()) {
    action.addFaces = walkB;
    action.canFlip = true;
  }

  return action;
}

// ----------------------------------------------------------------------------
// Phase 7: Apply flip
// Port of Sleeve.ComputeNewFaces() from C#
// ----------------------------------------------------------------------------
std::vector<Face> applyFlip(
    const std::vector<Face>& faces,
    const FlipAction& action) {

  if (!action.canFlip || action.removeFaces.empty()) {
    return faces;
  }

  // Build set of faces to remove
  std::set<Face> removeFaceSet(action.removeFaces.begin(), action.removeFaces.end());

  // Find splice boundaries
  int firstRemoveIdx = -1;
  int lastRemoveIdx = -1;
  for (size_t i = 0; i < faces.size(); i++) {
    if (removeFaceSet.count(faces[i])) {
      if (firstRemoveIdx < 0) firstRemoveIdx = static_cast<int>(i);
      lastRemoveIdx = static_cast<int>(i);
    }
  }

  if (firstRemoveIdx < 0) {
    return faces;  // Nothing to remove
  }

  std::vector<Face> newFaces;
  newFaces.reserve(faces.size() - action.removeFaces.size() + action.addFaces.size());

  // Prefix: faces before first remove
  for (int i = 0; i < firstRemoveIdx; i++) {
    newFaces.push_back(faces[i]);
  }

  // Add replacement faces
  for (Face f : action.addFaces) {
    newFaces.push_back(f);
  }

  // Suffix: faces after last remove
  for (size_t i = lastRemoveIdx + 1; i < faces.size(); i++) {
    newFaces.push_back(faces[i]);
  }

  return newFaces;
}

} // namespace funnel_internal

} // namespace surface
} // namespace geometrycentral
