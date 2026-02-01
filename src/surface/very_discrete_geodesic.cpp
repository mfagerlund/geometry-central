// VeryDiscreteGeodesicPathfinder - Line-by-line port from C#
//
// METHODOLOGY: Each C# line is commented, followed by its C++ equivalent.
// This ensures a 1:1 correspondence with the original C# implementation.
//
// Original C# files from Colonel.Meshing.GreedyFunnelRefinement:
//   - VeryDiscreteGeodesicPathfinder.cs
//   - VeryDiscreteGeodesicExplorer.cs
//   - VeryDiscreteGeodesicExplorerHelper.cs
//   - FaceStripWalker.cs
//   - VertexPathToFaceStripConverter.cs
//   - GeodesicPathJoiner.cs
//   - FaceStripUtils.cs

#include "geometrycentral/surface/very_discrete_geodesic.h"

#include <cmath>
#include <algorithm>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace geometrycentral {
namespace surface {
namespace very_discrete_geodesic {

// ============================================================================
// Helper: V2.ComputeTriangleApex
// C# from Colonel.Math.V2:
// public static Vector2 ComputeTriangleApex(Vector2 a, Vector2 b, float distA, float distB, bool pickPositiveY)
// ============================================================================
Vector2 computeTriangleApex(Vector2 a, Vector2 b, double distA, double distB, bool pickPositiveY) {
  // C#: var ab = b - a;
  Vector2 ab = b - a;

  // C#: var abLen = ab.Length();
  double abLen = norm(ab);

  // C#: if (abLen < 1e-10f) return new Vector2(float.NaN, float.NaN);
  if (abLen < 1e-10) return Vector2{std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN()};

  // C#: var cosA = (distA * distA + abLen * abLen - distB * distB) / (2 * distA * abLen);
  double cosA = (distA * distA + abLen * abLen - distB * distB) / (2.0 * distA * abLen);

  // C#: if (cosA > 1) cosA = 1;
  if (cosA > 1.0) cosA = 1.0;

  // C#: if (cosA < -1) cosA = -1;
  if (cosA < -1.0) cosA = -1.0;

  // C#: var sinA = MathF.Sqrt(1 - cosA * cosA);
  double sinA = std::sqrt(1.0 - cosA * cosA);

  // C#: var abNorm = ab / abLen;
  Vector2 abNorm = ab / abLen;

  // C#: var perpX = -abNorm.Y;
  // C#: var perpY = abNorm.X;
  double perpX = -abNorm.y;
  double perpY = abNorm.x;

  // C#: var localX = distA * cosA;
  double localX = distA * cosA;

  // C#: var localY = distA * sinA;
  double localY = distA * sinA;

  // C#: if (!pickPositiveY) localY = -localY;
  if (!pickPositiveY) localY = -localY;

  // C#: var worldX = a.X + localX * abNorm.X + localY * perpX;
  double worldX = a.x + localX * abNorm.x + localY * perpX;

  // C#: var worldY = a.Y + localX * abNorm.Y + localY * perpY;
  double worldY = a.y + localX * abNorm.y + localY * perpY;

  // C#: return new Vector2(worldX, worldY);
  return Vector2{worldX, worldY};
}

// ============================================================================
// FaceStripUtils - Line-by-line port
// C# from FaceStripUtils.cs
// ============================================================================

// C#: public static Hedge.Face? GetSharedFace(Hedge.Vertex v1, Hedge.Vertex v2)
Face getSharedFace(Vertex v1, Vertex v2) {
  // C#: foreach (var he in v1.GetAdjacentHalfedges())
  for (Halfedge he : v1.outgoingHalfedges()) {
    // C#: if (he.Next.Vertex == v2 || he.Prev.Vertex == v2)
    if (he.next().tipVertex() == v2 || he.next().next().tipVertex() == v2) {
      // C#: var face = he.Face;
      Face face = he.face();

      // C#: if (face != null && !face.IsBoundaryLoop() && !face.IsDeleted)
      if (face != Face() && !face.isBoundaryLoop()) {
        // C#: return face;
        return face;
      }
    }
  }
  // C#: return null;
  return Face();
}

// C#: public static List<Hedge.Face> GetAllSharedFaces(Hedge.Vertex v1, Hedge.Vertex v2)
std::vector<Face> getAllSharedFaces(Vertex v1, Vertex v2) {
  // C#: var result = new List<Hedge.Face>();
  std::vector<Face> result;

  // C#: foreach (var he in v1.GetAdjacentHalfedges())
  for (Halfedge he : v1.outgoingHalfedges()) {
    // C#: if (he.Next.Vertex == v2 || he.Prev.Vertex == v2)
    if (he.next().tipVertex() == v2 || he.next().next().tipVertex() == v2) {
      // C#: var face = he.Face;
      Face face = he.face();

      // C#: if (face != null && !face.IsBoundaryLoop() && !face.IsDeleted && !result.Contains(face))
      if (face != Face() && !face.isBoundaryLoop()) {
        bool found = false;
        for (Face f : result) {
          if (f == face) { found = true; break; }
        }
        if (!found) {
          // C#: result.Add(face);
          result.push_back(face);
        }
      }
    }
  }
  // C#: return result;
  return result;
}

// C#: public static bool AreFacesAdjacent(Hedge.Face f1, Hedge.Face f2)
bool areFacesAdjacent(Face f1, Face f2) {
  // C#: var halfedges = f1.AdjacentHalfedgesCached;
  // C#: for (int i = 0; i < halfedges.Length; i++)
  for (Halfedge he : f1.adjacentHalfedges()) {
    // C#: if (halfedges[i].Twin.Face == f2)
    if (he.twin().face() == f2) {
      // C#: return true;
      return true;
    }
  }
  // C#: return false;
  return false;
}

// C#: public static bool FaceContainsEdge(Hedge.Face face, Hedge.Vertex v1, Hedge.Vertex v2)
bool faceContainsEdge(Face face, Vertex v1, Vertex v2) {
  // C#: var halfedges = face.AdjacentHalfedgesCached;
  // C#: for (int i = 0; i < halfedges.Length; i++)
  for (Halfedge he : face.adjacentHalfedges()) {
    // C#: var he = halfedges[i];
    // C#: if ((he.Vertex == v1 && he.Next.Vertex == v2) ||
    // C#:     (he.Vertex == v2 && he.Next.Vertex == v1))
    Vertex heV = he.tailVertex();
    Vertex heNextV = he.tipVertex();
    if ((heV == v1 && heNextV == v2) || (heV == v2 && heNextV == v1)) {
      // C#: return true;
      return true;
    }
  }
  // C#: return false;
  return false;
}

// C#: public static bool FaceContainsVertex(Hedge.Face face, Hedge.Vertex vertex)
bool faceContainsVertex(Face face, Vertex vertex) {
  // C#: return face.GetAdjacentVertices().Contains(vertex);
  for (Vertex v : face.adjacentVertices()) {
    if (v == vertex) return true;
  }
  return false;
}

// C#: public static Hedge.Halfedge? FindHalfedgeToVertex(Hedge.Face face, Hedge.Vertex vertex)
Halfedge findHalfedgeToVertex(Face face, Vertex vertex) {
  // C#: var halfedges = face.AdjacentHalfedgesCached;
  // C#: for (int i = 0; i < halfedges.Length; i++)
  for (Halfedge he : face.adjacentHalfedges()) {
    // C#: if (halfedges[i].Next.Vertex == vertex)
    // Note: C# he.Next.Vertex is the tip of he, which is he.tipVertex() in geometry-central
    if (he.tipVertex() == vertex) {
      // C#: return halfedges[i];
      return he;
    }
  }
  // C#: return null;
  return Halfedge();
}

// C#: public static Hedge.Halfedge? FindHalfedgeFromVertex(Hedge.Face face, Hedge.Vertex vertex)
Halfedge findHalfedgeFromVertex(Face face, Vertex vertex) {
  // C#: var halfedges = face.AdjacentHalfedgesCached;
  // C#: for (int i = 0; i < halfedges.Length; i++)
  for (Halfedge he : face.adjacentHalfedges()) {
    // C#: if (halfedges[i].Vertex == vertex)
    // Note: C# he.Vertex is the tail of he, which is he.tailVertex() in geometry-central
    if (he.tailVertex() == vertex) {
      // C#: return halfedges[i];
      return he;
    }
  }
  // C#: return null;
  return Halfedge();
}

// ============================================================================
// VeryDiscreteGeodesicExplorerHelper - Line-by-line port
// C# from VeryDiscreteGeodesicExplorerHelper.cs
// ============================================================================

// C#: private static bool SegmentCrossesPortal(Vector2 v0, Vector2 target, Vector2 a, Vector2 b, float eps)
bool segmentCrossesPortal(Vector2 v0, Vector2 target, Vector2 a, Vector2 b, double eps) {
  // C#: // Direction from v0 to target
  // C#: var dx = target.X - v0.X;
  double dx = target.x - v0.x;

  // C#: var dy = target.Y - v0.Y;
  double dy = target.y - v0.y;

  // C#: // Cross products: which side of the line is each portal endpoint?
  // C#: var crossA = dx * (a.Y - v0.Y) - dy * (a.X - v0.X);
  double crossA = dx * (a.y - v0.y) - dy * (a.x - v0.x);

  // C#: var crossB = dx * (b.Y - v0.Y) - dy * (b.X - v0.X);
  double crossB = dx * (b.y - v0.y) - dy * (b.x - v0.x);

  // C#: // Opposite signs = crosses, product near zero = grazes endpoint (also valid)
  // C#: return crossA * crossB <= eps;
  return crossA * crossB <= eps;
}

// C#: internal static Vector2 FlattenVertex(Hedge.Vertex v, Hedge.Halfedge portal, Vector2 flatA, Vector2 flatB)
Vector2 flattenVertex(Vertex v, Halfedge portal, Vector2 flatA, Vector2 flatB,
                      VertexPositionGeometry& geom) {
  // C#: var d1 = v.Position.DistanceTo(portal.Vertex.Position);
  double d1 = norm(geom.vertexPositions[v] - geom.vertexPositions[portal.tailVertex()]);

  // C#: var d2 = v.Position.DistanceTo(portal.Twin.Vertex.Position);
  double d2 = norm(geom.vertexPositions[v] - geom.vertexPositions[portal.twin().tailVertex()]);

  // C#: return V2.ComputeTriangleApex(flatA, flatB, d1, d2, pickPositiveY: false);
  return computeTriangleApex(flatA, flatB, d1, d2, false);
}

// C#: internal static Vector2 GetFlatPosition(...)
Vector2 getFlatPosition(Vertex v, Halfedge portal, Vector2 flatP1, Vector2 flatP2, Vector2 flatApex) {
  // C#: if (v == portal.Vertex) return flatP1;
  if (v == portal.tailVertex()) return flatP1;

  // C#: if (v == portal.Twin.Vertex) return flatP2;
  if (v == portal.twin().tailVertex()) return flatP2;

  // C#: return flatApex;
  return flatApex;
}


// Helper for creating blocked candidates
static CandidateVertex makeBlockedCandidate(CandidateName name, BlockedReason reason) {
  CandidateVertex c;
  c.name = name;
  c.blocked = reason;
  c.isReachable = false;
  return c;
}

static CandidateVertex makeBlockedCandidateWithVertex(CandidateName name, Vertex v, Vector2 flatPos, BlockedReason reason, const std::string& portal = "") {
  CandidateVertex c;
  c.name = name;
  c.vertex = v;
  c.flatPosition = flatPos;
  c.blocked = reason;
  c.blockedPortal = portal;
  c.isReachable = false;
  return c;
}

// C#: internal static CandidateVertex CheckReachability1(...)
CandidateVertex checkReachability1(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, double scaleTolerance) {
  // C#: if (vertex == null || targetFace == null)
  if (vertex == Vertex() || targetFace == Face()) {
    // C#: return new CandidateVertex(name, null, VeryDiscreteGeodesicExplorer.NaN, float.NaN, false, BlockedReason.ParentFaceBoundary);
    return makeBlockedCandidate(name, BlockedReason::ParentFaceBoundary);
  }

  // C#: if (float.IsNaN(flatTarget.X))
  if (std::isnan(flatTarget.x)) {
    // C#: return new CandidateVertex(name, vertex, flatTarget, float.NaN, false, BlockedReason.FlatteningFailed);
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::FlatteningFailed);
  }

  // C#: if (targetFace.IsBoundaryLoop())
  if (targetFace.isBoundaryLoop()) {
    // C#: return new CandidateVertex(name, null, flatTarget, float.NaN, false, BlockedReason.TargetFaceBoundary);
    return makeBlockedCandidate(name, BlockedReason::TargetFaceBoundary);
  }

  // C#: var eps = scaleTolerance * 1e-6f;
  double eps = scaleTolerance * 1e-6;

  // C#: if (!SegmentCrossesPortal(flatV0, flatTarget, p1a, p1b, eps))
  if (!segmentCrossesPortal(flatV0, flatTarget, p1a, p1b, eps)) {
    // C#: return new CandidateVertex(name, vertex, flatTarget, float.NaN, false, BlockedReason.PortalBlocked, "portal1");
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal1");
  }

  // C#: return new CandidateVertex(name, vertex, flatTarget, (flatTarget - flatV0).Length(), true, BlockedReason.None);
  CandidateVertex c;
  c.name = name;
  c.vertex = vertex;
  c.flatPosition = flatTarget;
  c.distance = norm(flatTarget - flatV0);
  c.isReachable = true;
  c.blocked = BlockedReason::None_;
  return c;
}

// C#: internal static CandidateVertex CheckReachability2(...)
CandidateVertex checkReachability2(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, Vector2 p2a, Vector2 p2b, double scaleTolerance) {
  // Same pattern as checkReachability1 with 2 portals
  if (vertex == Vertex() || targetFace == Face()) {
    return makeBlockedCandidate(name, BlockedReason::ParentFaceBoundary);
  }
  if (std::isnan(flatTarget.x)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::FlatteningFailed);
  }
  if (targetFace.isBoundaryLoop()) {
    return makeBlockedCandidate(name, BlockedReason::TargetFaceBoundary);
  }

  double eps = scaleTolerance * 1e-6;
  if (!segmentCrossesPortal(flatV0, flatTarget, p1a, p1b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal1");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p2a, p2b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal2");
  }

  CandidateVertex c;
  c.name = name;
  c.vertex = vertex;
  c.flatPosition = flatTarget;
  c.distance = norm(flatTarget - flatV0);
  c.isReachable = true;
  c.blocked = BlockedReason::None_;
  return c;
}

// C#: internal static CandidateVertex CheckReachability3(...)
CandidateVertex checkReachability3(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, Vector2 p2a, Vector2 p2b, Vector2 p3a, Vector2 p3b,
                                   double scaleTolerance) {
  if (vertex == Vertex() || targetFace == Face()) {
    return makeBlockedCandidate(name, BlockedReason::ParentFaceBoundary);
  }
  if (std::isnan(flatTarget.x)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::FlatteningFailed);
  }
  if (targetFace.isBoundaryLoop()) {
    return makeBlockedCandidate(name, BlockedReason::TargetFaceBoundary);
  }

  double eps = scaleTolerance * 1e-6;
  if (!segmentCrossesPortal(flatV0, flatTarget, p1a, p1b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal1");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p2a, p2b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal2");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p3a, p3b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal3");
  }

  CandidateVertex c;
  c.name = name;
  c.vertex = vertex;
  c.flatPosition = flatTarget;
  c.distance = norm(flatTarget - flatV0);
  c.isReachable = true;
  c.blocked = BlockedReason::None_;
  return c;
}

// C#: internal static CandidateVertex CheckReachability4(...)
CandidateVertex checkReachability4(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, Vector2 p2a, Vector2 p2b, Vector2 p3a, Vector2 p3b,
                                   Vector2 p4a, Vector2 p4b, double scaleTolerance) {
  if (vertex == Vertex() || targetFace == Face()) {
    return makeBlockedCandidate(name, BlockedReason::ParentFaceBoundary);
  }
  if (std::isnan(flatTarget.x)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::FlatteningFailed);
  }
  if (targetFace.isBoundaryLoop()) {
    return makeBlockedCandidate(name, BlockedReason::TargetFaceBoundary);
  }

  double eps = scaleTolerance * 1e-6;
  if (!segmentCrossesPortal(flatV0, flatTarget, p1a, p1b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal1");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p2a, p2b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal2");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p3a, p3b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal3");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p4a, p4b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal4");
  }

  CandidateVertex c;
  c.name = name;
  c.vertex = vertex;
  c.flatPosition = flatTarget;
  c.distance = norm(flatTarget - flatV0);
  c.isReachable = true;
  c.blocked = BlockedReason::None_;
  return c;
}

// C#: internal static CandidateVertex CheckReachability5(...)
CandidateVertex checkReachability5(CandidateName name, Vertex vertex, Vector2 flatV0, Vector2 flatTarget, Face targetFace,
                                   Vector2 p1a, Vector2 p1b, Vector2 p2a, Vector2 p2b, Vector2 p3a, Vector2 p3b,
                                   Vector2 p4a, Vector2 p4b, Vector2 p5a, Vector2 p5b, double scaleTolerance) {
  if (vertex == Vertex() || targetFace == Face()) {
    return makeBlockedCandidate(name, BlockedReason::ParentFaceBoundary);
  }
  if (std::isnan(flatTarget.x)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::FlatteningFailed);
  }
  if (targetFace.isBoundaryLoop()) {
    return makeBlockedCandidate(name, BlockedReason::TargetFaceBoundary);
  }

  double eps = scaleTolerance * 1e-6;
  if (!segmentCrossesPortal(flatV0, flatTarget, p1a, p1b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal1");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p2a, p2b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal2");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p3a, p3b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal3");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p4a, p4b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal4");
  }
  if (!segmentCrossesPortal(flatV0, flatTarget, p5a, p5b, eps)) {
    return makeBlockedCandidateWithVertex(name, vertex, flatTarget, BlockedReason::PortalBlocked, "portal5");
  }

  CandidateVertex c;
  c.name = name;
  c.vertex = vertex;
  c.flatPosition = flatTarget;
  c.distance = norm(flatTarget - flatV0);
  c.isReachable = true;
  c.blocked = BlockedReason::None_;
  return c;
}

// ============================================================================
// VeryDiscreteGeodesicExplorer - L5 exploration (fixed depth)
// ============================================================================

ExplorationResult explore(Corner corner, VertexPositionGeometry& geom) {
  // Fixed L5 exploration depth
  constexpr bool computeL2 = true;
  constexpr bool computeL3 = true;
  constexpr bool computeL4 = true;
  constexpr bool computeL5 = true;

  Vertex v0 = corner.vertex();
  Face f0 = corner.face();
  Face fL1Raw = corner.halfedge().twin().face();

  // Create result
  ExplorationResult result;
  result.corner = corner;

  // NaN vector helper
  Vector2 NaN = {std::numeric_limits<double>::quiet_NaN(),
                 std::numeric_limits<double>::quiet_NaN()};

  // C#: if (fL1Raw.IsBoundaryLoop())
  if (fL1Raw.isBoundaryLoop()) {
    // All candidates blocked with F1Boundary
    result.L1 = makeBlockedCandidate(CandidateName::L1, BlockedReason::F1Boundary);
    result.L2L = makeBlockedCandidate(CandidateName::L2L, BlockedReason::F1Boundary);
    result.L2R = makeBlockedCandidate(CandidateName::L2R, BlockedReason::F1Boundary);
    result.L3L = makeBlockedCandidate(CandidateName::L3L, BlockedReason::F1Boundary);
    result.L3R = makeBlockedCandidate(CandidateName::L3R, BlockedReason::F1Boundary);
    result.L4L = makeBlockedCandidate(CandidateName::L4L, BlockedReason::F1Boundary);
    result.L4R = makeBlockedCandidate(CandidateName::L4R, BlockedReason::F1Boundary);
    result.L5L = makeBlockedCandidate(CandidateName::L5L, BlockedReason::F1Boundary);
    result.L5R = makeBlockedCandidate(CandidateName::L5R, BlockedReason::F1Boundary);
    result.L5LM = makeBlockedCandidate(CandidateName::L5LM, BlockedReason::F1Boundary);
    result.L5RM = makeBlockedCandidate(CandidateName::L5RM, BlockedReason::F1Boundary);
    return result;
  }

  // C#: var fL1 = fL1Raw;
  Face fL1 = fL1Raw;
  // C#: var heL1 = corner.Halfedge.Twin;
  Halfedge heL1 = corner.halfedge().twin();
  // C#: var vL1 = heL1.Prev.Vertex;
  // Note: heL1.Prev is the halfedge before heL1 in the face, and .Vertex is its tail
  // In geometry-central: heL1.next().next() goes around the triangle, or use the apex
  Vertex vL1 = heL1.next().next().vertex();  // The apex opposite to heL1

  // C#: var heL1_v1_3d = heL1.Vertex.Position;
  Vector3 heL1_v1_3d = geom.vertexPositions[heL1.tailVertex()];
  // C#: var heL1_v2_3d = heL1.Twin.Vertex.Position;
  Vector3 heL1_v2_3d = geom.vertexPositions[heL1.twin().tailVertex()];
  // C#: var portalLen = heL1_v1_3d.DistanceTo(heL1_v2_3d);
  double portalLen = norm(heL1_v1_3d - heL1_v2_3d);

  // C#: var flat_heL1_a = new Vector2(-portalLen / 2, 0);
  Vector2 flat_heL1_a = {-portalLen / 2.0, 0.0};
  // C#: var flat_heL1_b = new Vector2(portalLen / 2, 0);
  Vector2 flat_heL1_b = {portalLen / 2.0, 0.0};

  // C#: var flat_v0 = V2.ComputeTriangleApex(flat_heL1_a, flat_heL1_b,
  //         v0.Position.DistanceTo(heL1_v1_3d), v0.Position.DistanceTo(heL1_v2_3d), pickPositiveY: true);
  double d_v0_to_v1 = norm(geom.vertexPositions[v0] - heL1_v1_3d);
  double d_v0_to_v2 = norm(geom.vertexPositions[v0] - heL1_v2_3d);
  Vector2 flat_v0 = computeTriangleApex(flat_heL1_a, flat_heL1_b, d_v0_to_v1, d_v0_to_v2, true);

  if (std::isnan(flat_v0.x)) {
    // Return blocked result with V0FlatteningFailed
    result.L1 = makeBlockedCandidate(CandidateName::L1, BlockedReason::V0FlatteningFailed);
    result.L2L = makeBlockedCandidate(CandidateName::L2L, BlockedReason::V0FlatteningFailed);
    result.L2R = makeBlockedCandidate(CandidateName::L2R, BlockedReason::V0FlatteningFailed);
    result.L3L = makeBlockedCandidate(CandidateName::L3L, BlockedReason::V0FlatteningFailed);
    result.L3R = makeBlockedCandidate(CandidateName::L3R, BlockedReason::V0FlatteningFailed);
    result.L4L = makeBlockedCandidate(CandidateName::L4L, BlockedReason::V0FlatteningFailed);
    result.L4R = makeBlockedCandidate(CandidateName::L4R, BlockedReason::V0FlatteningFailed);
    result.L5L = makeBlockedCandidate(CandidateName::L5L, BlockedReason::V0FlatteningFailed);
    result.L5R = makeBlockedCandidate(CandidateName::L5R, BlockedReason::V0FlatteningFailed);
    result.L5LM = makeBlockedCandidate(CandidateName::L5LM, BlockedReason::V0FlatteningFailed);
    result.L5RM = makeBlockedCandidate(CandidateName::L5RM, BlockedReason::V0FlatteningFailed);
    return result;
  }

  // Compute L1 position
  bool computedL1 = true;
  double d_vL1_to_v1 = norm(geom.vertexPositions[vL1] - heL1_v1_3d);
  double d_vL1_to_v2 = norm(geom.vertexPositions[vL1] - heL1_v2_3d);
  Vector2 flat_L1 = computeTriangleApex(flat_heL1_a, flat_heL1_b, d_vL1_to_v1, d_vL1_to_v2, false);

  // ========== LEFT SIDE ==========
  Halfedge heL2L, heL3L, heL4L, heL5L, heL4LM, heL5LM;
  Face fL2L, fL3L, fL4L, fL5L, fL4LM, fL5LM;
  Vertex vL2L, vL3L, vL4L, vL5L, vL5LM;
  bool fL2LIsBoundary = true, fL3LIsBoundary = true, fL4LIsBoundary = true, fL5LIsBoundary = true;
  bool fL4LMIsBoundary = true, fL5LMIsBoundary = true;
  Vector2 flat_heL2L_a = NaN, flat_heL2L_b = NaN, flat_L2L = NaN;
  Vector2 flat_heL3L_a = NaN, flat_heL3L_b = NaN, flat_L3L = NaN;
  Vector2 flat_heL4L_a = NaN, flat_heL4L_b = NaN, flat_L4L = NaN;
  Vector2 flat_heL5L_a = NaN, flat_heL5L_b = NaN, flat_L5L = NaN;
  Vector2 flat_heL4LM_a = NaN, flat_heL4LM_b = NaN;
  Vector2 flat_heL5LM_a = NaN, flat_heL5LM_b = NaN, flat_L5LM = NaN;
  bool computedL2L = false, computedL3L = false, computedL4L = false, computedL5L = false;
  bool computedL5LM = false;

  // C#: if (computeL2)
  if (computeL2) {
    // C#: heL2L = heL1.Next.Twin; fL2L = heL2L.Face; fL2LIsBoundary = fL2L.IsBoundaryLoop(); vL2L = heL2L.Prev.Vertex;
    heL2L = heL1.next().twin();
    fL2L = heL2L.face();
    fL2LIsBoundary = fL2L.isBoundaryLoop();
    vL2L = heL2L.next().next().vertex();  // apex opposite to heL2L

    computedL2L = true;
    {
      // C#: flat_heL2L_a = VeryDiscreteGeodesicExplorerHelper.GetFlatPosition(heL2L.Vertex, heL1, flat_heL1_a, flat_heL1_b, flat_L1);
      flat_heL2L_a = getFlatPosition(heL2L.tailVertex(), heL1, flat_heL1_a, flat_heL1_b, flat_L1);
      // C#: flat_heL2L_b = VeryDiscreteGeodesicExplorerHelper.GetFlatPosition(heL2L.Twin.Vertex, heL1, flat_heL1_a, flat_heL1_b, flat_L1);
      flat_heL2L_b = getFlatPosition(heL2L.twin().tailVertex(), heL1, flat_heL1_a, flat_heL1_b, flat_L1);
      // C#: flat_L2L = VeryDiscreteGeodesicExplorerHelper.FlattenVertex(vL2L, heL2L, flat_heL2L_a, flat_heL2L_b);
      flat_L2L = flattenVertex(vL2L, heL2L, flat_heL2L_a, flat_heL2L_b, geom);

      // Continue to L3L if enabled and not boundary
      if (computeL3 && !fL2LIsBoundary) {
        // C#: heL3L = heL2L.Prev.Twin; fL3L = heL3L.Face; fL3LIsBoundary = fL3L.IsBoundaryLoop(); vL3L = heL3L.Prev.Vertex;
        heL3L = heL2L.next().next().twin();  // Prev = Next.Next in triangle
        fL3L = heL3L.face();
        fL3LIsBoundary = fL3L.isBoundaryLoop();
        vL3L = heL3L.next().next().vertex();

        computedL3L = true;
        {
          flat_heL3L_a = getFlatPosition(heL3L.tailVertex(), heL2L, flat_heL2L_a, flat_heL2L_b, flat_L2L);
          flat_heL3L_b = getFlatPosition(heL3L.twin().tailVertex(), heL2L, flat_heL2L_a, flat_heL2L_b, flat_L2L);
          flat_L3L = flattenVertex(vL3L, heL3L, flat_heL3L_a, flat_heL3L_b, geom);

          if (computeL4 && !fL3LIsBoundary) {
            // C#: heL4L = heL3L.Next.Twin;
            heL4L = heL3L.next().twin();
            fL4L = heL4L.face();
            fL4LIsBoundary = fL4L.isBoundaryLoop();
            vL4L = heL4L.next().next().vertex();

            computedL4L = true;
            {
              flat_heL4L_a = getFlatPosition(heL4L.tailVertex(), heL3L, flat_heL3L_a, flat_heL3L_b, flat_L3L);
              flat_heL4L_b = getFlatPosition(heL4L.twin().tailVertex(), heL3L, flat_heL3L_a, flat_heL3L_b, flat_L3L);
              flat_L4L = flattenVertex(vL4L, heL4L, flat_heL4L_a, flat_heL4L_b, geom);

              if (computeL5 && !fL4LIsBoundary) {
                // C#: heL5L = heL4L.Prev.Twin;
                heL5L = heL4L.next().next().twin();
                fL5L = heL5L.face();
                fL5LIsBoundary = fL5L.isBoundaryLoop();
                vL5L = heL5L.next().next().vertex();

                computedL5L = true;
                flat_heL5L_a = getFlatPosition(heL5L.tailVertex(), heL4L, flat_heL4L_a, flat_heL4L_b, flat_L4L);
                flat_heL5L_b = getFlatPosition(heL5L.twin().tailVertex(), heL4L, flat_heL4L_a, flat_heL4L_b, flat_L4L);
                flat_L5L = flattenVertex(vL5L, heL5L, flat_heL5L_a, flat_heL5L_b, geom);
              }
            }
          }

          // L5LM: Middle path from left (branches via Prev.Twin from L3)
          if (computeL5 && !fL3LIsBoundary) {
            // C#: heL4LM = heL3L.Prev.Twin;
            heL4LM = heL3L.next().next().twin();  // Different branch than L4L
            fL4LM = heL4LM.face();
            fL4LMIsBoundary = fL4LM.isBoundaryLoop();
            if (!fL4LMIsBoundary) {
              flat_heL4LM_a = getFlatPosition(heL4LM.tailVertex(), heL3L, flat_heL3L_a, flat_heL3L_b, flat_L3L);
              flat_heL4LM_b = getFlatPosition(heL4LM.twin().tailVertex(), heL3L, flat_heL3L_a, flat_heL3L_b, flat_L3L);

              // C#: heL5LM = heL4LM.Next.Twin;
              heL5LM = heL4LM.next().twin();
              fL5LM = heL5LM.face();
              fL5LMIsBoundary = fL5LM.isBoundaryLoop();
              // C#: vL5LM = heL5LM.ApexVertex;
              vL5LM = heL5LM.next().next().vertex();

              computedL5LM = true;
              {
                Vertex heL4LMApex = heL4LM.next().next().vertex();
                Vector2 flatL4LMApex = flattenVertex(heL4LMApex, heL4LM, flat_heL4LM_a, flat_heL4LM_b, geom);
                flat_heL5LM_a = getFlatPosition(heL5LM.tailVertex(), heL4LM, flat_heL4LM_a, flat_heL4LM_b, flatL4LMApex);
                flat_heL5LM_b = getFlatPosition(heL5LM.twin().tailVertex(), heL4LM, flat_heL4LM_a, flat_heL4LM_b, flatL4LMApex);
                flat_L5LM = flattenVertex(vL5LM, heL5LM, flat_heL5LM_a, flat_heL5LM_b, geom);
              }
            }
          }
        }
      }
    }
  }

  // ========== RIGHT SIDE ==========
  // Similar structure to left side
  // ========== RIGHT SIDE ==========
  Halfedge heL2R, heL3R, heL4R, heL5R, heL4RM, heL5RM;
  Face fL2R, fL3R, fL4R, fL5R, fL4RM, fL5RM;
  Vertex vL2R, vL3R, vL4R, vL5R, vL5RM;
  bool fL2RIsBoundary = true, fL3RIsBoundary = true, fL4RIsBoundary = true, fL5RIsBoundary = true;
  bool fL4RMIsBoundary = true, fL5RMIsBoundary = true;
  Vector2 flat_heL2R_a = NaN, flat_heL2R_b = NaN, flat_L2R = NaN;
  Vector2 flat_heL3R_a = NaN, flat_heL3R_b = NaN, flat_L3R = NaN;
  Vector2 flat_heL4R_a = NaN, flat_heL4R_b = NaN, flat_L4R = NaN;
  Vector2 flat_heL5R_a = NaN, flat_heL5R_b = NaN, flat_L5R = NaN;
  Vector2 flat_heL4RM_a = NaN, flat_heL4RM_b = NaN;
  Vector2 flat_heL5RM_a = NaN, flat_heL5RM_b = NaN, flat_L5RM = NaN;
  bool computedL2R = false, computedL3R = false, computedL4R = false, computedL5R = false;
  bool computedL5RM = false;

  if (computeL2) {
    // C#: heL2R = heL1.Prev.Twin;  (Prev = Next.Next in triangle)
    heL2R = heL1.next().next().twin();
    fL2R = heL2R.face();
    fL2RIsBoundary = fL2R.isBoundaryLoop();
    vL2R = heL2R.next().next().vertex();

    computedL2R = true;
    {
      flat_heL2R_a = getFlatPosition(heL2R.tailVertex(), heL1, flat_heL1_a, flat_heL1_b, flat_L1);
      flat_heL2R_b = getFlatPosition(heL2R.twin().tailVertex(), heL1, flat_heL1_a, flat_heL1_b, flat_L1);
      flat_L2R = flattenVertex(vL2R, heL2R, flat_heL2R_a, flat_heL2R_b, geom);

      if (computeL3 && !fL2RIsBoundary) {
        // C#: heL3R = heL2R.Next.Twin;
        heL3R = heL2R.next().twin();
        fL3R = heL3R.face();
        fL3RIsBoundary = fL3R.isBoundaryLoop();
        vL3R = heL3R.next().next().vertex();

        computedL3R = true;
        {
          flat_heL3R_a = getFlatPosition(heL3R.tailVertex(), heL2R, flat_heL2R_a, flat_heL2R_b, flat_L2R);
          flat_heL3R_b = getFlatPosition(heL3R.twin().tailVertex(), heL2R, flat_heL2R_a, flat_heL2R_b, flat_L2R);
          flat_L3R = flattenVertex(vL3R, heL3R, flat_heL3R_a, flat_heL3R_b, geom);

          if (computeL4 && !fL3RIsBoundary) {
            // C#: heL4R = heL3R.Prev.Twin;
            heL4R = heL3R.next().next().twin();
            fL4R = heL4R.face();
            fL4RIsBoundary = fL4R.isBoundaryLoop();
            vL4R = heL4R.next().next().vertex();

            computedL4R = true;
            {
              flat_heL4R_a = getFlatPosition(heL4R.tailVertex(), heL3R, flat_heL3R_a, flat_heL3R_b, flat_L3R);
              flat_heL4R_b = getFlatPosition(heL4R.twin().tailVertex(), heL3R, flat_heL3R_a, flat_heL3R_b, flat_L3R);
              flat_L4R = flattenVertex(vL4R, heL4R, flat_heL4R_a, flat_heL4R_b, geom);

              if (computeL5 && !fL4RIsBoundary) {
                // C#: heL5R = heL4R.Next.Twin;
                heL5R = heL4R.next().twin();
                fL5R = heL5R.face();
                fL5RIsBoundary = fL5R.isBoundaryLoop();
                vL5R = heL5R.next().next().vertex();

                computedL5R = true;
                flat_heL5R_a = getFlatPosition(heL5R.tailVertex(), heL4R, flat_heL4R_a, flat_heL4R_b, flat_L4R);
                flat_heL5R_b = getFlatPosition(heL5R.twin().tailVertex(), heL4R, flat_heL4R_a, flat_heL4R_b, flat_L4R);
                flat_L5R = flattenVertex(vL5R, heL5R, flat_heL5R_a, flat_heL5R_b, geom);
              }
            }
          }

          // L5RM: Middle path from right
          if (computeL5 && !fL3RIsBoundary) {
            // C#: heL4RM = heL3R.Next.Twin;
            heL4RM = heL3R.next().twin();
            fL4RM = heL4RM.face();
            fL4RMIsBoundary = fL4RM.isBoundaryLoop();
            if (!fL4RMIsBoundary) {
              flat_heL4RM_a = getFlatPosition(heL4RM.tailVertex(), heL3R, flat_heL3R_a, flat_heL3R_b, flat_L3R);
              flat_heL4RM_b = getFlatPosition(heL4RM.twin().tailVertex(), heL3R, flat_heL3R_a, flat_heL3R_b, flat_L3R);

              // C#: heL5RM = heL4RM.Prev.Twin;
              heL5RM = heL4RM.next().next().twin();
              fL5RM = heL5RM.face();
              fL5RMIsBoundary = fL5RM.isBoundaryLoop();
              vL5RM = heL5RM.next().next().vertex();

              computedL5RM = true;
              {
                Vertex heL4RMApex = heL4RM.next().next().vertex();
                Vector2 flatL4RMApex = flattenVertex(heL4RMApex, heL4RM, flat_heL4RM_a, flat_heL4RM_b, geom);
                flat_heL5RM_a = getFlatPosition(heL5RM.tailVertex(), heL4RM, flat_heL4RM_a, flat_heL4RM_b, flatL4RMApex);
                flat_heL5RM_b = getFlatPosition(heL5RM.twin().tailVertex(), heL4RM, flat_heL4RM_a, flat_heL4RM_b, flatL4RMApex);
                flat_L5RM = flattenVertex(vL5RM, heL5RM, flat_heL5RM_a, flat_heL5RM_b, geom);
              }
            }
          }
        }
      }
    }
  }

  // Check reachability for all candidates
  double portalLen1 = portalLen;

  // C#: var candidate_L1 = computedL1 ? ... : NotComputedL1;
  result.L1 = computedL1
    ? checkReachability1(CandidateName::L1, vL1, flat_v0, flat_L1, fL1, flat_heL1_a, flat_heL1_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L1, BlockedReason::NotComputed);

  result.L2L = computedL2L
    ? checkReachability2(CandidateName::L2L, vL2L, flat_v0, flat_L2L, fL2L, flat_heL1_a, flat_heL1_b, flat_heL2L_a, flat_heL2L_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L2L, BlockedReason::NotComputed);

  result.L3L = computedL3L
    ? checkReachability3(CandidateName::L3L, vL3L, flat_v0, flat_L3L, fL3L, flat_heL1_a, flat_heL1_b, flat_heL2L_a, flat_heL2L_b, flat_heL3L_a, flat_heL3L_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L3L, BlockedReason::NotComputed);

  result.L4L = computedL4L
    ? checkReachability4(CandidateName::L4L, vL4L, flat_v0, flat_L4L, fL4L, flat_heL1_a, flat_heL1_b, flat_heL2L_a, flat_heL2L_b, flat_heL3L_a, flat_heL3L_b, flat_heL4L_a, flat_heL4L_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L4L, BlockedReason::NotComputed);

  result.L5L = computedL5L
    ? checkReachability5(CandidateName::L5L, vL5L, flat_v0, flat_L5L, fL5L, flat_heL1_a, flat_heL1_b, flat_heL2L_a, flat_heL2L_b, flat_heL3L_a, flat_heL3L_b, flat_heL4L_a, flat_heL4L_b, flat_heL5L_a, flat_heL5L_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L5L, BlockedReason::NotComputed);

  result.L5LM = computedL5LM
    ? checkReachability5(CandidateName::L5LM, vL5LM, flat_v0, flat_L5LM, fL5LM, flat_heL1_a, flat_heL1_b, flat_heL2L_a, flat_heL2L_b, flat_heL3L_a, flat_heL3L_b, flat_heL4LM_a, flat_heL4LM_b, flat_heL5LM_a, flat_heL5LM_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L5LM, BlockedReason::NotComputed);

  result.L2R = computedL2R
    ? checkReachability2(CandidateName::L2R, vL2R, flat_v0, flat_L2R, fL2R, flat_heL1_a, flat_heL1_b, flat_heL2R_a, flat_heL2R_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L2R, BlockedReason::NotComputed);

  result.L3R = computedL3R
    ? checkReachability3(CandidateName::L3R, vL3R, flat_v0, flat_L3R, fL3R, flat_heL1_a, flat_heL1_b, flat_heL2R_a, flat_heL2R_b, flat_heL3R_a, flat_heL3R_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L3R, BlockedReason::NotComputed);

  result.L4R = computedL4R
    ? checkReachability4(CandidateName::L4R, vL4R, flat_v0, flat_L4R, fL4R, flat_heL1_a, flat_heL1_b, flat_heL2R_a, flat_heL2R_b, flat_heL3R_a, flat_heL3R_b, flat_heL4R_a, flat_heL4R_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L4R, BlockedReason::NotComputed);

  result.L5R = computedL5R
    ? checkReachability5(CandidateName::L5R, vL5R, flat_v0, flat_L5R, fL5R, flat_heL1_a, flat_heL1_b, flat_heL2R_a, flat_heL2R_b, flat_heL3R_a, flat_heL3R_b, flat_heL4R_a, flat_heL4R_b, flat_heL5R_a, flat_heL5R_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L5R, BlockedReason::NotComputed);

  result.L5RM = computedL5RM
    ? checkReachability5(CandidateName::L5RM, vL5RM, flat_v0, flat_L5RM, fL5RM, flat_heL1_a, flat_heL1_b, flat_heL2R_a, flat_heL2R_b, flat_heL3R_a, flat_heL3R_b, flat_heL4RM_a, flat_heL4RM_b, flat_heL5RM_a, flat_heL5RM_b, portalLen1)
    : makeBlockedCandidate(CandidateName::L5RM, BlockedReason::NotComputed);

  return result;
}

// ============================================================================
// FaceStripWalker - Line-by-line port
// C# from FaceStripWalker.cs
// ============================================================================

// C#: public static WalkDirection DetermineWalkDirection(Hedge.Vertex prev, Hedge.Vertex current, Hedge.Vertex next)
WalkDirection determineWalkDirection(Vertex prev, Vertex current, Vertex next,
                                     VertexPositionGeometry& geom) {
  // C#: var reference = prev.Position - current.Position;
  Vector3 reference = geom.vertexPositions[prev] - geom.vertexPositions[current];

  // C#: var outDir = next.Position - current.Position;
  Vector3 outDir = geom.vertexPositions[next] - geom.vertexPositions[current];

  // C#: var signedAngle = outDir.SignedAngleTo(reference, current.VertexNormal);
  // Compute vertex normal (average of adjacent face normals)
  Vector3 normal = {0, 0, 0};
  for (Face f : current.adjacentFaces()) {
    normal += geom.faceNormal(f);
  }
  normal = unit(normal);

  // Signed angle computation
  Vector3 crossProd = cross(outDir, reference);
  double sinAngle = dot(crossProd, normal);

  // C#: var isLeftTurn = signedAngle < 0;
  bool isLeftTurn = sinAngle < 0;

  // C#: return isLeftTurn ? WalkDirection.CounterClockwise : WalkDirection.Clockwise;
  return isLeftTurn ? WalkDirection::CounterClockwise : WalkDirection::Clockwise;
}

// C#: public static WalkResult WalkToOutgoingEdge(...)
WalkResult walkToOutgoingEdge(Face startFace, Vertex vertex, Vertex targetVertex,
                               WalkDirection direction) {
  // C#: var faces = new List<Hedge.Face>();
  WalkResult result;

  // C#: var startHe = FaceStripUtils.FindHalfedgeToVertex(startFace, vertex);
  Halfedge startHe = findHalfedgeToVertex(startFace, vertex);

  // C#: if (startHe == null)
  if (startHe == Halfedge()) {
    result.error = "No halfedge in face points to vertex";
    return result;
  }

  // C#: if (FaceStripUtils.FaceContainsEdge(startFace, vertex, targetVertex))
  if (faceContainsEdge(startFace, vertex, targetVertex)) {
    result.finalFace = startFace;
    result.reachedTarget = true;
    return result;
  }

  // C#: var currentHe = direction == WalkDirection.Clockwise ? startHe.Next.Twin : startHe.Twin;
  Halfedge currentHe = (direction == WalkDirection::Clockwise)
    ? startHe.next().twin()
    : startHe.twin();

  // C#: const int maxSteps = 100;
  const int maxSteps = 100;

  // C#: for (int i = 0; i < maxSteps; i++)
  for (int i = 0; i < maxSteps; i++) {
    // C#: if (currentHe.IsOnBoundary || currentHe.Face == null || currentHe.Face.IsBoundaryLoop())
    if (currentHe.isInterior() == false || currentHe.face().isBoundaryLoop()) {
      result.error = "Hit boundary while walking around vertex";
      return result;
    }

    // C#: var currentFace = currentHe.Face;
    Face currentFace = currentHe.face();

    // C#: if (currentFace == startFace)
    if (currentFace == startFace) {
      result.error = "Wrapped back to start face without finding edge to target";
      return result;
    }

    // C#: faces.Add(currentFace);
    result.faces.push_back(currentFace);

    // C#: if (FaceStripUtils.FaceContainsEdge(currentFace, vertex, targetVertex))
    if (faceContainsEdge(currentFace, vertex, targetVertex)) {
      result.finalFace = currentFace;
      result.reachedTarget = true;
      return result;
    }

    // C#: if (direction == WalkDirection.Clockwise)
    if (direction == WalkDirection::Clockwise) {
      // C#: var leaving = FaceStripUtils.FindHalfedgeFromVertex(currentFace, vertex);
      Halfedge leaving = findHalfedgeFromVertex(currentFace, vertex);
      if (leaving == Halfedge()) {
        result.error = "No halfedge in face leaves from vertex";
        return result;
      }
      // C#: currentHe = leaving.Twin;
      currentHe = leaving.twin();
    } else {
      // C#: var entering = FaceStripUtils.FindHalfedgeToVertex(currentFace, vertex);
      Halfedge entering = findHalfedgeToVertex(currentFace, vertex);
      if (entering == Halfedge()) {
        result.error = "No halfedge in face points to vertex";
        return result;
      }
      // C#: currentHe = entering.Twin;
      currentHe = entering.twin();
    }
  }

  result.error = "Exceeded max steps while walking around vertex";
  return result;
}

// C#: private static WalkResult WalkToFace(...)
WalkResult walkToFace(Face startFace, Face targetFace, Vertex vertex, WalkDirection direction) {
  WalkResult result;

  // C#: if (startFace == targetFace)
  if (startFace == targetFace) {
    result.finalFace = targetFace;
    result.reachedTarget = true;
    return result;
  }

  // C#: var startHe = FaceStripUtils.FindHalfedgeToVertex(startFace, vertex);
  Halfedge startHe = findHalfedgeToVertex(startFace, vertex);
  if (startHe == Halfedge()) {
    result.error = "No starting halfedge";
    return result;
  }

  // C#: var currentHe = direction == WalkDirection.Clockwise ? startHe.Next.Twin : startHe.Twin;
  Halfedge currentHe = (direction == WalkDirection::Clockwise)
    ? startHe.next().twin()
    : startHe.twin();

  const int maxSteps = 100;
  for (int i = 0; i < maxSteps; i++) {
    if (currentHe.isInterior() == false || currentHe.face().isBoundaryLoop()) {
      result.error = "Hit boundary";
      return result;
    }

    Face currentFace = currentHe.face();

    if (currentFace == startFace) {
      result.error = "Wrapped around";
      return result;
    }

    result.faces.push_back(currentFace);

    if (currentFace == targetFace) {
      result.finalFace = targetFace;
      result.reachedTarget = true;
      return result;
    }

    if (direction == WalkDirection::Clockwise) {
      Halfedge leaving = findHalfedgeFromVertex(currentFace, vertex);
      if (leaving == Halfedge()) {
        result.error = "No leaving halfedge";
        return result;
      }
      currentHe = leaving.twin();
    } else {
      Halfedge entering = findHalfedgeToVertex(currentFace, vertex);
      if (entering == Halfedge()) {
        result.error = "No entering halfedge";
        return result;
      }
      currentHe = entering.twin();
    }
  }

  result.error = "Max steps exceeded";
  return result;
}

// C#: public static Hedge.Face? SelectFirstFace(...)
Face selectFirstFace(Vertex v0, Vertex v1, Vertex v2,
                     VertexPositionGeometry& geom,
                     Face targetFace) {
  // C#: var facesOnFirstEdge = FaceStripUtils.GetAllSharedFaces(v0, v1);
  std::vector<Face> facesOnFirstEdge = getAllSharedFaces(v0, v1);

  // C#: if (facesOnFirstEdge.Count == 0) return null;
  if (facesOnFirstEdge.empty()) return Face();

  // C#: if (facesOnFirstEdge.Count == 1 || v2 == null) return facesOnFirstEdge[0];
  if (facesOnFirstEdge.size() == 1 || v2 == Vertex()) return facesOnFirstEdge[0];

  // C#: var turnDirection = DetermineWalkDirection(v0, v1, v2);
  WalkDirection turnDirection = determineWalkDirection(v0, v1, v2, geom);

  // C#: var oppositeDir = turnDirection == WalkDirection.Clockwise ? WalkDirection.CounterClockwise : WalkDirection.Clockwise;
  WalkDirection oppositeDir = (turnDirection == WalkDirection::Clockwise)
    ? WalkDirection::CounterClockwise
    : WalkDirection::Clockwise;

  Face bestFace;
  int bestScore = std::numeric_limits<int>::max();

  // C#: foreach (var face in facesOnFirstEdge)
  for (Face face : facesOnFirstEdge) {
    // C#: foreach (var dir in new[] { turnDirection, oppositeDir })
    for (WalkDirection dir : {turnDirection, oppositeDir}) {
      WalkResult testWalk;

      if (targetFace != Face()) {
        testWalk = walkToFace(face, targetFace, v1, dir);
      } else {
        testWalk = walkToOutgoingEdge(face, v1, v2, dir);
      }

      if (testWalk.reachedTarget) {
        // C#: int entryVertexFaces = testWalk.Faces.Count(f => FaceStripUtils.FaceContainsVertex(f, v0));
        int entryVertexFaces = 0;
        for (Face f : testWalk.faces) {
          if (faceContainsVertex(f, v0)) entryVertexFaces++;
        }

        if (entryVertexFaces < bestScore) {
          bestScore = entryVertexFaces;
          bestFace = face;

          if (entryVertexFaces == 0) return face;
        }
        break;
      }
    }
  }

  return (bestFace != Face()) ? bestFace : facesOnFirstEdge[0];
}

// ============================================================================
// VeryDiscreteGeodesicPathfinder - Line-by-line port (partial - main A* loop)
// C# from VeryDiscreteGeodesicPathfinder.cs
// ============================================================================

// ============================================================================
// GetCrossedFaces helper functions - ported line-by-line from C#
// ============================================================================

// C#: private static Hedge.Face[]? GetFacesL2L(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL2L(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2L = heL1.Next.Twin;
  Halfedge heL2L = heL1.next().twin();
  // C#: if (heL2L.Face.IsBoundaryLoop()) return null;
  if (heL2L.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2L.Face];
  return {f0, fL1, heL2L.face()};
}

// C#: private static Hedge.Face[]? GetFacesL2R(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL2R(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2R = heL1.Prev.Twin;
  Halfedge heL2R = heL1.prevOrbitFace().twin();
  // C#: if (heL2R.Face.IsBoundaryLoop()) return null;
  if (heL2R.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2R.Face];
  return {f0, fL1, heL2R.face()};
}

// C#: private static Hedge.Face[]? GetFacesL3L(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL3L(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2L = heL1.Next.Twin;
  Halfedge heL2L = heL1.next().twin();
  // C#: if (heL2L.Face.IsBoundaryLoop()) return null;
  if (heL2L.face().isBoundaryLoop()) return {};
  // C#: var heL3L = heL2L.Prev.Twin;
  Halfedge heL3L = heL2L.prevOrbitFace().twin();
  // C#: if (heL3L.Face.IsBoundaryLoop()) return null;
  if (heL3L.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2L.Face, heL3L.Face];
  return {f0, fL1, heL2L.face(), heL3L.face()};
}

// C#: private static Hedge.Face[]? GetFacesL3R(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL3R(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2R = heL1.Prev.Twin;
  Halfedge heL2R = heL1.prevOrbitFace().twin();
  // C#: if (heL2R.Face.IsBoundaryLoop()) return null;
  if (heL2R.face().isBoundaryLoop()) return {};
  // C#: var heL3R = heL2R.Next.Twin;
  Halfedge heL3R = heL2R.next().twin();
  // C#: if (heL3R.Face.IsBoundaryLoop()) return null;
  if (heL3R.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2R.Face, heL3R.Face];
  return {f0, fL1, heL2R.face(), heL3R.face()};
}

// C#: private static Hedge.Face[]? GetFacesL4L(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL4L(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2L = heL1.Next.Twin;
  Halfedge heL2L = heL1.next().twin();
  // C#: if (heL2L.Face.IsBoundaryLoop()) return null;
  if (heL2L.face().isBoundaryLoop()) return {};
  // C#: var heL3L = heL2L.Prev.Twin;
  Halfedge heL3L = heL2L.prevOrbitFace().twin();
  // C#: if (heL3L.Face.IsBoundaryLoop()) return null;
  if (heL3L.face().isBoundaryLoop()) return {};
  // C#: var heL4L = heL3L.Next.Twin;
  Halfedge heL4L = heL3L.next().twin();
  // C#: if (heL4L.Face.IsBoundaryLoop()) return null;
  if (heL4L.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2L.Face, heL3L.Face, heL4L.Face];
  return {f0, fL1, heL2L.face(), heL3L.face(), heL4L.face()};
}

// C#: private static Hedge.Face[]? GetFacesL4R(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL4R(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2R = heL1.Prev.Twin;
  Halfedge heL2R = heL1.prevOrbitFace().twin();
  // C#: if (heL2R.Face.IsBoundaryLoop()) return null;
  if (heL2R.face().isBoundaryLoop()) return {};
  // C#: var heL3R = heL2R.Next.Twin;
  Halfedge heL3R = heL2R.next().twin();
  // C#: if (heL3R.Face.IsBoundaryLoop()) return null;
  if (heL3R.face().isBoundaryLoop()) return {};
  // C#: var heL4R = heL3R.Prev.Twin;
  Halfedge heL4R = heL3R.prevOrbitFace().twin();
  // C#: if (heL4R.Face.IsBoundaryLoop()) return null;
  if (heL4R.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2R.Face, heL3R.Face, heL4R.Face];
  return {f0, fL1, heL2R.face(), heL3R.face(), heL4R.face()};
}

// C#: private static Hedge.Face[]? GetFacesL5L(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL5L(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2L = heL1.Next.Twin;
  Halfedge heL2L = heL1.next().twin();
  // C#: if (heL2L.Face.IsBoundaryLoop()) return null;
  if (heL2L.face().isBoundaryLoop()) return {};
  // C#: var heL3L = heL2L.Prev.Twin;
  Halfedge heL3L = heL2L.prevOrbitFace().twin();
  // C#: if (heL3L.Face.IsBoundaryLoop()) return null;
  if (heL3L.face().isBoundaryLoop()) return {};
  // C#: var heL4L = heL3L.Next.Twin;
  Halfedge heL4L = heL3L.next().twin();
  // C#: if (heL4L.Face.IsBoundaryLoop()) return null;
  if (heL4L.face().isBoundaryLoop()) return {};
  // C#: var heL5L = heL4L.Prev.Twin;
  Halfedge heL5L = heL4L.prevOrbitFace().twin();
  // C#: if (heL5L.Face.IsBoundaryLoop()) return null;
  if (heL5L.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2L.Face, heL3L.Face, heL4L.Face, heL5L.Face];
  return {f0, fL1, heL2L.face(), heL3L.face(), heL4L.face(), heL5L.face()};
}

// C#: private static Hedge.Face[]? GetFacesL5R(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL5R(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2R = heL1.Prev.Twin;
  Halfedge heL2R = heL1.prevOrbitFace().twin();
  // C#: if (heL2R.Face.IsBoundaryLoop()) return null;
  if (heL2R.face().isBoundaryLoop()) return {};
  // C#: var heL3R = heL2R.Next.Twin;
  Halfedge heL3R = heL2R.next().twin();
  // C#: if (heL3R.Face.IsBoundaryLoop()) return null;
  if (heL3R.face().isBoundaryLoop()) return {};
  // C#: var heL4R = heL3R.Prev.Twin;
  Halfedge heL4R = heL3R.prevOrbitFace().twin();
  // C#: if (heL4R.Face.IsBoundaryLoop()) return null;
  if (heL4R.face().isBoundaryLoop()) return {};
  // C#: var heL5R = heL4R.Next.Twin;
  Halfedge heL5R = heL4R.next().twin();
  // C#: if (heL5R.Face.IsBoundaryLoop()) return null;
  if (heL5R.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2R.Face, heL3R.Face, heL4R.Face, heL5R.Face];
  return {f0, fL1, heL2R.face(), heL3R.face(), heL4R.face(), heL5R.face()};
}

// C#: private static Hedge.Face[]? GetFacesL5LM(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL5LM(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2L = heL1.Next.Twin;
  Halfedge heL2L = heL1.next().twin();
  // C#: if (heL2L.Face.IsBoundaryLoop()) return null;
  if (heL2L.face().isBoundaryLoop()) return {};
  // C#: var heL3L = heL2L.Prev.Twin;
  Halfedge heL3L = heL2L.prevOrbitFace().twin();
  // C#: if (heL3L.Face.IsBoundaryLoop()) return null;
  if (heL3L.face().isBoundaryLoop()) return {};
  // C#: var heL4LM = heL3L.Prev.Twin;
  Halfedge heL4LM = heL3L.prevOrbitFace().twin();
  // C#: if (heL4LM.Face.IsBoundaryLoop()) return null;
  if (heL4LM.face().isBoundaryLoop()) return {};
  // C#: var heL5LM = heL4LM.Next.Twin;
  Halfedge heL5LM = heL4LM.next().twin();
  // C#: if (heL5LM.Face.IsBoundaryLoop()) return null;
  if (heL5LM.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2L.Face, heL3L.Face, heL4LM.Face, heL5LM.Face];
  return {f0, fL1, heL2L.face(), heL3L.face(), heL4LM.face(), heL5LM.face()};
}

// C#: private static Hedge.Face[]? GetFacesL5RM(Hedge.Face f0, Hedge.Face fL1, Hedge.Halfedge heL1)
static std::vector<Face> getFacesL5RM(Face f0, Face fL1, Halfedge heL1) {
  // C#: var heL2R = heL1.Prev.Twin;
  Halfedge heL2R = heL1.prevOrbitFace().twin();
  // C#: if (heL2R.Face.IsBoundaryLoop()) return null;
  if (heL2R.face().isBoundaryLoop()) return {};
  // C#: var heL3R = heL2R.Next.Twin;
  Halfedge heL3R = heL2R.next().twin();
  // C#: if (heL3R.Face.IsBoundaryLoop()) return null;
  if (heL3R.face().isBoundaryLoop()) return {};
  // C#: var heL4RM = heL3R.Next.Twin;
  Halfedge heL4RM = heL3R.next().twin();
  // C#: if (heL4RM.Face.IsBoundaryLoop()) return null;
  if (heL4RM.face().isBoundaryLoop()) return {};
  // C#: var heL5RM = heL4RM.Prev.Twin;
  Halfedge heL5RM = heL4RM.prevOrbitFace().twin();
  // C#: if (heL5RM.Face.IsBoundaryLoop()) return null;
  if (heL5RM.face().isBoundaryLoop()) return {};
  // C#: return [f0, fL1, heL2R.Face, heL3R.Face, heL4RM.Face, heL5RM.Face];
  return {f0, fL1, heL2R.face(), heL3R.face(), heL4RM.face(), heL5RM.face()};
}

// C#: private static Hedge.Face[]? GetCrossedFaces(Hedge.Corner corner, CandidateName candidateName)
static std::vector<Face> getCrossedFaces(Corner corner, CandidateName candidateName) {
  // C#: var f0 = corner.Face();
  Face f0 = corner.face();
  // C#: var heL1 = corner.Halfedge.Twin;
  Halfedge heL1 = corner.halfedge().twin();
  // C#: if (heL1.Face.IsBoundaryLoop()) return null;
  if (heL1.face().isBoundaryLoop()) return {};
  // C#: var fL1 = heL1.Face;
  Face fL1 = heL1.face();

  switch (candidateName) {
    case CandidateName::L1:   return {f0, fL1};
    case CandidateName::L2L:  return getFacesL2L(f0, fL1, heL1);
    case CandidateName::L2R:  return getFacesL2R(f0, fL1, heL1);
    case CandidateName::L3L:  return getFacesL3L(f0, fL1, heL1);
    case CandidateName::L3R:  return getFacesL3R(f0, fL1, heL1);
    case CandidateName::L4L:  return getFacesL4L(f0, fL1, heL1);
    case CandidateName::L4R:  return getFacesL4R(f0, fL1, heL1);
    case CandidateName::L5L:  return getFacesL5L(f0, fL1, heL1);
    case CandidateName::L5R:  return getFacesL5R(f0, fL1, heL1);
    case CandidateName::L5LM: return getFacesL5LM(f0, fL1, heL1);
    case CandidateName::L5RM: return getFacesL5RM(f0, fL1, heL1);
    default: return {};
  }
}

// C#: private static IEnumerable<(Hedge.Vertex neighbor, float distance, PathStep step)> GetNeighbors(...)
static std::vector<std::tuple<Vertex, double, PathStep>> getNeighbors(
    Vertex current, Vertex goal,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom) {

  std::vector<std::tuple<Vertex, double, PathStep>> result;
  std::unordered_set<size_t> seenVertices;

  // Edge neighbors
  for (Vertex neighbor : current.adjacentVertices()) {
    // Skip corners (2+ boundary edges) unless it's the goal
    int boundaryCount = 0;
    for (Halfedge he : neighbor.outgoingHalfedges()) {
      if (!he.isInterior()) boundaryCount++;
    }
    if (neighbor != goal && boundaryCount >= 2) continue;

    seenVertices.insert(neighbor.getIndex());

    double distance = norm(geom.vertexPositions[current] - geom.vertexPositions[neighbor]);

    PathStep step;
    step.from = current;
    step.to = neighbor;
    step.isExplorerJump = false;
    step.candidateName = CandidateName::None;
    step.distance = distance;

    result.push_back({neighbor, distance, step});
  }

  // Explorer candidates from L5 unfolding
  for (Corner corner : current.adjacentCorners()) {
    ExplorationResult exploreResult = explore(corner, geom);

    for (const CandidateVertex& candidate : exploreResult.getReachableCandidates()) {
      if (!candidate.hasVertex()) continue;
      if (seenVertices.count(candidate.vertex.getIndex())) continue;

      // Skip corners (2+ boundary edges) unless it's the goal
      int boundaryCount = 0;
      for (Halfedge he : candidate.vertex.outgoingHalfedges()) {
        if (!he.isInterior()) boundaryCount++;
      }
      if (candidate.vertex != goal && boundaryCount >= 2) continue;

      seenVertices.insert(candidate.vertex.getIndex());

      std::vector<Face> crossedFaces = getCrossedFaces(corner, candidate.name);

      PathStep step;
      step.from = current;
      step.to = candidate.vertex;
      step.isExplorerJump = true;
      step.candidateName = candidate.name;
      step.distance = candidate.distance;
      step.crossedFaces = crossedFaces;

      result.push_back({candidate.vertex, candidate.distance, step});
    }
  }

  return result;
}

// C#: private static (List<Hedge.Vertex> path, List<PathStep> steps, bool isComplete) ReconstructPath(...)
static std::tuple<std::vector<Vertex>, std::vector<PathStep>, bool> reconstructPath(
    Vertex from, Vertex to,
    const std::unordered_map<size_t, std::pair<Vertex, PathStep>>& cameFrom) {

  std::vector<Vertex> path;
  std::vector<PathStep> steps;
  bool isComplete = true;

  Vertex current = to;
  while (current != from) {
    path.push_back(current);

    auto it = cameFrom.find(current.getIndex());
    if (it != cameFrom.end()) {
      steps.push_back(it->second.second);
      current = it->second.first;
    } else {
      isComplete = false;
      break;
    }
  }

  path.push_back(from);

  std::reverse(path.begin(), path.end());
  std::reverse(steps.begin(), steps.end());

  return {path, steps, isComplete};
}

// C#: public static PathResult FindPath(Hedge.Vertex from, Hedge.Vertex to)
PathResult findPath(Vertex from, Vertex to,
                    ManifoldSurfaceMesh& mesh,
                    VertexPositionGeometry& geom) {

  PathResult result;

  // C#: if (from == to) return new PathResult([from], [], IsComplete: true, IsFallback: false);
  if (from == to) {
    result.path.push_back(from);
    result.isComplete = true;
    result.isFallback = false;
    return result;
  }

  Vector3 goalPos = geom.vertexPositions[to];

  std::unordered_map<size_t, double> gScore;
  gScore[from.getIndex()] = 0;

  std::unordered_map<size_t, std::pair<Vertex, PathStep>> cameFrom;

  auto cmp = [](const std::pair<double, Vertex>& a, const std::pair<double, Vertex>& b) {
    return a.first > b.first;  // Min-heap
  };
  std::priority_queue<std::pair<double, Vertex>,
                      std::vector<std::pair<double, Vertex>>,
                      decltype(cmp)> openSet(cmp);

  double startF = HEURISTIC_WEIGHT * norm(geom.vertexPositions[from] - goalPos);
  openSet.push(std::make_pair(startF, from));

  // C#: var closed = _closedPool ??= new HashSet<Hedge.Vertex>();
  std::unordered_set<size_t> closed;

  // C#: while (openSet.TryDequeue(out var current, out var dequeuedF))
  while (!openSet.empty()) {
    std::pair<double, Vertex> topPair = openSet.top();
    double dequeuedF = topPair.first;
    Vertex current = topPair.second;
    openSet.pop();

    // C#: if (current == to)
    if (current == to) {
      std::tuple<std::vector<Vertex>, std::vector<PathStep>, bool> recon = reconstructPath(from, to, cameFrom);
      result.path = std::get<0>(recon);
      result.steps = std::get<1>(recon);
      result.isComplete = std::get<2>(recon);
      result.isFallback = false;
      return result;
    }

    // C#: var currentG = gScore[current];
    double currentG = gScore[current.getIndex()];

    // C#: var currentF = currentG + HeuristicWeight * current.Position.DistanceTo(goalPos);
    double currentF = currentG + HEURISTIC_WEIGHT * norm(geom.vertexPositions[current] - goalPos);

    // C#: if (dequeuedF > currentF + 1e-6f) continue;
    if (dequeuedF > currentF + 1e-6) continue;

    // C#: if (!closed.Add(current)) continue;
    if (closed.count(current.getIndex())) continue;
    closed.insert(current.getIndex());

    // C#: var neighbors = GetNeighbors(current, to);
    auto neighbors = getNeighbors(current, to, mesh, geom);

    // C#: foreach (var (neighbor, distance, step) in neighbors)
    for (size_t ni = 0; ni < neighbors.size(); ni++) {
      Vertex neighbor = std::get<0>(neighbors[ni]);
      double distance = std::get<1>(neighbors[ni]);
      PathStep step = std::get<2>(neighbors[ni]);

      // C#: if (closed.Contains(neighbor)) continue;
      if (closed.count(neighbor.getIndex())) continue;

      // C#: var tentativeG = currentG + distance;
      double tentativeG = currentG + distance;

      // C#: if (!gScore.TryGetValue(neighbor, out var existingG) || tentativeG < existingG)
      std::unordered_map<size_t, double>::iterator it = gScore.find(neighbor.getIndex());
      if (it == gScore.end() || tentativeG < it->second) {
        cameFrom[neighbor.getIndex()] = std::make_pair(current, step);
        gScore[neighbor.getIndex()] = tentativeG;

        // C#: var f = tentativeG + HeuristicWeight * neighbor.Position.DistanceTo(goalPos);
        double f = tentativeG + HEURISTIC_WEIGHT * norm(geom.vertexPositions[neighbor] - goalPos);
        openSet.push(std::make_pair(f, neighbor));
      }
    }
  }

  // C#: throw new InvalidOperationException("A* failed...");
  // Return empty result instead of throwing
  result.isComplete = false;
  result.isFallback = true;
  return result;
}

// C#: public static (List<Hedge.Vertex> path, List<VertexPathStep> steps) FindGeodesicPath(...)
std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> findGeodesicPath(
    Vertex from, Vertex to,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom) {

  PathResult result = findPath(from, to, mesh, geom);

  std::vector<VertexPathStep> vertexSteps;
  for (const PathStep& step : result.steps) {
    VertexPathStep vs;
    vs.from = step.from;
    vs.to = step.to;
    vs.isApexJump = step.isExplorerJump;
    vs.crossedFaces = step.crossedFaces;
    vertexSteps.push_back(vs);
  }

  return {result.path, vertexSteps};
}

// C#: public static float ComputePathDistance(List<PathStep> steps)
double computePathDistance(const std::vector<PathStep>& steps) {
  double total = 0;
  for (const PathStep& s : steps) {
    total += s.distance;
  }
  return total;
}

// ============================================================================
// VertexPathToFaceStripConverter - Line-by-line port
// C# from VertexPathToFaceStripConverter.cs and FaceStripWalker.cs
// ============================================================================

// C#: Helper to find index of face in vector (-1 if not found)
static int findFaceIndex(const std::vector<Face>& faces, Face f) {
  for (size_t i = 0; i < faces.size(); i++) {
    if (faces[i] == f) return static_cast<int>(i);
  }
  return -1;
}

// C#: Helper to count faces containing a vertex
static int countFacesWithVertex(const std::vector<Face>& faces, Vertex v) {
  int count = 0;
  for (Face f : faces) {
    if (faceContainsVertex(f, v)) count++;
  }
  return count;
}

// C#: private static Hedge.Face? HandleApexJump(...)
static Face handleApexJump(
    std::vector<Face>& faces,
    Face currentFace,
    const VertexPathStep& step,
    Vertex entryVertex,
    Vertex exitVertex) {

  // C#: if (step.CrossedFaces == null || step.CrossedFaces.Length == 0)
  //     return currentFace;
  if (step.crossedFaces.empty()) {
    return currentFace;
  }

  // C#: var firstFace = step.CrossedFaces[0];
  Face firstFace = step.crossedFaces[0];

  // C#: // If we have a current face and it's not adjacent to the first crossed face, walk to connect
  // C#: if (currentFace != null && currentFace != firstFace && !FaceStripUtils.AreFacesAdjacent(currentFace, firstFace))
  if (currentFace != Face() && currentFace != firstFace && !areFacesAdjacent(currentFace, firstFace)) {
    // C#: var sharedVertex = step.From;
    Vertex sharedVertex = step.from;

    // C#: if (FaceStripUtils.FaceContainsVertex(currentFace, sharedVertex) && FaceStripUtils.FaceContainsVertex(firstFace, sharedVertex))
    if (faceContainsVertex(currentFace, sharedVertex) && faceContainsVertex(firstFace, sharedVertex)) {
      // C#: var walkCW = WalkToFace(currentFace, firstFace, sharedVertex, WalkDirection.Clockwise);
      WalkResult walkCW = walkToFace(currentFace, firstFace, sharedVertex, WalkDirection::Clockwise);
      // C#: var walkCCW = WalkToFace(currentFace, firstFace, sharedVertex, WalkDirection.CounterClockwise);
      WalkResult walkCCW = walkToFace(currentFace, firstFace, sharedVertex, WalkDirection::CounterClockwise);

      // C#: int CountBadFaces(WalkResult w) => w.Faces.Count(f =>
      //         FaceStripUtils.FaceContainsVertex(f, entryVertex) || FaceStripUtils.FaceContainsVertex(f, exitVertex));
      auto countBadFaces = [&](const WalkResult& w) -> int {
        int count = 0;
        for (Face f : w.faces) {
          if (faceContainsVertex(f, entryVertex) || faceContainsVertex(f, exitVertex)) {
            count++;
          }
        }
        return count;
      };

      // C#: var cwBad = walkCW.ReachedTarget ? CountBadFaces(walkCW) : int.MaxValue;
      int cwBad = walkCW.reachedTarget ? countBadFaces(walkCW) : std::numeric_limits<int>::max();
      // C#: var ccwBad = walkCCW.ReachedTarget ? CountBadFaces(walkCCW) : int.MaxValue;
      int ccwBad = walkCCW.reachedTarget ? countBadFaces(walkCCW) : std::numeric_limits<int>::max();

      // C#: WalkResult? walk = null;
      WalkResult* walk = nullptr;

      // C#: if (walkCW.ReachedTarget && walkCCW.ReachedTarget)
      if (walkCW.reachedTarget && walkCCW.reachedTarget) {
        // C#: if (cwBad < ccwBad) walk = walkCW;
        // C#: else if (ccwBad < cwBad) walk = walkCCW;
        // C#: else walk = walkCW.Faces.Count <= walkCCW.Faces.Count ? walkCW : walkCCW;
        if (cwBad < ccwBad) {
          walk = &walkCW;
        } else if (ccwBad < cwBad) {
          walk = &walkCCW;
        } else {
          walk = (walkCW.faces.size() <= walkCCW.faces.size()) ? &walkCW : &walkCCW;
        }
      }
      // C#: else if (walkCW.ReachedTarget) walk = walkCW;
      else if (walkCW.reachedTarget) {
        walk = &walkCW;
      }
      // C#: else if (walkCCW.ReachedTarget) walk = walkCCW;
      else if (walkCCW.reachedTarget) {
        walk = &walkCCW;
      }

      // C#: if (walk != null)
      if (walk != nullptr) {
        // C#: foreach (var f in walk.Faces)
        for (Face f : walk->faces) {
          // C#: var existingIdx = faces.IndexOf(f);
          int existingIdx = findFaceIndex(faces, f);
          // C#: if (existingIdx >= 0)
          if (existingIdx >= 0) {
            // C#: if (existingIdx + 1 < faces.Count)
            //         faces.RemoveRange(existingIdx + 1, faces.Count - existingIdx - 1);
            if (existingIdx + 1 < static_cast<int>(faces.size())) {
              faces.erase(faces.begin() + existingIdx + 1, faces.end());
            }
          }
          // C#: else { faces.Add(f); }
          else {
            faces.push_back(f);
          }
        }
      }
    }
  }

  // C#: // Add all crossed faces from the step
  // C#: Hedge.Face? lastFace = null;
  Face lastFace;

  // C#: foreach (var crossedFace in step.CrossedFaces)
  for (Face crossedFace : step.crossedFaces) {
    // C#: if (crossedFace.IsBoundaryLoop() || crossedFace.IsDeleted) continue;
    if (crossedFace.isBoundaryLoop()) continue;

    // C#: var existingIdx = faces.IndexOf(crossedFace);
    int existingIdx = findFaceIndex(faces, crossedFace);

    // C#: if (existingIdx >= 0)
    if (existingIdx >= 0) {
      // C#: if (existingIdx + 1 < faces.Count)
      //         faces.RemoveRange(existingIdx + 1, faces.Count - existingIdx - 1);
      if (existingIdx + 1 < static_cast<int>(faces.size())) {
        faces.erase(faces.begin() + existingIdx + 1, faces.end());
      }
    }
    // C#: else { faces.Add(crossedFace); }
    else {
      faces.push_back(crossedFace);
    }
    // C#: lastFace = crossedFace;
    lastFace = crossedFace;
  }

  // C#: return lastFace ?? currentFace;
  return (lastFace != Face()) ? lastFace : currentFace;
}

// C#: private static Hedge.Face? HandleEdgeStep(...)
static Face handleEdgeStep(
    std::vector<Face>& faces,
    Face currentFace,
    const VertexPathStep& step,
    const std::vector<VertexPathStep>& steps,
    size_t stepIndex,
    Vertex entryVertex,
    Vertex exitVertex,
    VertexPositionGeometry& geom) {

  // C#: var from = step.From;
  Vertex from = step.from;
  // C#: var to = step.To;
  Vertex to = step.to;
  // C#: var isFirst = stepIndex == 0;
  bool isFirst = (stepIndex == 0);

  // C#: if (isFirst)
  if (isFirst) {
    // C#: Hedge.Vertex? nextVertex = null;
    // C#: Hedge.Face? targetFace = null;
    Vertex nextVertex;
    Face targetFace;

    // C#: if (stepIndex + 1 < steps.Count)
    if (stepIndex + 1 < steps.size()) {
      // C#: var nextStep = steps[stepIndex + 1];
      const VertexPathStep& nextStep = steps[stepIndex + 1];
      // C#: nextVertex = nextStep.To;
      nextVertex = nextStep.to;

      // C#: if (nextStep.IsApexJump && nextStep.CrossedFaces is { Length: > 0 })
      //         targetFace = nextStep.CrossedFaces[0];
      if (nextStep.isApexJump && !nextStep.crossedFaces.empty()) {
        targetFace = nextStep.crossedFaces[0];
      }
    }

    // C#: currentFace = SelectFirstFace(from, to, nextVertex, targetFace);
    currentFace = selectFirstFace(from, to, nextVertex, geom, targetFace);

    // C#: if (currentFace != null && !faces.Contains(currentFace))
    //         faces.Add(currentFace);
    if (currentFace != Face() && findFaceIndex(faces, currentFace) < 0) {
      faces.push_back(currentFace);
    }

    // C#: return currentFace;
    return currentFace;
  }

  // C#: if (currentFace == null)
  if (currentFace == Face()) {
    // C#: currentFace = FaceStripUtils.GetSharedFace(from, to);
    currentFace = getSharedFace(from, to);
    // C#: if (currentFace != null && !faces.Contains(currentFace))
    //         faces.Add(currentFace);
    if (currentFace != Face() && findFaceIndex(faces, currentFace) < 0) {
      faces.push_back(currentFace);
    }
    // C#: return currentFace;
    return currentFace;
  }

  // C#: if (FaceStripUtils.FaceContainsEdge(currentFace, from, to))
  //         return currentFace;
  if (faceContainsEdge(currentFace, from, to)) {
    return currentFace;
  }

  // C#: var prevVertex = steps[stepIndex - 1].From;
  Vertex prevVertex = steps[stepIndex - 1].from;
  // C#: var direction = DetermineWalkDirection(prevVertex, from, to);
  WalkDirection direction = determineWalkDirection(prevVertex, from, to, geom);

  // C#: var walk = WalkToOutgoingEdge(currentFace, from, to, direction);
  WalkResult walk = walkToOutgoingEdge(currentFace, from, to, direction);

  // C#: if (!walk.ReachedTarget)
  if (!walk.reachedTarget) {
    // C#: var oppositeDir = direction == WalkDirection.Clockwise
    //         ? WalkDirection.CounterClockwise : WalkDirection.Clockwise;
    WalkDirection oppositeDir = (direction == WalkDirection::Clockwise)
      ? WalkDirection::CounterClockwise
      : WalkDirection::Clockwise;
    // C#: walk = WalkToOutgoingEdge(currentFace, from, to, oppositeDir);
    walk = walkToOutgoingEdge(currentFace, from, to, oppositeDir);
  }

  // C#: if (walk.ReachedTarget)
  if (walk.reachedTarget) {
    // C#: // Check if final face is already in list (backtracking)
    // C#: var finalIdx = walk.FinalFace != null ? faces.IndexOf(walk.FinalFace) : -1;
    int finalIdx = (walk.finalFace != Face()) ? findFaceIndex(faces, walk.finalFace) : -1;

    // C#: if (finalIdx >= 0)
    if (finalIdx >= 0) {
      // C#: if (finalIdx + 1 < faces.Count)
      //         faces.RemoveRange(finalIdx + 1, faces.Count - finalIdx - 1);
      if (finalIdx + 1 < static_cast<int>(faces.size())) {
        faces.erase(faces.begin() + finalIdx + 1, faces.end());
      }
      // C#: return walk.FinalFace;
      return walk.finalFace;
    }

    // C#: // Add intermediate faces
    // C#: foreach (var f in walk.Faces)
    for (Face f : walk.faces) {
      // C#: var existingIdx = faces.IndexOf(f);
      int existingIdx = findFaceIndex(faces, f);

      // C#: if (existingIdx >= 0)
      if (existingIdx >= 0) {
        // C#: if (existingIdx + 1 < faces.Count)
        //         faces.RemoveRange(existingIdx + 1, faces.Count - existingIdx - 1);
        if (existingIdx + 1 < static_cast<int>(faces.size())) {
          faces.erase(faces.begin() + existingIdx + 1, faces.end());
        }
      }
      // C#: else { faces.Add(f); }
      else {
        faces.push_back(f);
      }
    }
    // C#: return walk.FinalFace;
    return walk.finalFace;
  }

  // C#: // Walk failed - fall back to any face on the edge
  // C#: var fallbackFace = FaceStripUtils.GetSharedFace(from, to);
  Face fallbackFace = getSharedFace(from, to);
  // C#: if (fallbackFace != null && !faces.Contains(fallbackFace))
  //         faces.Add(fallbackFace);
  if (fallbackFace != Face() && findFaceIndex(faces, fallbackFace) < 0) {
    faces.push_back(fallbackFace);
  }
  // C#: return fallbackFace;
  return fallbackFace;
}

// C#: public static List<Hedge.Face> ConvertPath(...)
std::vector<Face> convertPath(const std::vector<Vertex>& path,
                               const std::vector<VertexPathStep>& steps,
                               Vertex entryVertex, Vertex exitVertex,
                               VertexPositionGeometry& geom) {

  // C#: if (path.Count < 2 || steps.Count == 0)
  //         return new List<Hedge.Face>();
  if (path.size() < 2 || steps.empty()) {
    return {};
  }

  // C#: var faces = new List<Hedge.Face>();
  std::vector<Face> faces;
  // C#: Hedge.Face? currentFace = null;
  Face currentFace;

  // C#: for (int i = 0; i < steps.Count; i++)
  for (size_t i = 0; i < steps.size(); i++) {
    // C#: var step = steps[i];
    const VertexPathStep& step = steps[i];

    // C#: if (step.IsApexJump)
    //         currentFace = HandleApexJump(faces, currentFace, step, entryVertex, exitVertex);
    // C#: else
    //         currentFace = HandleEdgeStep(faces, currentFace, step, steps, i, entryVertex, exitVertex);
    if (step.isApexJump) {
      currentFace = handleApexJump(faces, currentFace, step, entryVertex, exitVertex);
    } else {
      currentFace = handleEdgeStep(faces, currentFace, step, steps, i, entryVertex, exitVertex, geom);
    }

    // C#: if (currentFace == null) break;
    if (currentFace == Face()) break;
  }

  // C#: return faces;
  return faces;
}

// C#: public static (List<Hedge.Face> faces, List<Hedge.Vertex> vertexPath) FindFaceStripWithPath(...)
std::pair<std::vector<Face>, std::vector<Vertex>> findFaceStripWithPath(
    Vertex from, Vertex to,
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom) {

  std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> pathResult = findGeodesicPath(from, to, mesh, geom);
  std::vector<Vertex> path = pathResult.first;
  std::vector<VertexPathStep> steps = pathResult.second;

  if (path.size() < 2) {
    return std::make_pair(std::vector<Face>(), path);
  }

  std::vector<Face> faces = convertPath(path, steps, from, to, geom);

  return std::make_pair(faces, path);
}

// ============================================================================
// GeodesicPathJoiner - Simplified (single segment case)
// Full implementation would handle multi-segment paths
// ============================================================================

std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> joinPaths(
    const std::vector<std::pair<std::vector<Vertex>, std::vector<VertexPathStep>>>& segments,
    Vertex entryVertex, Vertex exitVertex,
    bool removeLoops) {

  if (segments.empty()) {
    return {{}, {}};
  }

  if (segments.size() == 1) {
    return segments[0];
  }

  // For multiple segments, concatenate (simplified - full impl handles loops)
  std::vector<Vertex> resultPath;
  std::vector<VertexPathStep> resultSteps;

  for (size_t i = 0; i < segments.size(); i++) {
    const std::vector<Vertex>& segPath = segments[i].first;
    const std::vector<VertexPathStep>& segSteps = segments[i].second;

    size_t startIdx = (i > 0 && !resultPath.empty() &&
                       !segPath.empty() && resultPath.back() == segPath[0]) ? 1 : 0;

    for (size_t j = startIdx; j < segPath.size(); j++) {
      resultPath.push_back(segPath[j]);
    }

    for (size_t k = 0; k < segSteps.size(); k++) {
      resultSteps.push_back(segSteps[k]);
    }
  }

  return std::make_pair(resultPath, resultSteps);
}

// ============================================================================
// CachedVeryDiscreteGeodesicPathfinder Implementation
// C# equivalent: CachedGeodesicPathfinder.cs
//
// Key optimizations:
// 1. Per-corner exploration cache (avoids recomputing L1-L7 unfolds)
// 2. Per-vertex boundary counts (avoids iterating halfedges repeatedly)
// 3. Per-vertex adjacency arrays (avoids creating iterators)
// 4. Reusable A* containers (avoids allocation per query)
// ============================================================================

CachedVeryDiscreteGeodesicPathfinder::CachedVeryDiscreteGeodesicPathfinder(
    ManifoldSurfaceMesh& mesh_,
    VertexPositionGeometry& geom_)
    : mesh(mesh_), geom(geom_),
      boundaryEdgeCount(mesh_, 0),
      adjacentVertices(mesh_),
      adjacentCorners(mesh_) {
  // Pre-compute vertex data eagerly for best performance
  precomputeVertexData();
}

void CachedVeryDiscreteGeodesicPathfinder::clearCache() {
  explorationCache.clear();
  cacheHits = 0;
  cacheMisses = 0;
}

void CachedVeryDiscreteGeodesicPathfinder::precomputeVertexData() {
  if (vertexDataComputed) return;

  // Compute boundary edge counts for all vertices
  // C#: public int BoundaryEdgeCountCached => _onBoundaryCount ??= GetAdjacentHalfedges().Count(h => h.IsOnBoundary);
  for (Vertex v : mesh.vertices()) {
    int count = 0;
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!he.isInterior()) count++;
    }
    boundaryEdgeCount[v] = count;
  }

  // Compute adjacent vertices for all vertices
  // C#: public Vertex[] AdjacentVerticesCached => _adjacentVertices ??= GetAdjacentVertices().ToArray();
  for (Vertex v : mesh.vertices()) {
    std::vector<Vertex> adj;
    for (Vertex neighbor : v.adjacentVertices()) {
      adj.push_back(neighbor);
    }
    adjacentVertices[v] = std::move(adj);
  }

  // Compute adjacent corners for all vertices
  // C#: public Corner[] AdjacentCornersCached => _adjacentCorners ??= GetAdjacentCorners().ToArray();
  for (Vertex v : mesh.vertices()) {
    std::vector<Corner> corners;
    for (Corner c : v.adjacentCorners()) {
      corners.push_back(c);
    }
    adjacentCorners[v] = std::move(corners);
  }

  vertexDataComputed = true;
}

void CachedVeryDiscreteGeodesicPathfinder::ensureVertexData() {
  if (!vertexDataComputed) {
    precomputeVertexData();
  }
}

const ExplorationResult& CachedVeryDiscreteGeodesicPathfinder::getOrComputeExploration(Corner corner) {
  size_t idx = corner.getIndex();

  auto it = explorationCache.find(idx);
  if (it != explorationCache.end()) {
    cacheHits++;
    return it->second;
  }

  cacheMisses++;

  // Compute L5 exploration
  ExplorationResult result = explore(corner, geom);

  // Insert into cache and return reference to the cached value
  auto insertResult = explorationCache.emplace(idx, std::move(result));
  return insertResult.first->second;
}

std::vector<Face> CachedVeryDiscreteGeodesicPathfinder::getCrossedFaces(Corner corner, CandidateName candidateName) {
  return ::geometrycentral::surface::very_discrete_geodesic::getCrossedFaces(corner, candidateName);
}

const std::vector<std::tuple<Vertex, double, PathStep>>&
CachedVeryDiscreteGeodesicPathfinder::getNeighbors(Vertex current, Vertex goal, double currentG) {
  ensureVertexData();

  neighbors.clear();
  seenVertices.clear();

  // 1. Standard edge neighbors using cached adjacency
  const std::vector<Vertex>& adjVerts = adjacentVertices[current];
  for (Vertex neighbor : adjVerts) {
    if (neighbor != goal && boundaryEdgeCount[neighbor] >= 2)
      continue;

    seenVertices.insert(neighbor.getIndex());

    double distance = norm(geom.vertexPositions[current] - geom.vertexPositions[neighbor]);

    PathStep step;
    step.from = current;
    step.to = neighbor;
    step.isExplorerJump = false;
    step.candidateName = CandidateName::None;
    step.distance = distance;

    neighbors.push_back({neighbor, distance, step});
  }

  // 2. L5 explorer candidates using cached explorations
  const std::vector<Corner>& corners = adjacentCorners[current];
  for (Corner corner : corners) {
    const ExplorationResult& exploreResult = getOrComputeExploration(corner);

    // All L5 candidates (11 total)
    std::array<const CandidateVertex*, 11> allCandidates = {{
      &exploreResult.L1,
      &exploreResult.L2L, &exploreResult.L3L, &exploreResult.L4L, &exploreResult.L5L, &exploreResult.L5LM,
      &exploreResult.L2R, &exploreResult.L3R, &exploreResult.L4R, &exploreResult.L5R, &exploreResult.L5RM
    }};

    for (const CandidateVertex* candidate : allCandidates) {
      if (!candidate->hasVertex() || !candidate->isReachable)
        continue;

      if (seenVertices.count(candidate->vertex.getIndex()))
        continue;
      seenVertices.insert(candidate->vertex.getIndex());

      if (candidate->vertex != goal && boundaryEdgeCount[candidate->vertex] >= 2)
        continue;

      PathStep step;
      step.from = current;
      step.to = candidate->vertex;
      step.isExplorerJump = true;
      step.candidateName = candidate->name;
      step.distance = candidate->distance;
      step.crossedFaces = getCrossedFaces(corner, candidate->name);

      neighbors.push_back({candidate->vertex, candidate->distance, std::move(step)});
    }
  }

  return neighbors;
}

std::tuple<std::vector<Vertex>, std::vector<PathStep>, bool>
CachedVeryDiscreteGeodesicPathfinder::reconstructPath(Vertex from, Vertex to) {
  std::vector<Vertex> path;
  std::vector<PathStep> steps;
  bool isComplete = true;

  Vertex current = to;
  while (current != from) {
    path.push_back(current);

    auto it = cameFrom.find(current.getIndex());
    if (it != cameFrom.end()) {
      steps.push_back(it->second.second);
      current = it->second.first;
    } else {
      isComplete = false;
      break;
    }
  }

  path.push_back(from);

  std::reverse(path.begin(), path.end());
  std::reverse(steps.begin(), steps.end());

  return {path, steps, isComplete};
}

PathResult CachedVeryDiscreteGeodesicPathfinder::findPath(Vertex from, Vertex to) {
  PathResult result;

  if (from == to) {
    result.path.push_back(from);
    result.isComplete = true;
    result.isFallback = false;
    return result;
  }

  ensureVertexData();

  Vector3 goalPos = geom.vertexPositions[to];

  // Clear and reuse containers
  gScore.clear();
  cameFrom.clear();
  closed.clear();

  gScore[from.getIndex()] = 0;

  // Priority queue
  auto cmp = [](const std::pair<double, Vertex>& a, const std::pair<double, Vertex>& b) {
    return a.first > b.first;
  };
  std::priority_queue<std::pair<double, Vertex>,
                      std::vector<std::pair<double, Vertex>>,
                      decltype(cmp)> openSet(cmp);

  double startF = HEURISTIC_WEIGHT * norm(geom.vertexPositions[from] - goalPos);
  openSet.push(std::make_pair(startF, from));

  while (!openSet.empty()) {
    std::pair<double, Vertex> topPair = openSet.top();
    double dequeuedF = topPair.first;
    Vertex current = topPair.second;
    openSet.pop();

    if (current == to) {
      std::tuple<std::vector<Vertex>, std::vector<PathStep>, bool> recon = reconstructPath(from, to);
      result.path = std::get<0>(recon);
      result.steps = std::get<1>(recon);
      result.isComplete = std::get<2>(recon);
      result.isFallback = false;
      return result;
    }

    auto gIt = gScore.find(current.getIndex());
    if (gIt == gScore.end()) continue;
    double currentG = gIt->second;

    double currentF = currentG + HEURISTIC_WEIGHT * norm(geom.vertexPositions[current] - goalPos);
    if (dequeuedF > currentF + 1e-6)
      continue;

    if (closed.count(current.getIndex()))
      continue;
    closed.insert(current.getIndex());

    const auto& neighborList = getNeighbors(current, to, currentG);

    for (size_t i = 0; i < neighborList.size(); i++) {
      Vertex neighbor = std::get<0>(neighborList[i]);
      double distance = std::get<1>(neighborList[i]);
      const PathStep& step = std::get<2>(neighborList[i]);

      if (closed.count(neighbor.getIndex()))
        continue;

      double tentativeG = currentG + distance;

      auto it = gScore.find(neighbor.getIndex());
      if (it == gScore.end() || tentativeG < it->second) {
        cameFrom[neighbor.getIndex()] = std::make_pair(current, step);
        gScore[neighbor.getIndex()] = tentativeG;

        double f = tentativeG + HEURISTIC_WEIGHT * norm(geom.vertexPositions[neighbor] - goalPos);
        openSet.push(std::make_pair(f, neighbor));
      }
    }
  }

  result.isComplete = false;
  result.isFallback = true;
  return result;
}

std::pair<std::vector<Vertex>, std::vector<VertexPathStep>>
CachedVeryDiscreteGeodesicPathfinder::findGeodesicPath(Vertex from, Vertex to) {
  PathResult result = findPath(from, to);

  std::vector<VertexPathStep> vertexSteps;
  for (const PathStep& step : result.steps) {
    VertexPathStep vs;
    vs.from = step.from;
    vs.to = step.to;
    vs.isApexJump = step.isExplorerJump;
    vs.crossedFaces = step.crossedFaces;
    vertexSteps.push_back(vs);
  }

  return std::make_pair(result.path, vertexSteps);
}

std::pair<std::vector<Face>, std::vector<Vertex>>
CachedVeryDiscreteGeodesicPathfinder::findFaceStripWithPath(Vertex from, Vertex to) {
  std::pair<std::vector<Vertex>, std::vector<VertexPathStep>> geodesicResult = findGeodesicPath(from, to);
  std::vector<Vertex>& path = geodesicResult.first;
  std::vector<VertexPathStep>& steps = geodesicResult.second;

  if (path.size() < 2) {
    return std::make_pair(std::vector<Face>(), path);
  }

  std::vector<Face> faces = convertPath(path, steps, from, to, geom);
  return std::make_pair(faces, path);
}

} // namespace very_discrete_geodesic
} // namespace surface
} // namespace geometrycentral
