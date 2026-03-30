// VeryDiscreteGeodesic - A* pathfinder with L5 multi-face exploration
// Uses 2D unfolding to explore vertices up to 5 faces deep, finding shorter
// discrete path approximations than standard edge-only Dijkstra.

#include "geometrycentral/surface/very_discrete_geodesic.h"

#include <cmath>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace geometrycentral {
namespace surface {
namespace very_discrete_geodesic {

// ============================================================================
// Helper: V2.ComputeTriangleApex
// public static Vector2 ComputeTriangleApex(Vector2 a, Vector2 b, float distA, float distB, bool pickPositiveY)
// ============================================================================
Vector2 computeTriangleApex(Vector2 a, Vector2 b, double distA, double distB, bool pickPositiveY) {
  // Handle degenerate cases where apex coincides with a base point
  if (distA < 1e-10) return a;  // Apex at point a
  if (distB < 1e-10) return b;  // Apex at point b

  Vector2 ab = b - a;

  double abLen = norm(ab);

  if (abLen < 1e-10) return Vector2{std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN()};

  double cosA = (distA * distA + abLen * abLen - distB * distB) / (2.0 * distA * abLen);

  if (cosA > 1.0) cosA = 1.0;

  if (cosA < -1.0) cosA = -1.0;

  double sinA = std::sqrt(1.0 - cosA * cosA);

  Vector2 abNorm = ab / abLen;

  double perpX = -abNorm.y;
  double perpY = abNorm.x;

  double localX = distA * cosA;

  double localY = distA * sinA;

  if (!pickPositiveY) localY = -localY;

  double worldX = a.x + localX * abNorm.x + localY * perpX;

  double worldY = a.y + localX * abNorm.y + localY * perpY;

  return Vector2{worldX, worldY};
}

// ============================================================================
// Face strip utilities
// ============================================================================

Face getSharedFace(Vertex v1, Vertex v2) {
  for (Halfedge he : v1.outgoingHalfedges()) {
    // Check if v2 is any vertex in the same face as this halfedge
    if (he.tipVertex() == v2 || he.next().tipVertex() == v2) {
      Face face = he.face();

      if (face != Face() && !face.isBoundaryLoop()) {
        return face;
      }
    }
  }
  return Face();
}

std::vector<Face> getAllSharedFaces(Vertex v1, Vertex v2) {
  std::vector<Face> result;

  for (Halfedge he : v1.outgoingHalfedges()) {
    // Check if v2 is any vertex in the same face as this halfedge:
    // - he.tipVertex() == v2: v2 is directly connected to v1 via this edge
    // - he.next().tipVertex() == v2: v2 is the third vertex of the triangle
    // NOTE: he.next().next().tipVertex() is v1 itself, so we don't check that
    if (he.tipVertex() == v2 || he.next().tipVertex() == v2) {
      Face face = he.face();

      if (face != Face() && !face.isBoundaryLoop()) {
        bool found = false;
        for (Face f : result) {
          if (f == face) { found = true; break; }
        }
        if (!found) {
          result.push_back(face);
        }
      }
    }
  }
  return result;
}

bool areFacesAdjacent(Face f1, Face f2) {
  for (Halfedge he : f1.adjacentHalfedges()) {
    if (he.twin().face() == f2) {
      return true;
    }
  }
  return false;
}

bool faceContainsEdge(Face face, Vertex v1, Vertex v2) {
  for (Halfedge he : face.adjacentHalfedges()) {
    Vertex heV = he.tailVertex();
    Vertex heNextV = he.tipVertex();
    if ((heV == v1 && heNextV == v2) || (heV == v2 && heNextV == v1)) {
      return true;
    }
  }
  return false;
}

bool faceContainsVertex(Face face, Vertex vertex) {
  for (Vertex v : face.adjacentVertices()) {
    if (v == vertex) return true;
  }
  return false;
}

Halfedge findHalfedgeToVertex(Face face, Vertex vertex) {
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.tipVertex() == vertex) {
      return he;
    }
  }
  return Halfedge();
}

Halfedge findHalfedgeFromVertex(Face face, Vertex vertex) {
  for (Halfedge he : face.adjacentHalfedges()) {
    if (he.tailVertex() == vertex) {
      return he;
    }
  }
  return Halfedge();
}

// ============================================================================
// Explorer geometry helpers
// ============================================================================

bool segmentCrossesPortal(Vector2 v0, Vector2 target, Vector2 a, Vector2 b, double eps) {
  double dx = target.x - v0.x;

  double dy = target.y - v0.y;

  double crossA = dx * (a.y - v0.y) - dy * (a.x - v0.x);

  double crossB = dx * (b.y - v0.y) - dy * (b.x - v0.x);

  // Opposite signs = crosses (product negative)
  // Product near zero = grazes endpoint (also valid, handles numerical precision)
  return crossA * crossB <= eps;
}

Vector2 flattenVertex(Vertex v, Halfedge portal, Vector2 flatA, Vector2 flatB,
                      VertexPositionGeometry& geom) {
  double d1 = norm(geom.vertexPositions[v] - geom.vertexPositions[portal.tailVertex()]);

  double d2 = norm(geom.vertexPositions[v] - geom.vertexPositions[portal.twin().tailVertex()]);

  return computeTriangleApex(flatA, flatB, d1, d2, false);
}

// Version that uses a reference point to pick the correct side (matches flattenSleeve)
Vector2 flattenVertexWithRef(Vertex v, Halfedge portal, Vector2 flatA, Vector2 flatB,
                             Vector2 refPoint, VertexPositionGeometry& geom) {
  double d1 = norm(geom.vertexPositions[v] - geom.vertexPositions[portal.tailVertex()]);
  double d2 = norm(geom.vertexPositions[v] - geom.vertexPositions[portal.twin().tailVertex()]);

  // Use circleIntersect with reference point (same as flattenSleeve)
  double dx = flatB.x - flatA.x;
  double dy = flatB.y - flatA.y;
  double d = std::sqrt(dx * dx + dy * dy);

  if (d < 1e-14) return Vector2{flatA.x + d1, flatA.y};

  if (d > d1 + d2 + 1e-10 || d < std::abs(d1 - d2) - 1e-10)
    return Vector2{flatA.x + dx / d * d1, flatA.y + dy / d * d1};

  double a = (d1 * d1 - d2 * d2 + d * d) / (2.0 * d);
  double hSq = d1 * d1 - a * a;
  if (hSq < 0) hSq = 0;
  double h = std::sqrt(hSq);

  double invD = 1.0 / d;
  double dirX = dx * invD;
  double dirY = dy * invD;
  double midX = flatA.x + dirX * a;
  double midY = flatA.y + dirY * a;

  if (h < 1e-10) return Vector2{midX, midY};

  double perpX = -dirY * h;
  double perpY = dirX * h;

  double p1X = midX + perpX;
  double p1Y = midY + perpY;

  // Pick the side OPPOSITE to refPoint (matching flattenSleeve)
  double refSide = dx * (refPoint.y - flatA.y) - dy * (refPoint.x - flatA.x);
  double p1Side = dx * (p1Y - flatA.y) - dy * (p1X - flatA.x);

  if (p1Side * refSide < 0) {
    return Vector2{p1X, p1Y};
  } else {
    return Vector2{midX - perpX, midY - perpY};
  }
}

Vector2 getFlatPosition(Vertex v, Halfedge portal, Vector2 flatP1, Vector2 flatP2, Vector2 flatApex) {
  if (v == portal.tailVertex()) return flatP1;

  if (v == portal.twin().tailVertex()) return flatP2;

  return flatApex;
}


// Helper for creating unreachable candidates
static CandidateVertex makeUnreachable(CandidateName name) {
  CandidateVertex c;
  c.name = name;
  c.isReachable = false;
  return c;
}

// Block all candidates in an ExplorationResult
static void blockAllCandidates(ExplorationResult& result) {
  for (size_t i = 0; i < NUM_CANDIDATES; i++) {
    result.candidates[i] = makeUnreachable(static_cast<CandidateName>(i));
  }
}

// Portal pair for reachability checking
struct Portal { Vector2 a, b; };

// Halfedge transition direction: N = he.next().twin(), P = he.next().next().twin()
enum class Tr : uint8_t { N, P };

// Unified reachability check for 1-5 portals
static CandidateVertex checkReachability(CandidateName name, Vertex vertex, Vector2 flatV0,
                                         Vector2 flatTarget, Face targetFace,
                                         const Portal* portals, size_t numPortals, double scaleTolerance) {
  if (vertex == Vertex() || targetFace == Face() || std::isnan(flatTarget.x) || targetFace.isBoundaryLoop()) {
    return makeUnreachable(name);
  }

  double eps = scaleTolerance * 1e-6;
  for (size_t i = 0; i < numPortals; i++) {
    if (!segmentCrossesPortal(flatV0, flatTarget, portals[i].a, portals[i].b, eps)) {
      return makeUnreachable(name);
    }
  }

  CandidateVertex c;
  c.name = name;
  c.vertex = vertex;
  c.flatPosition = flatTarget;
  c.distance = norm(flatTarget - flatV0);
  c.isReachable = true;
  return c;
}

// State for one level of the trident unfolding tree.
// Accumulates portals so reachability can be checked directly.
struct UnfoldLevel {
  Halfedge he;
  Face face;
  Vertex apex;
  Vector2 flat_he_a, flat_he_b, flat_apex;
  Portal portals[5];
  uint8_t portalCount = 0;
  bool computed = false;
  bool isBoundary = true;
};

// Unfold one face in direction N (he.next().twin()) or P (he.next().next().twin()).
// If prev is boundary or not computed, returns an uncommitted level.
static UnfoldLevel unfoldNext(const UnfoldLevel& prev, Tr direction, VertexPositionGeometry& geom) {
  UnfoldLevel next;
  if (!prev.computed || prev.isBoundary) return next;

  next.he = (direction == Tr::N) ? prev.he.next().twin() : prev.he.next().next().twin();
  next.face = next.he.face();
  next.isBoundary = next.face.isBoundaryLoop();
  next.apex = next.he.next().next().vertex();

  next.flat_he_a = getFlatPosition(next.he.tailVertex(), prev.he, prev.flat_he_a, prev.flat_he_b, prev.flat_apex);
  next.flat_he_b = getFlatPosition(next.he.twin().tailVertex(), prev.he, prev.flat_he_a, prev.flat_he_b, prev.flat_apex);

  // The "fell off" vertex: N transition → tailVertex, P transition → tipVertex
  Vertex fellOff = (direction == Tr::N) ? prev.he.tailVertex() : prev.he.tipVertex();
  Vector2 ref = getFlatPosition(fellOff, prev.he, prev.flat_he_a, prev.flat_he_b, prev.flat_apex);
  next.flat_apex = flattenVertexWithRef(next.apex, next.he, next.flat_he_a, next.flat_he_b, ref, geom);

  // Accumulate parent portals + our own
  for (uint8_t i = 0; i < prev.portalCount; i++) next.portals[i] = prev.portals[i];
  next.portals[prev.portalCount] = {next.flat_he_a, next.flat_he_b};
  next.portalCount = prev.portalCount + 1;
  next.computed = true;
  return next;
}

// ============================================================================
// VeryDiscreteGeodesicExplorer - L5 exploration (fixed depth)
// ============================================================================

ExplorationResult explore(Corner corner, VertexPositionGeometry& geom) {
  ExplorationResult result;
  result.corner = corner;

  Vertex v0 = corner.vertex();
  Face fL1Raw = corner.halfedge().twin().face();

  if (fL1Raw.isBoundaryLoop()) {
    blockAllCandidates(result);
    return result;
  }

  Halfedge heL1 = corner.halfedge().twin();
  Vertex vL1 = heL1.next().next().vertex();

  Vector3 v1_3d = geom.vertexPositions[heL1.tailVertex()];
  Vector3 v2_3d = geom.vertexPositions[heL1.twin().tailVertex()];
  double portalLen = norm(v1_3d - v2_3d);

  Vector2 flat_a = {-portalLen / 2.0, 0.0};
  Vector2 flat_b = { portalLen / 2.0, 0.0};

  Vector2 flat_v0 = computeTriangleApex(flat_a, flat_b,
      norm(geom.vertexPositions[v0] - v1_3d),
      norm(geom.vertexPositions[v0] - v2_3d), true);

  if (std::isnan(flat_v0.x)) {
    blockAllCandidates(result);
    return result;
  }

  // Reference for L1: the third vertex of f0 that "falls off"
  Vertex v0_third = corner.halfedge().next().tipVertex();
  Vector2 flat_ref = computeTriangleApex(flat_a, flat_b,
      norm(geom.vertexPositions[v0_third] - v1_3d),
      norm(geom.vertexPositions[v0_third] - v2_3d), true);
  Vector2 flat_L1 = flattenVertexWithRef(vL1, heL1, flat_a, flat_b, flat_ref, geom);

  // Base level (L1)
  UnfoldLevel base;
  base.he = heL1;
  base.face = heL1.face();
  base.apex = vL1;
  base.flat_he_a = flat_a;
  base.flat_he_b = flat_b;
  base.flat_apex = flat_L1;
  base.portals[0] = {flat_a, flat_b};
  base.portalCount = 1;
  base.computed = true;
  base.isBoundary = false;

  // Unfold the trident: left (N,P,N,P), right (P,N,P,N), middle branches from L3
  UnfoldLevel L2L  = unfoldNext(base, Tr::N, geom);
  UnfoldLevel L3L  = unfoldNext(L2L,  Tr::P, geom);
  UnfoldLevel L4L  = unfoldNext(L3L,  Tr::N, geom);
  UnfoldLevel L5L  = unfoldNext(L4L,  Tr::P, geom);
  UnfoldLevel L4LM = unfoldNext(L3L,  Tr::P, geom);  // Middle branch from L3L
  UnfoldLevel L5LM = unfoldNext(L4LM, Tr::N, geom);

  UnfoldLevel L2R  = unfoldNext(base, Tr::P, geom);
  UnfoldLevel L3R  = unfoldNext(L2R,  Tr::N, geom);
  UnfoldLevel L4R  = unfoldNext(L3R,  Tr::P, geom);
  UnfoldLevel L5R  = unfoldNext(L4R,  Tr::N, geom);
  UnfoldLevel L4RM = unfoldNext(L3R,  Tr::N, geom);  // Middle branch from L3R
  UnfoldLevel L5RM = unfoldNext(L4RM, Tr::P, geom);

  // Check reachability using accumulated portals
  auto check = [&](CandidateName name, const UnfoldLevel& lvl) -> CandidateVertex {
    if (!lvl.computed) return makeUnreachable(name);
    return checkReachability(name, lvl.apex, flat_v0, lvl.flat_apex,
                             lvl.face, lvl.portals, lvl.portalCount, portalLen);
  };

  result[CandidateName::L1]   = check(CandidateName::L1,   base);
  result[CandidateName::L2L]  = check(CandidateName::L2L,  L2L);
  result[CandidateName::L3L]  = check(CandidateName::L3L,  L3L);
  result[CandidateName::L4L]  = check(CandidateName::L4L,  L4L);
  result[CandidateName::L5L]  = check(CandidateName::L5L,  L5L);
  result[CandidateName::L5LM] = check(CandidateName::L5LM, L5LM);
  result[CandidateName::L2R]  = check(CandidateName::L2R,  L2R);
  result[CandidateName::L3R]  = check(CandidateName::L3R,  L3R);
  result[CandidateName::L4R]  = check(CandidateName::L4R,  L4R);
  result[CandidateName::L5R]  = check(CandidateName::L5R,  L5R);
  result[CandidateName::L5RM] = check(CandidateName::L5RM, L5RM);

  return result;
}

// ============================================================================
// Face strip walker
// ============================================================================

WalkDirection determineWalkDirection(Vertex prev, Vertex current, Vertex next,
                                     VertexPositionGeometry& geom) {
  Vector3 reference = geom.vertexPositions[prev] - geom.vertexPositions[current];

  Vector3 outDir = geom.vertexPositions[next] - geom.vertexPositions[current];

  // Compute vertex normal (average of adjacent face normals)
  Vector3 normal = {0, 0, 0};
  for (Face f : current.adjacentFaces()) {
    normal += geom.faceNormal(f);
  }
  normal = unit(normal);

  // Signed angle computation
  Vector3 crossProd = cross(outDir, reference);
  double sinAngle = dot(crossProd, normal);

  bool isLeftTurn = sinAngle < 0;

  return isLeftTurn ? WalkDirection::CounterClockwise : WalkDirection::Clockwise;
}

// Compute corner angle at vertex in face (in radians)
static double cornerAngleAtVertex(Face face, Vertex vertex, VertexPositionGeometry& geom) {
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
  return std::acos(cosAngle);
}

// Compute total angle sum for faces in a walk result
static double computeWalkAngleSum(const WalkResult& walk, Vertex vertex, VertexPositionGeometry& geom) {
  double angleSum = 0.0;
  for (Face f : walk.faces) {
    angleSum += cornerAngleAtVertex(f, vertex, geom);
  }
  return angleSum;
}

WalkResult walkToOutgoingEdge(Face startFace, Vertex vertex, Vertex targetVertex,
                               WalkDirection direction) {
  WalkResult result;

  Halfedge startHe = findHalfedgeToVertex(startFace, vertex);

  if (startHe == Halfedge()) {
    result.error = "No halfedge in face points to vertex";
    return result;
  }

  if (faceContainsEdge(startFace, vertex, targetVertex)) {
    result.finalFace = startFace;
    result.reachedTarget = true;
    return result;
  }

  Halfedge currentHe = (direction == WalkDirection::Clockwise)
    ? startHe.next().twin()
    : startHe.twin();

  const int maxSteps = 100;

  for (int i = 0; i < maxSteps; i++) {
    if (currentHe.isInterior() == false || currentHe.face().isBoundaryLoop()) {
      result.error = "Hit boundary while walking around vertex";
      return result;
    }

    Face currentFace = currentHe.face();

    if (currentFace == startFace) {
      result.error = "Wrapped back to start face without finding edge to target";
      return result;
    }

    result.faces.push_back(currentFace);

    if (faceContainsEdge(currentFace, vertex, targetVertex)) {
      result.finalFace = currentFace;
      result.reachedTarget = true;
      return result;
    }

    if (direction == WalkDirection::Clockwise) {
      Halfedge leaving = findHalfedgeFromVertex(currentFace, vertex);
      if (leaving == Halfedge()) {
        result.error = "No halfedge in face leaves from vertex";
        return result;
      }
      currentHe = leaving.twin();
    } else {
      Halfedge entering = findHalfedgeToVertex(currentFace, vertex);
      if (entering == Halfedge()) {
        result.error = "No halfedge in face points to vertex";
        return result;
      }
      currentHe = entering.twin();
    }
  }

  result.error = "Exceeded max steps while walking around vertex";
  return result;
}

WalkResult walkToFace(Face startFace, Face targetFace, Vertex vertex, WalkDirection direction) {
  WalkResult result;

  if (startFace == targetFace) {
    result.finalFace = targetFace;
    result.reachedTarget = true;
    return result;
  }

  Halfedge startHe = findHalfedgeToVertex(startFace, vertex);
  if (startHe == Halfedge()) {
    result.error = "No starting halfedge";
    return result;
  }

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

Face selectFirstFace(Vertex v0, Vertex v1, Vertex v2,
                     VertexPositionGeometry& geom,
                     Face targetFace) {
  std::vector<Face> facesOnFirstEdge = getAllSharedFaces(v0, v1);

  if (facesOnFirstEdge.empty()) return Face();

  if (facesOnFirstEdge.size() == 1 || v2 == Vertex()) return facesOnFirstEdge[0];

  WalkDirection turnDirection = determineWalkDirection(v0, v1, v2, geom);

  WalkDirection oppositeDir = (turnDirection == WalkDirection::Clockwise)
    ? WalkDirection::CounterClockwise
    : WalkDirection::Clockwise;

  Face bestFace;
  int bestScore = std::numeric_limits<int>::max();

  for (Face face : facesOnFirstEdge) {
    for (WalkDirection dir : {turnDirection, oppositeDir}) {
      WalkResult testWalk;

      if (targetFace != Face()) {
        testWalk = walkToFace(face, targetFace, v1, dir);
      } else {
        testWalk = walkToOutgoingEdge(face, v1, v2, dir);
      }

      if (testWalk.reachedTarget) {
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
// Standalone A* pathfinder (for testing)
// ============================================================================

// Table-driven face traversal for getCrossedFaces
struct TraversalSpec { uint8_t depth; Tr steps[4]; };

static constexpr TraversalSpec TRAVERSALS[NUM_CANDIDATES] = {
  {0, {}},                                  // L1
  {1, {Tr::N}},                             // L2L
  {2, {Tr::N, Tr::P}},                      // L3L
  {3, {Tr::N, Tr::P, Tr::N}},               // L4L
  {4, {Tr::N, Tr::P, Tr::N, Tr::P}},        // L5L
  {4, {Tr::N, Tr::P, Tr::P, Tr::N}},        // L5LM
  {1, {Tr::P}},                             // L2R
  {2, {Tr::P, Tr::N}},                      // L3R
  {3, {Tr::P, Tr::N, Tr::P}},               // L4R
  {4, {Tr::P, Tr::N, Tr::P, Tr::N}},        // L5R
  {4, {Tr::P, Tr::N, Tr::N, Tr::P}},        // L5RM
};

std::vector<Face> getCrossedFaces(Corner corner, CandidateName candidateName) {
  if (candidateName == CandidateName::None) return {};

  Face f0 = corner.face();
  Halfedge he = corner.halfedge().twin();
  if (he.face().isBoundaryLoop()) return {};

  const TraversalSpec& spec = TRAVERSALS[static_cast<size_t>(candidateName)];
  std::vector<Face> faces = {f0, he.face()};
  faces.reserve(2 + spec.depth);

  for (uint8_t i = 0; i < spec.depth; i++) {
    he = (spec.steps[i] == Tr::N) ? he.next().twin() : he.next().next().twin();
    if (he.face().isBoundaryLoop()) return {};
    faces.push_back(he.face());
  }

  return faces;
}

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

    for (const CandidateVertex& candidate : exploreResult.candidates) {
      if (!candidate.hasVertex() || !candidate.isReachable) continue;
      if (seenVertices.count(candidate.vertex.getIndex())) continue;

      // Skip corners (2+ boundary edges) unless it's the goal
      int boundaryCount = 0;
      for (Halfedge he : candidate.vertex.outgoingHalfedges()) {
        if (!he.isInterior()) boundaryCount++;
      }
      if (candidate.vertex != goal && boundaryCount >= 2) continue;

      seenVertices.insert(candidate.vertex.getIndex());

      PathStep step;
      step.from = current;
      step.to = candidate.vertex;
      step.isExplorerJump = true;
      step.candidateName = candidate.name;
      step.distance = candidate.distance;
      step.sourceCorner = corner;

      result.push_back({candidate.vertex, candidate.distance, step});
    }
  }

  return result;
}

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

PathResult findPath(Vertex from, Vertex to,
                    ManifoldSurfaceMesh& mesh,
                    VertexPositionGeometry& geom) {

  PathResult result;

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

  std::unordered_set<size_t> closed;

  while (!openSet.empty()) {
    std::pair<double, Vertex> topPair = openSet.top();
    double dequeuedF = topPair.first;
    Vertex current = topPair.second;
    openSet.pop();

    if (current == to) {
      std::tuple<std::vector<Vertex>, std::vector<PathStep>, bool> recon = reconstructPath(from, to, cameFrom);
      result.path = std::get<0>(recon);
      result.steps = std::get<1>(recon);
      result.isComplete = std::get<2>(recon);
      result.isFallback = false;
      return result;
    }

    double currentG = gScore[current.getIndex()];

    double currentF = currentG + HEURISTIC_WEIGHT * norm(geom.vertexPositions[current] - goalPos);

    if (dequeuedF > currentF + 1e-6) continue;

    if (closed.count(current.getIndex())) continue;
    closed.insert(current.getIndex());

    auto neighbors = getNeighbors(current, to, mesh, geom);

    for (size_t ni = 0; ni < neighbors.size(); ni++) {
      Vertex neighbor = std::get<0>(neighbors[ni]);
      double distance = std::get<1>(neighbors[ni]);
      PathStep step = std::get<2>(neighbors[ni]);

      if (closed.count(neighbor.getIndex())) continue;

      double tentativeG = currentG + distance;

      std::unordered_map<size_t, double>::iterator it = gScore.find(neighbor.getIndex());
      if (it == gScore.end() || tentativeG < it->second) {
        cameFrom[neighbor.getIndex()] = std::make_pair(current, step);
        gScore[neighbor.getIndex()] = tentativeG;

        double f = tentativeG + HEURISTIC_WEIGHT * norm(geom.vertexPositions[neighbor] - goalPos);
        openSet.push(std::make_pair(f, neighbor));
      }
    }
  }

  // Return empty result instead of throwing
  result.isComplete = false;
  result.isFallback = true;
  return result;
}

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
    if (step.isExplorerJump) {
      vs.crossedFaces = getCrossedFaces(step.sourceCorner, step.candidateName);
    }
    vertexSteps.push_back(vs);
  }

  return {result.path, vertexSteps};
}

double computePathDistance(const std::vector<PathStep>& steps) {
  double total = 0;
  for (const PathStep& s : steps) {
    total += s.distance;
  }
  return total;
}

// ============================================================================
// Vertex path to face strip converter
// ============================================================================

static int findFaceIndex(const std::vector<Face>& faces, Face f) {
  for (size_t i = 0; i < faces.size(); i++) {
    if (faces[i] == f) return static_cast<int>(i);
  }
  return -1;
}

static Face handleApexJump(
    std::vector<Face>& faces,
    Face currentFace,
    const VertexPathStep& step,
    Vertex entryVertex,
    Vertex exitVertex) {

  if (step.crossedFaces.empty()) {
    return currentFace;
  }

  Face firstFace = step.crossedFaces[0];

  if (currentFace != Face() && currentFace != firstFace && !areFacesAdjacent(currentFace, firstFace)) {
    Vertex sharedVertex = step.from;

    if (faceContainsVertex(currentFace, sharedVertex) && faceContainsVertex(firstFace, sharedVertex)) {
      WalkResult walkCW = walkToFace(currentFace, firstFace, sharedVertex, WalkDirection::Clockwise);
      WalkResult walkCCW = walkToFace(currentFace, firstFace, sharedVertex, WalkDirection::CounterClockwise);

      auto countBadFaces = [&](const WalkResult& w) -> int {
        int count = 0;
        for (Face f : w.faces) {
          if (faceContainsVertex(f, entryVertex) || faceContainsVertex(f, exitVertex)) {
            count++;
          }
        }
        return count;
      };

      int cwBad = walkCW.reachedTarget ? countBadFaces(walkCW) : std::numeric_limits<int>::max();
      int ccwBad = walkCCW.reachedTarget ? countBadFaces(walkCCW) : std::numeric_limits<int>::max();

      WalkResult* walk = nullptr;

      if (walkCW.reachedTarget && walkCCW.reachedTarget) {
        if (cwBad < ccwBad) {
          walk = &walkCW;
        } else if (ccwBad < cwBad) {
          walk = &walkCCW;
        } else {
          walk = (walkCW.faces.size() <= walkCCW.faces.size()) ? &walkCW : &walkCCW;
        }
      }
      else if (walkCW.reachedTarget) {
        walk = &walkCW;
      }
      else if (walkCCW.reachedTarget) {
        walk = &walkCCW;
      }

      if (walk != nullptr) {
        for (Face f : walk->faces) {
          int existingIdx = findFaceIndex(faces, f);
          if (existingIdx >= 0) {
            if (existingIdx + 1 < static_cast<int>(faces.size())) {
              faces.erase(faces.begin() + existingIdx + 1, faces.end());
            }
          }
          else {
            faces.push_back(f);
          }
        }
      }
    }
  }

  Face lastFace;

  for (Face crossedFace : step.crossedFaces) {
    if (crossedFace.isBoundaryLoop()) continue;

    int existingIdx = findFaceIndex(faces, crossedFace);

    if (existingIdx >= 0) {
      if (existingIdx + 1 < static_cast<int>(faces.size())) {
        faces.erase(faces.begin() + existingIdx + 1, faces.end());
      }
    }
    else {
      faces.push_back(crossedFace);
    }
    lastFace = crossedFace;
  }

  return (lastFace != Face()) ? lastFace : currentFace;
}

static Face handleEdgeStep(
    std::vector<Face>& faces,
    Face currentFace,
    const VertexPathStep& step,
    const std::vector<VertexPathStep>& steps,
    size_t stepIndex,
    Vertex entryVertex,
    Vertex exitVertex,
    VertexPositionGeometry& geom) {

  Vertex from = step.from;
  Vertex to = step.to;
  bool isFirst = (stepIndex == 0);

  if (isFirst) {
    Vertex nextVertex;
    Face targetFace;

    if (stepIndex + 1 < steps.size()) {
      const VertexPathStep& nextStep = steps[stepIndex + 1];
      nextVertex = nextStep.to;

      if (nextStep.isApexJump && !nextStep.crossedFaces.empty()) {
        targetFace = nextStep.crossedFaces[0];
      }
    }

    currentFace = selectFirstFace(from, to, nextVertex, geom, targetFace);

    if (currentFace != Face() && findFaceIndex(faces, currentFace) < 0) {
      faces.push_back(currentFace);
    }

    return currentFace;
  }

  if (currentFace == Face()) {
    currentFace = getSharedFace(from, to);
    if (currentFace != Face() && findFaceIndex(faces, currentFace) < 0) {
      faces.push_back(currentFace);
    }
    return currentFace;
  }

  if (faceContainsEdge(currentFace, from, to)) {
    return currentFace;
  }

  // Try both directions and pick the one with smaller angle sum
  WalkResult walkCW = walkToOutgoingEdge(currentFace, from, to, WalkDirection::Clockwise);
  WalkResult walkCCW = walkToOutgoingEdge(currentFace, from, to, WalkDirection::CounterClockwise);

  WalkResult* walk = nullptr;

  if (walkCW.reachedTarget && walkCCW.reachedTarget) {
    // Both reached - pick the one with smaller angle sum
    double angleCW = computeWalkAngleSum(walkCW, from, geom);
    double angleCCW = computeWalkAngleSum(walkCCW, from, geom);
    walk = (angleCW <= angleCCW) ? &walkCW : &walkCCW;
  } else if (walkCW.reachedTarget) {
    walk = &walkCW;
  } else if (walkCCW.reachedTarget) {
    walk = &walkCCW;
  }

  if (walk != nullptr && walk->reachedTarget) {
    int finalIdx = (walk->finalFace != Face()) ? findFaceIndex(faces, walk->finalFace) : -1;

    if (finalIdx >= 0) {
      if (finalIdx + 1 < static_cast<int>(faces.size())) {
        faces.erase(faces.begin() + finalIdx + 1, faces.end());
      }
      return walk->finalFace;
    }

    for (Face f : walk->faces) {
      int existingIdx = findFaceIndex(faces, f);

      if (existingIdx >= 0) {
        if (existingIdx + 1 < static_cast<int>(faces.size())) {
          faces.erase(faces.begin() + existingIdx + 1, faces.end());
        }
      }
      else {
        faces.push_back(f);
      }
    }
    return walk->finalFace;
  }

  Face fallbackFace = getSharedFace(from, to);
  if (fallbackFace != Face() && findFaceIndex(faces, fallbackFace) < 0) {
    faces.push_back(fallbackFace);
  }
  return fallbackFace;
}

std::vector<Face> convertPath(const std::vector<Vertex>& path,
                               const std::vector<VertexPathStep>& steps,
                               Vertex entryVertex, Vertex exitVertex,
                               VertexPositionGeometry& geom) {

  if (path.size() < 2 || steps.empty()) {
    return {};
  }

  std::vector<Face> faces;
  Face currentFace;

  for (size_t i = 0; i < steps.size(); i++) {
    const VertexPathStep& step = steps[i];

    if (step.isApexJump) {
      currentFace = handleApexJump(faces, currentFace, step, entryVertex, exitVertex);
    } else {
      currentFace = handleEdgeStep(faces, currentFace, step, steps, i, entryVertex, exitVertex, geom);
    }

    if (currentFace == Face()) break;
  }

  return faces;
}

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
// CachedVeryDiscreteGeodesicPathfinder
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
  for (Vertex v : mesh.vertices()) {
    int count = 0;
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!he.isInterior()) count++;
    }
    boundaryEdgeCount[v] = count;
  }

  // Compute adjacent vertices for all vertices
  for (Vertex v : mesh.vertices()) {
    std::vector<Vertex> adj;
    for (Vertex neighbor : v.adjacentVertices()) {
      adj.push_back(neighbor);
    }
    adjacentVertices[v] = std::move(adj);
  }

  // Compute adjacent corners for all vertices
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

    for (const CandidateVertex& candidate : exploreResult.candidates) {
      if (!candidate.hasVertex() || !candidate.isReachable)
        continue;

      if (seenVertices.count(candidate.vertex.getIndex()))
        continue;
      seenVertices.insert(candidate.vertex.getIndex());

      if (candidate.vertex != goal && boundaryEdgeCount[candidate.vertex] >= 2)
        continue;

      PathStep step;
      step.from = current;
      step.to = candidate.vertex;
      step.isExplorerJump = true;
      step.candidateName = candidate.name;
      step.distance = candidate.distance;
      step.sourceCorner = corner;

      neighbors.push_back({candidate.vertex, candidate.distance, step});
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
    if (step.isExplorerJump) {
      vs.crossedFaces = getCrossedFaces(step.sourceCorner, step.candidateName);
    }
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
