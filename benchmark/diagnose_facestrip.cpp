// Diagnose face strip difference between VeryDiscreteGeodesic and Dijkstra
//
// Usage: diagnose_facestrip <mesh.obj> <v0> <v1>

#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/very_discrete_geodesic.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <iomanip>
#include <iostream>
#include <set>

using namespace geometrycentral;
using namespace geometrycentral::surface;

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <mesh.obj> <v0> <v1>" << std::endl;
    return 1;
  }

  std::string meshPath = argv[1];
  size_t v0idx = std::stoul(argv[2]);
  size_t v1idx = std::stoul(argv[3]);

  // Load mesh
  std::cout << "Loading mesh: " << meshPath << std::endl;
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshPath);

  std::cout << "Mesh: " << mesh->nVertices() << " vertices, "
            << mesh->nFaces() << " faces" << std::endl;
  std::cout << std::endl;

  Vertex v0 = mesh->vertex(v0idx);
  Vertex v1 = mesh->vertex(v1idx);

  std::cout << "============ PATH: V" << v0idx << " -> V" << v1idx << " ============" << std::endl;

  // Check adjacency of V1156 (v0) to key vertices
  std::cout << std::endl << "=== Edge Neighbors of V" << v0idx << " ===" << std::endl;
  std::cout << "Adjacent vertices: ";
  for (Vertex neighbor : v0.adjacentVertices()) {
    std::cout << "V" << neighbor.getIndex() << " ";
  }
  std::cout << std::endl;

  // Check if V1174 is adjacent (it's on Dijkstra path)
  bool v1174Adjacent = false;
  for (Vertex neighbor : v0.adjacentVertices()) {
    if (neighbor.getIndex() == 1174) {
      v1174Adjacent = true;
      double dist = norm(geometry->vertexPositions[v0] - geometry->vertexPositions[neighbor]);
      std::cout << "V1174 IS adjacent to V" << v0idx << ", edge distance = " << dist << std::endl;
      break;
    }
  }
  if (!v1174Adjacent) {
    std::cout << "V1174 is NOT adjacent to V" << v0idx << std::endl;
  }

  // Method 1: VeryDiscreteGeodesic (current GFR default)
  std::cout << std::endl << "=== VeryDiscreteGeodesic Face Strip ===" << std::endl;

  // Get the FULL path info including A* costs
  very_discrete_geodesic::CachedVeryDiscreteGeodesicPathfinder pathfinder(*mesh, *geometry);
  auto pathResult = pathfinder.findPath(v0, v1);

  // Also get face strip
  auto vdgResult = pathfinder.findFaceStripWithPath(v0, v1);
  auto& facesVDG = vdgResult.first;
  auto& vdgVertexPath = vdgResult.second;

  // Compute A* total cost
  double astarCost = very_discrete_geodesic::computePathDistance(pathResult.steps);
  std::cout << "A* path cost (what A* used): " << std::fixed << std::setprecision(6) << astarCost << std::endl;

  std::cout << "Vertex path: ";
  for (const auto& v : vdgVertexPath) std::cout << "V" << v.getIndex() << " ";
  std::cout << std::endl;

  // Show each step with its cost and type
  std::cout << "Path steps:" << std::endl;
  for (const auto& step : pathResult.steps) {
    std::cout << "  V" << step.from.getIndex() << " -> V" << step.to.getIndex()
              << ": distance=" << step.distance
              << (step.isExplorerJump ? " (L5 jump)" : " (edge)")
              << ", crossedFaces=[";
    for (size_t i = 0; i < step.crossedFaces.size(); i++) {
      if (i > 0) std::cout << ",";
      std::cout << "F" << step.crossedFaces[i].getIndex();
    }
    std::cout << "]" << std::endl;
  }

  std::cout << "Final sleeve faces: " << facesVDG.size() << std::endl;
  std::cout << "Face IDs: ";
  for (const auto& f : facesVDG) std::cout << "F" << f.getIndex() << " ";
  std::cout << std::endl;

  // Flatten and compute funnel
  auto flatPosVDG = funnel_internal::flattenSleeve(facesVDG, v0, *geometry);
  auto portalsVDG = funnel_internal::buildPortals(facesVDG, flatPosVDG);
  auto funnelVDG = funnel_internal::runFunnel(portalsVDG, flatPosVDG[v0], flatPosVDG[v1]);
  std::cout << "Funnel distance: " << std::fixed << std::setprecision(6) << funnelVDG.distance << std::endl;
  std::cout << "Waypoints: " << funnelVDG.waypointVertices.size() << std::endl;
  for (const auto& v : funnelVDG.waypointVertices) {
    std::cout << "  V" << v.getIndex() << std::endl;
  }

  // Show flat positions of key vertices
  std::cout << "Flat positions in sleeve:" << std::endl;
  std::cout << "  V" << v0idx << " (start): " << flatPosVDG[v0] << std::endl;
  for (const auto& v : vdgVertexPath) {
    if (v != v0) std::cout << "  V" << v.getIndex() << ": " << flatPosVDG[v] << std::endl;
  }
  // Show waypoint positions
  for (const auto& v : funnelVDG.waypointVertices) {
    std::cout << "  V" << v.getIndex() << " (waypoint): " << flatPosVDG[v] << std::endl;
  }

  // Check if direct line from v0 to intermediate vertex is blocked
  if (vdgVertexPath.size() >= 2) {
    Vertex intermediate = vdgVertexPath[1];
    Vector2 start2D = flatPosVDG[v0];
    Vector2 intermediate2D = flatPosVDG[intermediate];
    double directDist = (intermediate2D - start2D).norm();
    std::cout << "Direct 2D distance V" << v0idx << " -> V" << intermediate.getIndex()
              << ": " << directDist << std::endl;

    // Find the first waypoint that blocks this line
    for (const auto& wp : funnelVDG.waypointVertices) {
      Vector2 wp2D = flatPosVDG[wp];
      // Check if waypoint is between start and intermediate
      std::cout << "  Waypoint V" << wp.getIndex() << " at " << wp2D
                << " (dist from start: " << (wp2D - start2D).norm() << ")" << std::endl;
    }
  }

  // Verify isometry: check that flat edge lengths match 3D edge lengths
  std::cout << "Isometry check (flat vs 3D edge lengths):" << std::endl;
  for (const auto& f : facesVDG) {
    Halfedge he = f.halfedge();
    for (int j = 0; j < 3; j++) {
      Vertex va = he.tailVertex();
      Vertex vb = he.tipVertex();
      double len3D = (geometry->vertexPositions[va] - geometry->vertexPositions[vb]).norm();
      double len2D = (flatPosVDG[va] - flatPosVDG[vb]).norm();
      if (std::abs(len3D - len2D) > 1e-6) {
        std::cout << "  MISMATCH F" << f.getIndex() << ": V" << va.getIndex() << "-V" << vb.getIndex()
                  << " 3D=" << len3D << " 2D=" << len2D << " diff=" << (len2D - len3D) << std::endl;
      }
      he = he.next();
    }
  }

  // Show all portals
  std::cout << "Portals:" << std::endl;
  for (size_t i = 0; i < portalsVDG.size(); i++) {
    const auto& p = portalsVDG[i];
    std::cout << "  Portal " << i << ": L=V" << p.leftVert.getIndex() << " " << p.left
              << ", R=V" << p.rightVert.getIndex() << " " << p.right << std::endl;
  }

  // Check if direct line V1156->V705 crosses all portals correctly
  if (vdgVertexPath.size() >= 2) {
    Vertex intermediate = vdgVertexPath[1];
    Vector2 start2D = flatPosVDG[v0];
    Vector2 target2D = flatPosVDG[intermediate];
    std::cout << "Checking if line V" << v0idx << "->V" << intermediate.getIndex()
              << " crosses all portals:" << std::endl;

    for (size_t i = 0; i < portalsVDG.size(); i++) {
      const auto& p = portalsVDG[i];
      // Check if the target vertex is reached by this portal
      if (p.leftVert == intermediate || p.rightVert == intermediate) {
        std::cout << "  Portal " << i << ": Target vertex reached" << std::endl;
        break;
      }

      // Cross product test: does line cross portal?
      Vector2 a = p.left, b = p.right;
      double dx = target2D.x - start2D.x;
      double dy = target2D.y - start2D.y;
      double crossA = dx * (a.y - start2D.y) - dy * (a.x - start2D.x);
      double crossB = dx * (b.y - start2D.y) - dy * (b.x - start2D.x);
      bool crosses = (crossA * crossB <= 0);
      std::cout << "  Portal " << i << ": crossA=" << crossA << ", crossB=" << crossB
                << " -> " << (crosses ? "CROSSES" : "BLOCKED") << std::endl;
    }
  }

  // Method 2: Dijkstra (same as FlipOut uses)
  std::cout << std::endl << "=== Dijkstra Face Strip ===" << std::endl;

  // Get the Dijkstra edge path
  std::vector<Halfedge> dijkstraEdgePath = shortestEdgePath(*geometry, v0, v1);
  std::cout << "Dijkstra vertex path: V" << v0idx;
  for (const auto& he : dijkstraEdgePath) {
    std::cout << " V" << he.tipVertex().getIndex();
  }
  std::cout << std::endl;

  auto facesDijkstra = funnel_internal::buildFaceStrip(*mesh, *geometry, v0, v1);
  std::cout << "Faces: " << facesDijkstra.size() << std::endl;
  std::cout << "Face IDs: ";
  for (const auto& f : facesDijkstra) std::cout << "F" << f.getIndex() << " ";
  std::cout << std::endl;

  // Flatten and compute funnel
  auto flatPosDij = funnel_internal::flattenSleeve(facesDijkstra, v0, *geometry);
  auto portalsDij = funnel_internal::buildPortals(facesDijkstra, flatPosDij);
  auto funnelDij = funnel_internal::runFunnel(portalsDij, flatPosDij[v0], flatPosDij[v1]);
  std::cout << "Funnel distance: " << std::fixed << std::setprecision(6) << funnelDij.distance << std::endl;
  std::cout << "Waypoints: " << funnelDij.waypointVertices.size() << std::endl;
  for (const auto& v : funnelDij.waypointVertices) {
    std::cout << "  V" << v.getIndex() << std::endl;
  }

  // FlipOut Dijkstra length for reference
  std::cout << std::endl << "=== FlipOut Dijkstra ===" << std::endl;
  auto flipNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, v0, v1);
  double dijkstraLength = flipNetwork->length();
  std::cout << "Dijkstra edge path length: " << dijkstraLength << std::endl;

  // Calculate vertex path lengths (Euclidean distances)
  double vdgEuclidean = 0;
  std::cout << "VDG path segment distances:" << std::endl;
  for (size_t i = 1; i < vdgVertexPath.size(); i++) {
    double segDist = norm(geometry->vertexPositions[vdgVertexPath[i]] -
                         geometry->vertexPositions[vdgVertexPath[i-1]]);
    std::cout << "  V" << vdgVertexPath[i-1].getIndex() << " -> V" << vdgVertexPath[i].getIndex()
              << ": " << segDist << std::endl;
    vdgEuclidean += segDist;
  }

  // Calculate funnel waypoint distances
  std::cout << "VDG funnel waypoint distances:" << std::endl;
  std::cout << "  (funnel finds actual shortest path through corridor)" << std::endl;

  double dijkstraVertexLen = 0;
  Vertex prev = v0;
  for (const auto& he : dijkstraEdgePath) {
    dijkstraVertexLen += norm(geometry->vertexPositions[he.tipVertex()] -
                              geometry->vertexPositions[prev]);
    prev = he.tipVertex();
  }

  // Summary
  std::cout << std::endl << "=== SUMMARY ===" << std::endl;
  std::cout << "VeryDiscreteGeodesic:" << std::endl;
  std::cout << "  A* path cost:            " << astarCost << std::endl;
  std::cout << "  Vertex path 3D Euclidean:" << vdgEuclidean << std::endl;
  std::cout << "  Funnel through corridor: " << funnelVDG.distance << std::endl;
  std::cout << "Dijkstra:" << std::endl;
  std::cout << "  Vertex path edge length: " << dijkstraLength << std::endl;
  std::cout << "  Funnel through corridor: " << funnelDij.distance << std::endl;
  std::cout << std::endl;

  // THE KEY INVARIANT CHECK
  if (astarCost <= dijkstraLength + 1e-9) {
    std::cout << "A* found a path <= Dijkstra (" << astarCost << " <= " << dijkstraLength << ")" << std::endl;
    if (funnelVDG.distance > astarCost + 1e-9) {
      std::cout << "BUG: Funnel (" << funnelVDG.distance << ") > A* cost (" << astarCost << ")!" << std::endl;
      std::cout << "The sleeve does NOT properly contain the A* path!" << std::endl;
    } else {
      std::cout << "OK: Funnel <= A* cost (sleeve contains A* path)" << std::endl;
    }
  } else {
    std::cout << "A* found LONGER path than Dijkstra (" << astarCost << " > " << dijkstraLength << ")" << std::endl;
    std::cout << "A* should have found the edge path!" << std::endl;
  }

  if (funnelVDG.distance > dijkstraLength + 1e-9) {
    std::cout << "WARNING: VeryDiscreteGeodesic funnel is LONGER than Dijkstra edge path!" << std::endl;
    std::cout << "  Difference: " << (funnelVDG.distance - dijkstraLength)
              << " (+" << 100.0*(funnelVDG.distance - dijkstraLength)/dijkstraLength << "%)" << std::endl;
  }

  if (funnelDij.distance > dijkstraLength + 1e-9) {
    std::cout << "WARNING: Dijkstra funnel is LONGER than Dijkstra edge path!" << std::endl;
    std::cout << "  Difference: " << (funnelDij.distance - dijkstraLength)
              << " (+" << 100.0*(funnelDij.distance - dijkstraLength)/dijkstraLength << "%)" << std::endl;
  }

  if (funnelDij.distance <= dijkstraLength + 1e-9 && funnelVDG.distance > dijkstraLength + 1e-9) {
    std::cout << std::endl << "DIAGNOSIS: VeryDiscreteGeodesic is finding a different corridor" << std::endl;
    std::cout << "that doesn't contain the optimal Dijkstra path." << std::endl;

    // Find faces that are in Dijkstra but not VDG
    std::set<size_t> vdgFaceSet, dijFaceSet;
    for (const auto& f : facesVDG) vdgFaceSet.insert(f.getIndex());
    for (const auto& f : facesDijkstra) dijFaceSet.insert(f.getIndex());

    std::cout << std::endl << "Face difference analysis:" << std::endl;
    std::cout << "  Faces in Dijkstra but NOT in VDG: ";
    for (size_t fIdx : dijFaceSet) {
      if (vdgFaceSet.find(fIdx) == vdgFaceSet.end()) {
        std::cout << "F" << fIdx << " ";
      }
    }
    std::cout << std::endl;

    std::cout << "  Faces in VDG but NOT in Dijkstra: ";
    for (size_t fIdx : vdgFaceSet) {
      if (dijFaceSet.find(fIdx) == dijFaceSet.end()) {
        std::cout << "F" << fIdx << " ";
      }
    }
    std::cout << std::endl;

    std::cout << "  Shared faces: ";
    for (size_t fIdx : vdgFaceSet) {
      if (dijFaceSet.find(fIdx) != dijFaceSet.end()) {
        std::cout << "F" << fIdx << " ";
      }
    }
    std::cout << std::endl;
  }

  return 0;
}
