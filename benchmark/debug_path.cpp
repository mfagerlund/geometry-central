// Debug tool for investigating GFR vs FlipOut on specific vertex pairs
//
// Usage: debug_path <mesh.obj> <v0> <v1>

#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <iomanip>
#include <iostream>

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

  Vector3 pos0 = geometry->vertexPositions[v0];
  Vector3 pos1 = geometry->vertexPositions[v1];
  double euclidean = (pos1 - pos0).norm();

  std::cout << "============ PATH: V" << v0idx << " -> V" << v1idx << " ============" << std::endl;
  std::cout << "Euclidean distance: " << std::fixed << std::setprecision(6) << euclidean << std::endl;

  // Check if vertices share an edge
  bool sharesEdge = false;
  double edgeLength = 0;
  for (Halfedge he : v0.outgoingHalfedges()) {
    if (he.tipVertex() == v1) {
      sharesEdge = true;
      edgeLength = geometry->edgeLengths[he.edge()];
      break;
    }
  }
  if (sharesEdge) {
    std::cout << "Vertices share edge with length: " << edgeLength << std::endl;
  } else {
    std::cout << "Vertices do NOT share an edge" << std::endl;
  }
  std::cout << std::endl;

  // Run GFR
  std::cout << "=== GFR ===" << std::endl;
  auto gfrPath = computeFunnelGeodesic(*mesh, *geometry, v0, v1);
  std::cout << "Initial funnel length: " << gfrPath->initialFunnelLength() << std::endl;
  std::cout << "Final path length: " << gfrPath->length() << std::endl;
  std::cout << "Iterations: " << gfrPath->iterationCount() << std::endl;
  std::cout << "Faces in sleeve: " << gfrPath->faceCount() << std::endl;
  const auto& pathPts = gfrPath->getPath();
  std::cout << "Path points: " << pathPts.size() << std::endl;

  // Print path vertices
  std::cout << "Path: ";
  for (const auto& sp : pathPts) {
    if (sp.type == SurfacePointType::Vertex) {
      std::cout << "V" << sp.vertex.getIndex() << " ";
    } else if (sp.type == SurfacePointType::Edge) {
      std::cout << "E" << sp.edge.getIndex() << " ";
    } else {
      std::cout << "F" << sp.face.getIndex() << " ";
    }
  }
  std::cout << std::endl;

  // Analyze corners
  const auto& sleeveFaces = gfrPath->getSleeveFaces();
  auto flatPos = funnel_internal::flattenSleeve(sleeveFaces, v0, *geometry);
  auto portals = funnel_internal::buildPortals(sleeveFaces, flatPos);
  auto funnel = funnel_internal::runFunnel(portals, flatPos[v0], flatPos[v1]);
  auto corners = funnel_internal::analyzeCorners(sleeveFaces, funnel, flatPos);

  std::cout << "Corner analysis (" << corners.size() << " corners):" << std::endl;
  for (const auto& c : corners) {
    std::cout << "  V" << c.vertex.getIndex() << ": angleErrorDeg=" << c.angleErrorDeg
              << " wantsToFlip=" << (c.wantsToFlip() ? "YES" : "no") << std::endl;
  }
  std::cout << std::endl;

  // Run FlipOut
  std::cout << "=== FlipOut ===" << std::endl;
  auto flipNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*mesh, *geometry, v0, v1);

  // Get initial length before straightening
  double initialLength = flipNetwork->length();
  std::cout << "Initial (Dijkstra) length: " << initialLength << std::endl;

  flipNetwork->iterativeShorten();
  std::cout << "Final length: " << flipNetwork->length() << std::endl;
  std::cout << "Iterations: " << flipNetwork->nShortenIters << std::endl;
  std::cout << std::endl;

  // Compare
  std::cout << "=== COMPARISON ===" << std::endl;

  // Initial length comparison
  double initDiff = gfrPath->initialFunnelLength() - initialLength;
  double initPctDiff = 100.0 * initDiff / initialLength;
  std::cout << "Initial: GFR funnel=" << gfrPath->initialFunnelLength()
            << ", FlipOut Dijkstra=" << initialLength
            << " (" << std::showpos << initPctDiff << "%)" << std::noshowpos << std::endl;
  if (initDiff > 1e-9) {
    std::cout << "  WARNING: GFR initial funnel is LONGER than Dijkstra!" << std::endl;
  }

  // Final length comparison
  double diff = gfrPath->length() - flipNetwork->length();
  double pctDiff = 100.0 * diff / flipNetwork->length();
  std::cout << "Final: GFR=" << gfrPath->length()
            << ", FlipOut=" << flipNetwork->length()
            << " (" << std::showpos << pctDiff << "%)" << std::noshowpos << std::endl;

  if (gfrPath->length() < flipNetwork->length()) {
    std::cout << "Winner: GFR" << std::endl;
  } else if (flipNetwork->length() < gfrPath->length()) {
    std::cout << "Winner: FlipOut" << std::endl;
  } else {
    std::cout << "Winner: Tie" << std::endl;
  }

  return 0;
}
