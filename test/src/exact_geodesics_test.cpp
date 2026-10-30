#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ExactGeodesicSuite : public MeshAssetSuite {};

// The distance field does not depend on whether the caller holds the curvatures used by the saddle test
TEST_F(ExactGeodesicSuite, ResultIndependentOfRequiredQuantities) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("bob_small.ply", true)}) {
    a.printThyName();
    SurfaceMesh& mesh = *a.mesh;
    VertexPositionGeometry& geometry = *a.geometry;
    Vertex source = mesh.vertex(0);

    GeodesicAlgorithmExact mmp(mesh, geometry);
    mmp.propagate(source);
    VertexData<double> dist = mmp.getDistanceFunction();

    geometry.requireVertexGaussianCurvatures();
    GeodesicAlgorithmExact mmpHeld(mesh, geometry);
    mmpHeld.propagate(source);
    VertexData<double> distHeld = mmpHeld.getDistanceFunction();
    geometry.unrequireVertexGaussianCurvatures();

    for (Vertex v : mesh.vertices()) {
      EXPECT_EQ(dist[v], distHeld[v]);
    }
  }
}
