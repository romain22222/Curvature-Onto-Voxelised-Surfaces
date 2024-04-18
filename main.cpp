#include <DGtal/base/Common.h>

#include <utility>
#include "DGtal/geometry/meshes/NormalCycleComputer.h"
#include "DGtal/shapes/SurfaceMeshHelper.h"
#include "DGtal/io/writers/SurfaceMeshWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/QuantifiedColorMap.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace DGtal;
typedef SurfaceMesh< Z3i::RealPoint, Z3i::RealVector > SM;
typedef NormalCycleComputer< Z3i::RealPoint, Z3i::RealVector > NC;
typedef SurfaceMeshHelper< Z3i::RealPoint, Z3i::RealVector > SMH;
typedef SurfaceMeshWriter< Z3i::RealPoint, Z3i::RealVector > SMW;
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;
typedef polyscope::SurfaceMesh PolyMesh;
#define R 0.5

GradientColorMap<double> makeColorMap(double minv, double maxv) {
    GradientColorMap<double> gcm(minv, maxv);
    gcm.addColor(Color(0,0,255));
    gcm.addColor(Color(0,255,255));
    gcm.addColor(Color(255,255,255));
    gcm.addColor(Color(255,255,255));
    gcm.addColor(Color(255,0,255));
    return gcm;
}

void testDGtal() {
    DGtal::trace.info() << "Helloworld from DGtal ";
    DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;
}

void paper2022Checks() {
    SM smesh = SMH::makeTorus(3.0, 1.0, Z3i::RealPoint(), 20, 20, 0, SMH::NormalsType::VERTEX_NORMALS);
    NC nc(smesh);
    auto mu0 = nc.computeMu0();
    auto mu1 = nc.computeMu1();
    auto mu2 = nc.computeMu2();

    std::vector<double> H(smesh.nbFaces());
    std::vector<double> G(smesh.nbFaces());
    for (auto f = 0; f < smesh.nbFaces(); ++f) {
        const auto b = smesh.faceCentroid(f);
        const auto area = mu0.measure(b, R, f);
        H[f] = NC::meanCurvature(area, mu1.measure(b, R, f));
        G[f] = NC::GaussianCurvature(area, mu2.measure(b, R, f));
    }

    auto mmm = std::minmax_element(H.begin(), H.end());
    auto mmg = std::minmax_element(G.begin(), G.end());

    trace.info() << "Expected mean curvatures: min=0.25 max=0.625\n"
                << "Computed mean curvatures: min=" << *mmm.first << " max=" << *mmm.second << "\n"
                << "Expected Gaussian curvatures: min=-0.5 max=0.25\n"
                << "Computed mean curvatures: min=" << *mmg.first << " max=" << *mmg.second << "\n";

    trace.info() << "\n\n";


    auto muXY = nc.computeMuXY();

    std::vector<double> K1(smesh.nbFaces());
    std::vector<double> K2(smesh.nbFaces());
    std::vector<Z3i::RealVector> D1(smesh.nbFaces());
    std::vector<Z3i::RealVector> D2(smesh.nbFaces());
    smesh.computeFaceNormalsFromPositions();
    for (auto f = 0; f < smesh.nbFaces(); ++f) {
        const auto b = smesh.faceCentroid(f);
        const auto N = smesh.faceNormals()[f];
        const auto area = mu0.measure(b,R,f);
        const auto M = muXY.measure(b,R,f);
        std::tie(K1[f], K2[f], D1[f], D2[f]) = NC::principalCurvatures(area, M, N);
    }

    auto mmk1 = std::minmax_element(K1.begin(), K1.end());
    auto mmk2 = std::minmax_element(K2.begin(), K2.end());

    trace.info() << "Expected k1 curvatures: min=-0.5 max=0.25\n"
                << "Computed k1 curvatures: min="<< *mmk1.first <<" max="<< *mmk1.second <<"\n"
                << "Expected k2 curvatures: min=1 max=1\n"
                << "Computed k2 curvatures: min="<< *mmk2.first <<" max="<< *mmk2.second <<"\n";

    const auto colormapH = makeQuantifiedColorMap(makeColorMap(-0.625, 0.625));
    const auto colormapG = makeQuantifiedColorMap(makeColorMap(-0.625, 0.625));
    auto colorsH = SMW::Colors(smesh.nbFaces());
    auto colorsG = SMW::Colors(smesh.nbFaces());
    for ( auto i = 0; i < smesh.nbFaces(); i++) {
        colorsH[i] = colormapH(H[i]);
        colorsG[i] = colormapG(G[i]);
    }

    SMW::writeOBJ("ex-nc-H", smesh, colorsH);
    SMW::writeOBJ("ex-nc-G", smesh, colorsG);

}

PolyMesh* registerSurface(const CountedPtr<SH3::SurfaceMesh>& primalSurface, std::string name) {
    std::vector<std::vector<size_t>> faces;
    std::vector<Z3i::RealPoint> positions;

    for (auto f = 0; f < primalSurface->nbFaces(); ++f) {
        faces.push_back(primalSurface->incidentVertices(f));
    }
    positions = primalSurface->positions();
    return polyscope::registerSurfaceMesh(std::move(name), positions, faces);
}

PolyMesh* registerDual(const CountedPtr<SH3::SurfaceMesh>& primalSurface, std::string name) {
    std::vector<std::vector<size_t>> faces;
    std::vector<Z3i::RealPoint> positions;

    for (auto f = 0; f < primalSurface->nbFaces(); ++f) {
        faces.push_back(primalSurface->incidentVertices(f));
    }
    positions = primalSurface->positions();
    return polyscope::registerSurfaceMesh(std::move(name), positions, faces);
}


int main(int argc, char** argv)
{
    polyscope::init();

    auto params = SH3::defaultParameters() | SHG3::defaultParameters();
    std::string filename = "../DGtalObjects/bunny66.vol";
    auto binImage = SH3::makeBinaryImage(filename, params);
    auto K = SH3::getKSpace(binImage);
    auto surface = SH3::makeDigitalSurface(binImage, K, params);
    auto primalSurface = SH3::makePrimalSurfaceMesh(surface);
    /*auto polyBunny =*/ registerSurface(primalSurface, "bunny");
    // Create the dual surface mesh and register it
    /*auto polyDualBunny =*/ registerDual(primalSurface, "dual bunny");
    polyscope::show();
    return 0;
}