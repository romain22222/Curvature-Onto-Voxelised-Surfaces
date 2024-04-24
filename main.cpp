#include <utility>

#include "DGtal/base/Common.h"
#include "DGtal/geometry/meshes/NormalCycleComputer.h"
#include "DGtal/geometry/meshes/CorrectedNormalCurrentComputer.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "DGtal/io/writers/SurfaceMeshWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/QuantifiedColorMap.h"
#include "DGtal/shapes/SurfaceMeshHelper.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace DGtal;
using namespace DGtal::Z3i;
typedef SurfaceMesh< RealPoint, RealVector > SM;
typedef NormalCycleComputer< RealPoint, RealVector > NC;
typedef CorrectedNormalCurrentComputer< RealPoint, RealVector > CNC;
typedef SurfaceMeshHelper< RealPoint, RealVector > SMH;
typedef SurfaceMeshWriter< RealPoint, RealVector > SMW;
typedef Shortcuts<KSpace>         SH3;
typedef ShortcutsGeometry<KSpace> SHG3;
typedef polyscope::SurfaceMesh PolyMesh;
#define R 0.5

#define MAX_STEPS_L2_MEASURE 5

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
    SM smesh = SMH::makeTorus(3.0, 1.0, RealPoint(), 20, 20, 0, SMH::NormalsType::VERTEX_NORMALS);
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
    std::vector<RealVector> D1(smesh.nbFaces());
    std::vector<RealVector> D2(smesh.nbFaces());
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

PolyMesh* registerSurface(const SH3::SurfaceMesh& surface, std::string name) {
    std::vector<std::vector<size_t>> faces;
    std::vector<RealPoint> positions;
    auto smP = CountedPtr<SH3::SurfaceMesh>(new SH3::SurfaceMesh(surface));

    for (auto f = 0; f < surface.nbFaces(); ++f) {
        faces.push_back(surface.incidentVertices(f));
    }
    positions = surface.positions();
    return polyscope::registerSurfaceMesh(std::move(name), positions, faces);
}

PolyMesh* registerDual(PolygonalSurface<RealPoint> surface, std::string name) {
    std::vector<std::vector<size_t>> faces;

    for (auto f = 0; f < surface.nbFaces(); ++f) {
        faces.push_back(surface.verticesAroundFace(f));
    }
    auto positions = surface.positions();
    return polyscope::registerSurfaceMesh(std::move(name), positions, faces);
}

class Varifold {
public:
    Varifold(const RealPoint& position, const RealVector& dirPlaneX, const RealVector& dirPlaneY)
        : position(position), dirPlaneX(dirPlaneX), dirPlaneY(dirPlaneY) {}
    RealPoint position;
    RealVector dirPlaneX;
    RealVector dirPlaneY;
};

std::pair<RealVector, RealVector> computeDirectionalPlaneFromNormal(const RealVector &N) {
    const auto e1 = RealVector(1, 0, 0);
    const auto v = N.crossProduct(e1)/N.crossProduct(e1).norm();
    const auto w = N.crossProduct(v);
    return std::make_pair(v, w);
}

double oldComputeL2Measure(SH3::SurfaceMesh& mesh, int f) {
    // 1 - Create a queue of faces to process, each element is a pair of face index and the previous step + 1
    std::queue<std::pair<unsigned long, int>> queue;
    queue.emplace(f, 0);
    // 2 - Create a set of visited faces
    std::set<unsigned long> visited;
    visited.insert(f);
    // On each step, transfer one face from the queue to the visited set, and add its neighbors to the queue
    double l2measure = 0;
    while (!queue.empty()) {
        auto top = queue.front();
        queue.pop();
        l2measure += 1;
        if (top.second < MAX_STEPS_L2_MEASURE) {
            for (auto neighbor : mesh.neighborFaces(top.first)) {
                if (visited.find(neighbor) == visited.end()) {
                    queue.emplace(neighbor, top.second + 1);
                    visited.insert(neighbor);
                }
            }
        }
    }
    return l2measure;
}


typedef std::pair<int, double> WeightedFace;

typedef enum {
    FlatDisc,
    Cone,
    HalfSphere
} DistributionType;

class RadialDistance {
public:
    RadialDistance(const RealPoint& center, const double radius, const DistributionType& distribution)
        : center(center), radius(radius) {
        switch (distribution) {
            case DistributionType::FlatDisc:
                measureFunction = [&radius](const SH3::SurfaceMesh& mesh, int face, double distance) {
                    // Not sure about if we should weight the radius here or not
                    return WeightedFace(face, 3./(4*M_1_PI));
//                    return WeightedFace(face, 3./(4*M_1_PI*radius*radius*radius));
                };
                break;
            case DistributionType::Cone:
                measureFunction = [&radius](const SH3::SurfaceMesh& mesh, int face, double distance) {
                    return WeightedFace(face, 1.); // TODO: Implement the cone distribution
                };
                break;
            case DistributionType::HalfSphere:
                measureFunction = [&radius](const SH3::SurfaceMesh& mesh, int face, double distance) {
                    return WeightedFace(face, 1.); // TODO: Implement the half sphere distribution
                };
                break;
        }
    }
    RealPoint center;
    double radius;
    std::function<WeightedFace(const SH3::SurfaceMesh&, int face, double distance)> measureFunction;

    std::vector<WeightedFace> operator()(const SH3::SurfaceMesh& mesh) const {
        std::vector<WeightedFace> wf;
        for (auto f = 0; f < mesh.nbFaces(); ++f) {
            // If the face is inside the radius, compute the weight
            const auto b = mesh.faceCentroid(f);
            const auto d = (b - center).norm();
            wf.push_back(d < radius ? measureFunction(mesh, f, d) : WeightedFace(f, 0));
        }
        return wf;
    }
};

std::pair<RealVector, RealVector> computeDirectionalPlane(
        const SH3::SurfaceMesh& mesh,
        int f) {
    return computeDirectionalPlaneFromNormal(mesh.faceNormal(f));
}

std::vector<Varifold> computeVarifolds(SH3::SurfaceMesh& surface) {
    std::vector<Varifold> varifolds;

    surface.computeFaceNormalsFromPositions();
    const CNC cnc(surface);
    int percent = 0;
    for (auto f = 0; f < surface.nbFaces(); ++f) {
        if(f*100/surface.nbFaces() > percent) {
            percent += 1;
            DGtal::trace.info() << "Computing varifolds: " << percent << "%\n";
        }
        const auto l2measure = oldComputeL2Measure(surface, f);
        const auto dirPlane = computeDirectionalPlane(surface, f);
        varifolds.emplace_back(surface.faceCentroid(f), dirPlane.first, dirPlane.second);
    }

    return varifolds;
}


int main(int argc, char** argv)
{
    polyscope::init();

    auto params = SH3::defaultParameters() | SHG3::defaultParameters();
    std::string filename = "../DGtalObjects/bunny66.vol";
    auto binImage = SH3::makeBinaryImage(filename, params);
    auto K = SH3::getKSpace(binImage);
    auto surface = SH3::makeDigitalSurface(binImage, K, params);

    auto primalSurface = *SH3::makePrimalSurfaceMesh(surface);
//    SH3::SurfaceMesh primalSurface = SMH::makeTorus(5.0, 2.0, RealPoint(), 60, 60, 0, SMH::NormalsType::VERTEX_NORMALS);


    auto polyBunny = registerSurface(primalSurface, "bunny");

    auto dualSurface = *SH3::makeDualPolygonalSurface(surface);
    /*auto polyDualBunny =*/ registerDual(dualSurface, "dual bunny");

    // Compute a heat map of the L2 measure
    DGtal::trace.info() << primalSurface.nbFaces() << std::endl;
    auto varifolds = computeVarifolds(primalSurface);
    DGtal::trace.info() << "Computed " << varifolds.size() << " varifolds" << std::endl;

    // Create 2 vector fields on the surface, one for each direction of the principal curvatures
    std::vector<RealVector> dir1;
    std::vector<RealVector> dir2;
    for (auto i = 0; i < primalSurface.nbFaces(); i++) {
        dir1.push_back(varifolds[i].dirPlaneX);
        dir2.push_back(varifolds[i].dirPlaneY);
    }
    polyBunny->addFaceVectorQuantity("Principal Curvature 1", dir1);
    polyBunny->addFaceVectorQuantity("Principal Curvature 2", dir2);



    polyscope::show();
    return 0;
}