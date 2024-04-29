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
    gcm.addColor(Color(255,255,255));
    gcm.addColor(Color(255,0,0));
    gcm.addColor(Color(0,0,0));
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
    Varifold(const RealPoint& position, const RealVector& planeNormal, const RealVector& curvature)
        : position(position), planeNormal(planeNormal), curvature(curvature) {
    }

    RealPoint position;
    RealVector planeNormal;
    RealVector curvature;
};

std::pair<RealVector, RealVector> computeDirectionalPlaneFromNormal(const RealVector &N) {
    const auto e1 = RealVector(1, 0, 0);
    if (N.crossProduct(e1).norm() == 0) {
        return std::make_pair(RealVector(0, 1, 0), RealVector(0, 0, 1));
    }
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


typedef enum {
    FlatDisc,
    Cone,
    HalfSphere
} DistributionType;

class RadialDistance {
public:
    RadialDistance(): center(0,0,0), radius(1) {};
    RadialDistance(const RealPoint& center, const double radius, const DistributionType& distribution)
        : center(center), radius(radius) {
        switch (distribution) {
            case DistributionType::FlatDisc:
                measureFunction = [](double dRatio) {
                    return 3./(4*M_PI);
                };
                measureFunctionDerivate = [](double dRatio) {
                    return 0;
                };
                break;
            case DistributionType::Cone:
                measureFunction = [](double dRatio) {
                    return (1-dRatio) * M_PI/12.;
                };
                measureFunctionDerivate = [](double dRatio) {
                    return -M_PI/12.;
                };
                break;
            case DistributionType::HalfSphere:
                measureFunction = [](double dRatio) {
                    return (1-dRatio*dRatio)/(M_PI * 2);
                };
                measureFunctionDerivate = [](double dRatio) {
                    return -dRatio/(M_PI);
                };
                break;
        }
    }
    RealPoint center;
    double radius;
    std::function<double(double)> measureFunction;
    std::function<double(double)> measureFunctionDerivate;

    std::vector<std::pair<double,double>> operator()(const SH3::RealPoints& mesh) const {
        std::vector<std::pair<double,double>> wf;
        for (const auto& b : mesh) {
            // If the face is inside the radius, compute the weight
            const auto d = (b - center).norm();
            wf.push_back(d < radius ? std::make_pair(measureFunction(d/radius), measureFunctionDerivate(d/radius)) : std::make_pair(0., 0.));
        }
        return wf;
    }
};

typedef enum {
    TrivialNormalFaceCentroid,
    DualNormalFaceCentroid,
    CorrectedNormalFaceCentroid,
    ProbabilisticOfTrivials,
    VertexInterpolation
} Method;

RealVector projection(const RealVector& toProject, const RealVector& planeNormal) {
    return toProject - planeNormal * (toProject.dot(planeNormal)/planeNormal.squaredNorm());
}

std::vector<RealVector> computeLocalCurvature(const CountedPtr<SH3::BinaryImage>& bimage, const CountedPtr<SH3::DigitalSurface>& surface, const double cRadius, const DistributionType cDistribType, const Method method) {
    std::vector<RealVector> curvatures;
    const CountedPtr<SH3::SurfaceMesh> pSurface = SH3::makePrimalSurfaceMesh(surface);
    RadialDistance rd;
    std::vector<std::pair<double, double>> weights;
    RealVector tmpSumTop;
    double tmpSumBottom;
    RealVector tmpVector;
    int percent = 0;

    auto positions = SH3::RealPoints();

    auto normals = SH3::RealVectors();

    unsigned long nbElements;

    pSurface->computeFaceNormalsFromPositions();
    pSurface->computeVertexNormalsFromFaceNormals();

    switch (method) {
        case TrivialNormalFaceCentroid:
            nbElements = pSurface->nbFaces();
            for (auto f = 0; f < nbElements; ++f) {
                positions.push_back(pSurface->faceCentroid(f));
                normals.push_back(pSurface->faceNormal(f));
            }
            break;
        case DualNormalFaceCentroid:
            nbElements = pSurface->nbVertices();
            for (auto v = 0; v < nbElements; ++v) {
                positions.push_back(pSurface->position(v));
                normals.push_back(pSurface->vertexNormal(v));
            }
            break;
        case CorrectedNormalFaceCentroid:
            nbElements = pSurface->nbFaces();
            normals = SHG3::getIINormalVectors(bimage, SH3::getSurfelRange(surface), SHG3::defaultParameters());
            for (auto f = 0; f < nbElements; ++f) {
                positions.push_back(pSurface->faceCentroid(f));
            }
            break;
        default:
            return curvatures;
    }

    for (auto f = 0; f < nbElements; ++f) {
        if(f*100/nbElements > percent) {
            percent += 1;
            DGtal::trace.info() << "Computing curvatures: " << percent << "%\n";
        }
        tmpSumTop = RealVector();
        tmpSumBottom = 0;
        const auto b = positions[f];
        rd = RadialDistance(b, cRadius, cDistribType);
        weights = rd(positions);
        for (auto otherF = 0; otherF < nbElements; ++otherF) {
            if (weights[otherF].first > 0) {
                // Add to tmpSum the awaited result
                if (f != otherF) {
                    tmpVector = positions[otherF] - b;
                    tmpSumTop += weights[otherF].first * projection(tmpVector, normals[otherF])/tmpVector.norm();
                }
                tmpSumBottom += weights[otherF].first;
            }
        }
        curvatures.push_back(-tmpSumTop/(tmpSumBottom*cRadius));
    }

    return curvatures;
}

std::vector<Varifold> computeVarifolds(const CountedPtr<SH3::BinaryImage>& bimage, const CountedPtr<SH3::DigitalSurface>& surface, const double cRadius, const DistributionType cDistribType, const Method method) {
    std::vector<Varifold> varifolds;

    auto ps = *SH3::makePrimalSurfaceMesh(surface);

    int percent = 0;

    auto curvatures = computeLocalCurvature(bimage, surface, cRadius, cDistribType, method);

    switch (method) {
        case TrivialNormalFaceCentroid:
        case CorrectedNormalFaceCentroid:
            ps.computeFaceNormalsFromPositions();
            for (auto f = 0; f < ps.nbFaces(); ++f) {
                if(f*100/ps.nbFaces() > percent) {
                    percent += 1;
                    DGtal::trace.info() << "Computing varifolds: " << percent << "%\n";
                }
                varifolds.emplace_back(ps.faceCentroid(f), ps.faceNormal(f), curvatures[f]);
            }
            break;
        case DualNormalFaceCentroid:
            ps.computeFaceNormalsFromPositions();
            ps.computeVertexNormalsFromFaceNormals();
            for (auto v = 0; v < ps.nbVertices(); ++v) {
                if(v*100/ps.nbVertices() > percent) {
                    percent += 1;
                    DGtal::trace.info() << "Computing varifolds: " << percent << "%\n";
                }
                varifolds.emplace_back(ps.position(v), ps.vertexNormal(v), curvatures[v]);
            }
            break;
        default:
            break;
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
    auto varifolds = computeVarifolds(binImage, surface, 10, DistributionType::HalfSphere, Method::TrivialNormalFaceCentroid);
    DGtal::trace.info() << "Computed " << varifolds.size() << " varifolds" << std::endl;

    std::vector<RealVector> normals;
    std::vector<RealVector> lcs;
    for (auto i = 0; i < primalSurface.nbFaces(); i++) {
        normals.push_back(varifolds[i].planeNormal);
        lcs.push_back(varifolds[i].curvature);
    }
//    polyBunny->addFaceVectorQuantity("Trivial Normals", normals);
    polyBunny->addFaceVectorQuantity("TNFC Local Curvatures", lcs);

    // Create a heatmap of the norm of the local curvatures
    std::vector<double> lcsNorm;
    for (auto i = 0; i < primalSurface.nbFaces(); i++) {
        lcsNorm.push_back(lcs[i].norm());
    }

    auto minmax = std::minmax_element(lcsNorm.begin(), lcsNorm.end());
    const auto colormap = makeColorMap(*minmax.first, *minmax.second);
    std::vector<std::vector<double>> colorLcsNorm;
    for (auto i = 0; i < primalSurface.nbFaces(); i++) {
        const auto color = colormap(lcsNorm[i]);
        colorLcsNorm.push_back({static_cast<double>(color.red())/255, static_cast<double>(color.green())/255, static_cast<double>(color.blue())/255});
    }
    polyBunny->addFaceColorQuantity("TNFC Local Curvatures Norm", colorLcsNorm);

    auto varifolds2 = computeVarifolds(binImage, surface, 10, DistributionType::HalfSphere, Method::DualNormalFaceCentroid);

    std::vector<RealVector> normals2;
    std::vector<RealVector> lcs2;

    for (auto i = 0; i < primalSurface.nbVertices(); i++) {
        normals2.push_back(varifolds2[i].planeNormal);
        lcs2.push_back(varifolds2[i].curvature);
    }
//    polyBunny->addVertexVectorQuantity("Dual Normals", normals2);
    polyBunny->addVertexVectorQuantity("DNFC Local Curvatures", lcs2);

    // Create a heatmap of the norm of the local curvatures
    std::vector<double> lcsNorm2;
    for (auto i = 0; i < primalSurface.nbVertices(); i++) {
        lcsNorm2.push_back(lcs2[i].norm());
    }

    auto minmax2 = std::minmax_element(lcsNorm2.begin(), lcsNorm2.end());
    const auto colormap2 = makeColorMap(*minmax2.first, *minmax2.second);
    std::vector<std::vector<double>> colorLcsNorm2;
    for (auto i = 0; i < primalSurface.nbVertices(); i++) {
        const auto color = colormap2(lcsNorm2[i]);
        colorLcsNorm2.push_back({static_cast<double>(color.red())/255, static_cast<double>(color.green())/255, static_cast<double>(color.blue())/255});
    }
    polyBunny->addVertexColorQuantity("DNFC Local Curvatures Norm", colorLcsNorm2);

    auto varifolds3 = computeVarifolds(binImage, surface, 10, DistributionType::HalfSphere, Method::CorrectedNormalFaceCentroid);

    std::vector<RealVector> normals3;
    std::vector<RealVector> lcs3;

    for (auto i = 0; i < primalSurface.nbFaces(); i++) {
        normals3.push_back(varifolds3[i].planeNormal);
        lcs3.push_back(varifolds3[i].curvature);
    }
//    polyBunny->addFaceVectorQuantity("Corrected Normals", normals3);
    polyBunny->addFaceVectorQuantity("CNFC Local Curvatures", lcs3);

    // Create a heatmap of the norm of the local curvatures
    std::vector<double> lcsNorm3;
    for (auto i = 0; i < primalSurface.nbFaces(); i++) {
        lcsNorm3.push_back(lcs3[i].norm());
    }

    auto minmax3 = std::minmax_element(lcsNorm3.begin(), lcsNorm3.end());
    const auto colormap3 = makeColorMap(*minmax3.first, *minmax3.second);
    std::vector<std::vector<double>> colorLcsNorm3;
    for (auto i = 0; i < primalSurface.nbFaces(); i++) {
        const auto color = colormap3(lcsNorm3[i]);
        colorLcsNorm3.push_back({static_cast<double>(color.red())/255, static_cast<double>(color.green())/255, static_cast<double>(color.blue())/255});
    }
    polyBunny->addFaceColorQuantity("CNFC Local Curvatures Norm", colorLcsNorm3);

    polyscope::show();
    return 0;
}