#include <utility>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "DGtal/io/writers/SurfaceMeshWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/QuantifiedColorMap.h"
#include "DGtal/shapes/SurfaceMeshHelper.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "externalLibs/LinearKDTree.h"

using namespace DGtal;
using namespace DGtal::Z3i;
typedef Shortcuts<KSpace>         SH3;
typedef ShortcutsGeometry<KSpace> SHG3;
typedef polyscope::SurfaceMesh PolyMesh;

std::pair<GradientColorMap<double>, GradientColorMap<double>> makeColorMap(double minv, double maxv) {
    if (maxv < 0) {
        GradientColorMap<double> gcm(minv, maxv);
        gcm.addColor(Color(0,0,255));
        gcm.addColor(Color(255,255,255));
        return {gcm, gcm};
    }
    if (minv > 0) {
        GradientColorMap<double> gcm(minv, maxv);
        gcm.addColor(Color(255,255,255));
        gcm.addColor(Color(255,0,0));
        gcm.addColor(Color(0,0,0));
        return {gcm, gcm};
    }
    GradientColorMap<double> gcm(minv, 0.);
    gcm.addColor(Color(0,0,255));
    gcm.addColor(Color(255,255,255));
    GradientColorMap<double> gcm2(0., maxv);
    gcm2.addColor(Color(255,255,255));
    gcm2.addColor(Color(255,0,0));
    gcm2.addColor(Color(0,0,0));
    return {gcm, gcm2};
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

class Varifold {
public:
    Varifold(const RealPoint& position, const RealVector& planeNormal, const RealVector& curvature)
            : position(position), planeNormal(planeNormal), curvature(curvature) {
    }

    RealPoint position;
    RealVector planeNormal;
    RealVector curvature;
};

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
                    return (1-dRatio) * M_PI/3.;
                };
                measureFunctionDerivate = [](double dRatio) {
                    return -M_PI/3.;
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

    std::vector<std::pair<double,double>> operator()(const SH3::RealPoints& mesh, const std::vector<size_t>& poi) const {
        std::vector<std::pair<double,double>> wf;
        for (const auto& b : poi) {
            // If the face is inside the radius, compute the weight
            const auto d = (mesh[b] - center).norm();
            if (d < radius) {
                wf.emplace_back(measureFunction(d / radius), measureFunctionDerivate(d / radius));
            } else {
                wf.emplace_back(0., 0.);
            }
        }
        return wf;
    }
};

typedef enum {
    TrivialNormalFaceCentroid,
    DualNormalVertexPosition,
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
        case DualNormalVertexPosition:
            nbElements = pSurface->nbVertices();
            for (auto v = 0; v < nbElements; ++v) {
                positions.push_back(pSurface->position(v));
                normals.push_back(pSurface->vertexNormal(v));
            }
            break;
        case CorrectedNormalFaceCentroid:
            nbElements = pSurface->nbFaces();
            normals = SHG3::getIINormalVectors(bimage, SH3::getSurfelRange(surface), SHG3::defaultParameters()("verbose", 0));
            for (auto f = 0; f < nbElements; ++f) {
                positions.push_back(pSurface->faceCentroid(f));
            }
            break;
        default:
            return curvatures;
    }

    auto kdTree = LinearKDTree<RealPoint, 3>(positions);
    std::vector<size_t> indices;
    for (auto f = 0; f < nbElements; ++f) {
        tmpSumTop = RealVector();
        tmpSumBottom = 0;
        const auto b = kdTree.position(f);
        rd = RadialDistance(b, cRadius, cDistribType);
        indices = kdTree.pointsInBall(b, cRadius);
        weights = rd(positions, indices);
        for (auto otherF = 0; otherF < weights.size(); ++otherF) {
            if (weights[otherF].first > 0) {
                if (f != indices[otherF]) {
                    tmpVector = positions[indices[otherF]] - b;
                    tmpSumTop += weights[otherF].second * projection(tmpVector, normals[indices[otherF]])/tmpVector.norm();
                }
                tmpSumBottom += weights[otherF].first;
            }
        }
        curvatures.push_back(-tmpSumTop/(tmpSumBottom*cRadius));
    }

    return curvatures;
}

std::vector<Varifold> computeVarifolds(const CountedPtr<SH3::BinaryImage>& bimage, const CountedPtr<SH3::DigitalSurface>& surface, const double cRadius, const DistributionType cDistribType, const Method method, const double gridStep = 1.0) {
    std::vector<Varifold> varifolds;

    auto ps = *SH3::makePrimalSurfaceMesh(surface);

    auto curvatures = computeLocalCurvature(bimage, surface, cRadius, cDistribType, method);

    SH3::RealVectors normals;

    switch (method) {
        case TrivialNormalFaceCentroid:
        case CorrectedNormalFaceCentroid:
            ps.computeFaceNormalsFromPositions();
            if (method == Method::CorrectedNormalFaceCentroid) {
                normals = SHG3::getIINormalVectors(bimage, SH3::getSurfelRange(surface), SHG3::defaultParameters()("verbose", 0));
            } else {
                ps.computeVertexNormalsFromFaceNormals();
                normals = ps.faceNormals();
            }

            for (auto f = 0; f < ps.nbFaces(); ++f) {
                varifolds.emplace_back(ps.faceCentroid(f), normals[f], curvatures[f]);
            }
            break;
        case DualNormalVertexPosition:
            ps.computeFaceNormalsFromPositions();
            ps.computeVertexNormalsFromFaceNormals();
            for (auto v = 0; v < ps.nbVertices(); ++v) {
                varifolds.emplace_back(ps.position(v), ps.vertexNormal(v), curvatures[v]);
            }
            break;
        default:
            break;
    }
    std::transform(varifolds.begin(), varifolds.end(), varifolds.begin(), [&gridStep](const Varifold& v) {
        return Varifold(v.position, v.planeNormal, 0.5*v.curvature / gridStep);
    });
    return varifolds;
}

DistributionType argToDistribType(const std::string& arg) {
    if (arg == "fd") {
        return DistributionType::FlatDisc;
    } else if (arg == "c") {
        return DistributionType::Cone;
    } else {
        return DistributionType::HalfSphere;
    }
}

Method argToMethod(const std::string& arg) {
    if (arg == "tnfc") {
        return Method::TrivialNormalFaceCentroid;
    } else if (arg == "dnvp") {
        return Method::DualNormalVertexPosition;
    } else if (arg == "cnfc") {
        return Method::CorrectedNormalFaceCentroid;
    } else if (arg == "pot") {
        return Method::ProbabilisticOfTrivials;
    } else {
        return Method::VertexInterpolation;
    }
}

std::string methodToString(const Method& method) {
    switch (method) {
        case TrivialNormalFaceCentroid:
            return "Trivial Normal Face Centroid";
        case DualNormalVertexPosition:
            return "Dual Normal Face Centroid";
        case CorrectedNormalFaceCentroid:
            return "Corrected Normal Face Centroid";
        case ProbabilisticOfTrivials:
            return "Probabilistic Of Trivials";
        case VertexInterpolation:
            return "Vertex Interpolation";
        default:
            return "Unknown";
    }
}


std::vector<double> computeSignedNorms(const SH3::SurfaceMesh& primalSurface, const std::vector<Varifold>& varifolds, const Method& m)
{
    std::vector<double> lcsNorm;
    for (const auto & varifold : varifolds) {
        lcsNorm.push_back(varifold.planeNormal.dot(varifold.curvature) > 0 ? -varifold.curvature.norm() : varifold.curvature.norm());
    }
    if (m == Method::DualNormalVertexPosition) {
        for (auto i = 0; i < varifolds.size(); i++) {
            const auto& position = primalSurface.position(i);
            auto sum = 0.;
            for (auto f = 0; f < varifolds.size(); f++) {
                if (f != i && primalSurface.vertexInclusionRatio(position, 1, f) > 0) {
                    sum += lcsNorm[f];
                }
            }
            lcsNorm[i] = abs(lcsNorm[i]) * (sum < 0 ? -1 : 1);
        }
    } else {
        for (auto i = 0; i < varifolds.size(); i++) {
            auto sum = 0.;
            for (auto f: primalSurface.computeFacesInclusionsInBall(1, i)) {
                if (f.second > 0) {
                    sum += lcsNorm[f.first];
                }
            }
            lcsNorm[i] = abs(lcsNorm[i]) * (sum < 0 ? -1 : 1);
        }
    }
    return lcsNorm;
}