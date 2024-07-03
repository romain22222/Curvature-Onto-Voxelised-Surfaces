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
    Linear,
    Polynomial,
    Exponential,
    CNCLike
} DistributionType;

class RadialDistance {
public:
    RadialDistance(): center(0,0,0), radius(1), modifier(5.0) {};
    RadialDistance(const RealPoint& center, const double radius, const DistributionType& distribution, const double modifier = 5.0)
            : center(center), radius(radius), modifier(modifier) {
        switch (distribution) {
            case DistributionType::Exponential:
                // TODO : comprendre pourquoi le modifier merde lors du calcul de courbure, MAIS PAS quand on affiche le kernel
                measureFunction = [](double dRatio) {
                    return exp(-5.0/(1-dRatio*dRatio));
                };
                measureFunctionDerivate = [](double dRatio) {
                    double d = 1-dRatio*dRatio;
                    return -2*5.0*dRatio*exp(-5.0/d)/(d*d);
                };
                break;
            case DistributionType::Linear:
                measureFunction = [](double dRatio) {
                    return (1-dRatio);
                };
                measureFunctionDerivate = [](double dRatio) {
                    return -1.0;
                };
                break;
            case DistributionType::Polynomial:
                measureFunction = [](double dRatio) {
                    return (1-dRatio*dRatio)/(M_PI * 2);
                };
                measureFunctionDerivate = [](double dRatio) {
                    return -dRatio/(M_PI);
                };
                break;
            case DistributionType::CNCLike:
                measureFunction = [](double dRatio) {
                    return (1.-tanh(50.*(dRatio-0.9)))/2.;
                };
                measureFunctionDerivate = [](double dRatio) {
                    return 25.*tanh(50.*(dRatio-0.9))*tanh(50.*(dRatio-0.9)) - 25.;
                };
                break;
        }
    }
    RealPoint center;
    double radius;
    double modifier;
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

std::vector<Varifold> computeVarifoldsFromPositionsAndNormals(
        const SH3::RealPoints& positions,
        const SH3::RealVectors& normals,
        const double cRadius,
        const DistributionType cDistribType,
        const double modifier) {
    auto kdTree = LinearKDTree<RealPoint, 3>(positions);
    std::vector<size_t> indices;
    RealVector tmpSumTop;
    double tmpSumBottom;
    RadialDistance rd;
    std::vector<std::pair<double, double>> weights;
    RealVector tmpVector;
    std::vector<RealVector> curvatures (positions.size());

    for (auto f = 0; f < positions.size(); ++f) {
        tmpSumTop = RealVector();
        tmpSumBottom = 0;
        const auto b = kdTree.position(f);
        rd = RadialDistance(b, cRadius, cDistribType, modifier);
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
        curvatures[f] = -tmpSumTop/(tmpSumBottom*cRadius);
    }
    std::vector<Varifold> varifolds;
    for (auto f = 0; f < positions.size(); ++f) {
        varifolds.emplace_back(positions[f], normals[f], .5*curvatures[f]);
    }
    return varifolds;
}

std::vector<Varifold> computeVarifoldsV2(const CountedPtr<SH3::BinaryImage>& bimage, const CountedPtr<SH3::DigitalSurface>& surface, const double cRadius, const DistributionType cDistribType, const Method method, const double gridStep = 1.0, const double modifier = 5.0, const Parameters& params = SHG3::defaultParameters(), SH3::RealVectors normals = std::vector<RealVector>()) {
    auto ps = *SH3::makePrimalSurfaceMesh(surface);
    SH3::RealPoints positions;
    switch (method) {
        case TrivialNormalFaceCentroid:
        case CorrectedNormalFaceCentroid:
            for (auto f = 0; f < ps.nbFaces(); ++f) {
                positions.push_back(ps.faceCentroid(f) * gridStep);
            }
            break;
        case DualNormalVertexPosition:
            for (auto v = 0; v < ps.nbVertices(); ++v) {
                positions.push_back(ps.position(v) * gridStep);
            }
            break;
        default:
            break;
    }
    if (normals.empty()) {
        std::cout << "Computing normals" << std::endl;
        switch (method) {
            case TrivialNormalFaceCentroid:
                ps.computeFaceNormalsFromPositions();
                normals = ps.faceNormals();
                break;
            case DualNormalVertexPosition:
                ps.computeFaceNormalsFromPositions();
                ps.computeVertexNormalsFromFaceNormals();
                for (auto v = 0; v < ps.nbVertices(); ++v) {
                    normals.push_back(ps.vertexNormal(v));
                }
                break;
            case CorrectedNormalFaceCentroid:
                normals = SHG3::getIINormalVectors(bimage, SH3::getSurfelRange(surface), params);
                break;
            default:
                break;
        }
    } else {
        std::cout << "Using normals from input" << std::endl;
    }
    auto varifolds = computeVarifoldsFromPositionsAndNormals(positions, normals, cRadius, cDistribType, modifier);
    if (method == CorrectedNormalFaceCentroid) {
        for (int i = 0; i < varifolds.size(); ++i) {
            varifolds[i].curvature = normals[i] * dotProduct(varifolds[i].curvature,normals[i]) / normals[i].squaredNorm();
        }
    }
    return varifolds;
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
    if (arg == "e") {
        return DistributionType::Exponential;
    } else if (arg == "l") {
        return DistributionType::Linear;
    } else if (arg == "p") {
        return DistributionType::Polynomial;
    }
    return DistributionType::CNCLike;
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