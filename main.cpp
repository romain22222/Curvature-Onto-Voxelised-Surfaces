#include "core.cpp"

int main(int argc, char** argv)
{
    polyscope::init();

    auto params = SH3::defaultParameters() | SHG3::defaultParameters();
    std::string filename = argc > 1 ? argv[1] : "../DGtalObjects/bunny66.vol";
    double radius = argc > 2 ? std::atof( argv[2] ) : 10.0;

    auto distribType = argc > 3 ? argToDistribType(argv[3]) : DistributionType::Exponential;

    auto binImage = SH3::makeBinaryImage(filename, params);
    auto K = SH3::getKSpace(binImage);
    auto surface = SH3::makeDigitalSurface(binImage, K, params);
    auto primalSurface = *SH3::makePrimalSurfaceMesh(surface);

    auto polyBunny = registerSurface(primalSurface, "bunny");

    for (auto m: {Method::TrivialNormalFaceCentroid, Method::DualNormalVertexPosition, Method::CorrectedNormalFaceCentroid}) {
        auto varifolds = computeVarifoldsV2(binImage, surface, radius, distribType, m);

        auto nbElements = m == Method::DualNormalVertexPosition ? primalSurface.nbVertices() : primalSurface.nbFaces();

        std::vector<RealVector> lcs;
        for (auto i = 0; i < nbElements; i++) {
            lcs.push_back(varifolds[i].curvature);
        }
        if (m == Method::DualNormalVertexPosition) {
            polyBunny->addVertexVectorQuantity(methodToString(m) + " Local Curvatures", lcs);
        } else {
            polyBunny->addFaceVectorQuantity(methodToString(m) + " Local Curvatures", lcs);
        }

        const auto lcsNorm = computeSignedNorms(primalSurface, varifolds, m);


        if (m == Method::DualNormalVertexPosition) {
            polyBunny->addVertexScalarQuantity(methodToString(m) + " Local Curvatures Norm", lcsNorm);
        } else {
            polyBunny->addFaceScalarQuantity(methodToString(m) + " Local Curvatures Norm", lcsNorm);
        }
    }

    polyscope::show();
    return 0;
}
