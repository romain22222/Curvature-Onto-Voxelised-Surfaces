#include "core.cpp"

int main(int argc, char** argv)
{
    polyscope::init();

    auto params = SH3::defaultParameters() | SHG3::defaultParameters();
    std::string filename = argc > 1 ? argv[1] : "../DGtalObjects/bunny66.vol";
    double radius = argc > 2 ? std::atof( argv[2] ) : 10.0;

    auto distribType = argc > 3 ? argToDistribType(argv[3]) : DistributionType::HalfSphere;

    auto binImage = SH3::makeBinaryImage(filename, params);
    auto K = SH3::getKSpace(binImage);
    auto surface = SH3::makeDigitalSurface(binImage, K, params);
    auto primalSurface = *SH3::makePrimalSurfaceMesh(surface);

    auto polyBunny = registerSurface(primalSurface, "bunny");

    for (auto m: {Method::TrivialNormalFaceCentroid, Method::DualNormalVertexPosition, Method::CorrectedNormalFaceCentroid}) {
        auto varifolds = computeVarifolds(binImage, surface, radius, distribType, m);

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


        auto minmax = std::minmax_element(lcsNorm.begin(), lcsNorm.end());
        DGtal::trace.info() << "Min: " << *minmax.first << " Max: " << *minmax.second << std::endl;
        const auto colormap = makeColorMap(*minmax.first, *minmax.second);
        std::vector<std::vector<double>> colorLcsNorm;
        for (auto i = 0; i < nbElements; i++) {
            const auto color = lcsNorm[i] < 0 ? colormap.first(lcsNorm[i]) : colormap.second(lcsNorm[i]);
            colorLcsNorm.push_back({static_cast<double>(color.red())/255, static_cast<double>(color.green())/255, static_cast<double>(color.blue())/255});
        }
        if (m == Method::DualNormalVertexPosition) {
            polyBunny->addVertexColorQuantity(methodToString(m) + " Local Curvatures Norm", colorLcsNorm);
        } else {
            polyBunny->addFaceColorQuantity(methodToString(m) + " Local Curvatures Norm", colorLcsNorm);
        }
    }

    polyscope::show();
    return 0;
}
