
#include <iostream>
#include <algorithm>
#include "DGtal/base/Common.h"
#include "DGtal/shapes/SurfaceMesh.h"
#include "DGtal/geometry/meshes/CorrectedNormalCurrentComputer.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "DGtal/io/writers/SurfaceMeshWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/QuantifiedColorMap.h"
#include "core.cpp"

void usage( int argc, char* argv[] )
{
    using namespace DGtal;
    using namespace DGtal::Z3i;
    typedef Shortcuts< KSpace >          SH;
    std::cout << "Usage: " << std::endl
              << "\t" << argv[ 0 ] << " <P> <B> <h> <R> <kernel> <method>" << std::endl
              << std::endl
              << "Computation of mean and Gaussian curvatures on an "      << std::endl
              << "digitized implicit shape using constant or "             << std::endl
              << "interpolated corrected curvature measures (based "       << std::endl
              << "on the theory of corrected normal currents)."            << std::endl
              << "- builds the surface mesh from polynomial <P>"           << std::endl
              << "- <B> defines the digitization space size [-B,B]^3"      << std::endl
              << "- <h> is the gridstep digitization"                      << std::endl
              << "- <R> is the radius of the measuring balls"              << std::endl
              << "- <kernel> is the kernel used to sample the surface ('fd': flat disc, 'c': cone, 'hs': half sphere)" << std::endl
              << "- <method> is the method used to compute the curvature ('tnfc': trivial normal face centroid, 'cnfc': corrected normal face centroid)" << std::endl
              << std::endl
              << "It produces several OBJ files to display mean and"       << std::endl
              << "Gaussian curvature estimation results: `example-cnc-H.obj`" << std::endl
              << "and `example-cnc-G.obj` as well as the associated MTL file." << std::endl;
    std::cout << "You may either write your own polynomial as 3*x^2*y-z^2*x*y+1" << std::endl
              <<"or use a predefined polynomial in the following list:" << std::endl;
    auto L = SH::getPolynomialList();
    for ( const auto& p : L )
        std::cout << p.first << " : " << p.second << std::endl;
}

int main( int argc, char* argv[] )
{
    polyscope::init();
    if ( argc <= 1 )
    {
        usage( argc, argv );
        return 0;
    }
    using namespace DGtal;
    using namespace DGtal::Z3i;
    typedef SurfaceMesh< RealPoint, RealVector >                    SM;
    typedef CorrectedNormalCurrentComputer< RealPoint, RealVector > CNC;
    typedef Shortcuts< KSpace >          SH;
    typedef ShortcutsGeometry< KSpace > SHG;
    std::string  poly = argv[ 1 ]; // polynomial
    const double    B = argc > 2 ? atof( argv[ 2 ] ) : 1.0; // max ||_oo bbox
    const double    h = argc > 3 ? atof( argv[ 3 ] ) : 1.0; // gridstep
    const double    R = argc > 4 ? atof( argv[ 4 ] ) : 2.0; // radius of measuring ball
    const auto kernel = argc > 5 ? argToDistribType( argv[ 5 ] ) : DistributionType::HalfSphere;
    const auto method = argc > 6 ? argToMethod( argv[ 6 ] ) : Method::CorrectedNormalFaceCentroid;
    const auto checkCNC = argc > 7 && std::strcmp(argv[7],"y") == 0;

    // Read polynomial and build digital surface
    auto params = SH::defaultParameters() | SHG::defaultParameters();
    params( "t-ring", 6 )( "surfaceTraversal", "Default" );
    params( "polynomial", poly )( "gridstep", h );
    params( "minAABB", -B )( "maxAABB", B );
    params( "offset", 3.0 );
    auto shape       = SH::makeImplicitShape3D( params );
    auto K           = SH::getKSpace( params );
    auto dshape      = SH::makeDigitizedImplicitShape3D( shape, params );
    auto bimage      = SH::makeBinaryImage( dshape, params );
    if ( bimage == nullptr )
    {
        trace.error() <<  "Unable to read polynomial <"
                      << poly.c_str() << ">" << std::endl;
        return 1;
    }
    auto sembedder   = SH::getSCellEmbedder( K );
    auto embedder    = SH::getCellEmbedder( K );
    auto surface     = SH::makeDigitalSurface( bimage, K, params );
    auto surfels     = SH::getSurfelRange( surface, params );
    trace.info() << "- surface has " << surfels.size()<< " surfels." << std::endl;

    SM smesh;
    std::vector< SM::Vertices > faces;
    SH::Cell2Index c2i;
    auto pointels = SH::getPointelRange( c2i, surface );
    auto vertices = SH::RealPoints( pointels.size() );
    std::transform( pointels.cbegin(), pointels.cend(), vertices.begin(),
                    [&] (const SH::Cell& c) { return h * embedder( c ); } );
    for ( auto&& surfel : *surface )
    {
        const auto primal_surfel_vtcs = SH::getPointelRange( K, surfel );
        SM::Vertices face;
        for ( auto&& primal_vtx : primal_surfel_vtcs )
            face.push_back( c2i[ primal_vtx ] );
        faces.push_back( face );
    }
    smesh.init( vertices.cbegin(), vertices.cend(),
                faces.cbegin(),    faces.cend() );
    trace.info() << smesh << std::endl;

    auto polysurf = registerSurface(smesh, "studied mesh");


    std::vector<Varifold> varifolds = computeVarifolds(bimage, surface, R, kernel, method, h);

    auto exp_H = SHG::getMeanCurvatures( shape, K, surfels, params );
    auto exp_G = SHG::getGaussianCurvatures( shape, K, surfels, params );

    std::vector< double > H = computeSignedNorms(smesh, varifolds, method);
    std::vector< double > G( varifolds.size() );

    auto H_min_max = std::minmax_element( H.cbegin(), H.cend() );
    auto G_min_max = std::minmax_element( G.cbegin(), G.cend() );
    auto exp_H_min_max = std::minmax_element( exp_H.cbegin(), exp_H.cend() );
    auto exp_G_min_max = std::minmax_element( exp_G.cbegin(), exp_G.cend() );
    std::cout << "Expected mean curvatures:"
              << " min=" << *exp_H_min_max.first << " max=" << *exp_H_min_max.second
              << std::endl;
    std::cout << "Computed mean curvatures:"
              << " min=" << *H_min_max.first << " max=" << *H_min_max.second
              << std::endl;
    std::cout << "Expected Gaussian curvatures:"
              << " min=" << *exp_G_min_max.first << " max=" << *exp_G_min_max.second
              << std::endl;
    std::cout << "Computed Gaussian curvatures:"
              << " min=" << *G_min_max.first << " max=" << *G_min_max.second
              << std::endl;

    const auto      error_H = SHG::getScalarsAbsoluteDifference( H, exp_H );
    const auto stat_error_H = SHG::getStatistic( error_H );
    const auto   error_H_l2 = SHG::getScalarsNormL2( H, exp_H );
    trace.info() << "|He-H|_oo = " << stat_error_H.max() << std::endl;
    trace.info() << "|He-H|_2  = " << error_H_l2 << std::endl;
    const auto      error_G = SHG::getScalarsAbsoluteDifference( G, exp_G );
    const auto stat_error_G = SHG::getStatistic( error_G );
    const auto   error_G_l2 = SHG::getScalarsNormL2( G, exp_G );
    trace.info() << "|Ge-G|_oo = " << stat_error_G.max() << std::endl;
    trace.info() << "|Ge-G|_2  = " << error_G_l2 << std::endl;

    // Remove normals for better blocky display.
    smesh.vertexNormals() = SH::RealVectors();
    smesh.faceNormals()   = SH::RealVectors();
    typedef SurfaceMeshWriter< RealPoint, RealVector > SMW;
    const double    Hmax = std::max( fabs( *exp_H_min_max.first ),
                                     fabs( *exp_H_min_max.second ) );
    const double    Gmax = std::max( fabs( *exp_G_min_max.first ),
                                     fabs( *exp_G_min_max.second ) );
    const auto colormapH = makeQuantifiedColorMap( makeColorMap( -Hmax, Hmax ).first );
    const auto colormapG = makeQuantifiedColorMap( makeColorMap( -Gmax, Gmax ).first );
    auto colorsH = SMW::Colors( varifolds.size() );
    auto colorsG = SMW::Colors( varifolds.size() );
    for ( auto i = 0; i < varifolds.size(); i++ )
    {
        colorsH[ i ] = colormapH( H[ i ] );
        colorsG[ i ] = colormapG( G[ i ] );
    }

    SMW::writeOBJ( "example-cnc-H", smesh, colorsH );
    SMW::writeOBJ( "example-cnc-G", smesh, colorsG );


    if (checkCNC) {
        // Builds a CorrectedNormalCurrentComputer object onto the SurfaceMesh object
        CNC cnc(smesh);
        // Estimates normal vectors using Convolved Trivial Normal estimator
        auto face_normals = SHG::getCTrivialNormalVectors(surface, surfels, params);
        // Set corrected face normals => Corrected Normal Current with
        // constant per face corrected vector field.
        smesh.setFaceNormals(face_normals.cbegin(), face_normals.cend()); // CCNC
        // computes area, mean and Gaussian curvature measures
        auto mu0 = cnc.computeMu0();
        auto mu1 = cnc.computeMu1();
        auto mu2 = cnc.computeMu2();
        // estimates mean (H) and Gaussian (G) curvatures by measure normalization.
        std::vector<double> H_CNC(varifolds.size());
        std::vector<double> G_CNC(varifolds.size());

        for (auto f = 0; f < varifolds.size(); ++f) {
            const auto b = smesh.faceCentroid(f);
            const auto area = mu0.measure(b, R, f);
            H_CNC[f] = cnc.meanCurvature(area, mu1.measure(b, R, f));
            G_CNC[ f ] = cnc.GaussianCurvature( area, mu2.measure( b, R, f ) );
        }
        auto H_CNC_min_max = std::minmax_element( H_CNC.cbegin(), H_CNC.cend() );
        auto G_CNC_min_max = std::minmax_element( G_CNC.cbegin(), G_CNC.cend() );
        std::cout << "CNC computed mean curvatures:"
                  << " min=" << *H_CNC_min_max.first << " max=" << *H_CNC_min_max.second
                  << std::endl;
        std::cout << "CNC computed Gaussian curvatures:"
                  << " min=" << *G_CNC_min_max.first << " max=" << *G_CNC_min_max.second
                  << std::endl;
        const auto      error_H_CNC = SHG::getScalarsAbsoluteDifference( H_CNC, exp_H );
        const auto stat_error_H_CNC = SHG::getStatistic( error_H_CNC );
        const auto   error_H_CNC_l2 = SHG::getScalarsNormL2( H_CNC, exp_H );
        trace.info() << "|He-H_CNC|_oo = " << stat_error_H_CNC.max() << std::endl;
        trace.info() << "|He-H_CNC|_2  = " << error_H_CNC_l2 << std::endl;
        const auto      error_G_CNC = SHG::getScalarsAbsoluteDifference( G_CNC, exp_G );
        const auto stat_error_G_CNC = SHG::getStatistic( error_G_CNC );
        const auto   error_G_CNC_l2 = SHG::getScalarsNormL2( G_CNC, exp_G );
        trace.info() << "|Ge-G_CNC|_oo = " << stat_error_G_CNC.max() << std::endl;
        trace.info() << "|Ge-G_CNC|_2  = " << error_G_CNC_l2 << std::endl;
        polysurf->addFaceScalarQuantity("CNC H", H_CNC );
    }

    //polysurf->addFaceColorQuantity("Error H He-H", colorErrorH);
    polysurf->addFaceScalarQuantity("Computed H", H );
    polysurf->addFaceScalarQuantity("True H", exp_H );
    polysurf->addFaceScalarQuantity("Error H He-H", error_H );
    polyscope::show();



    return 0;
}
