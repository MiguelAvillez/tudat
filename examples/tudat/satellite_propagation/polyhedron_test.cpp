/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/simulation.h>

#include "tudat/io/applicationOutput.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::input_output;
using namespace tudat::reference_frames;
using namespace tudat::gravitation;

using namespace tudat::orbital_element_conversions;
using namespace tudat::numerical_integrators;
using namespace tudat::propagators;
using namespace tudat::ephemerides;

//! Execute propagation of orbit of Asterix around the Earth.
int main( )
{
    const std::string associatedReferenceFrame = "IAU_Sun";

    Eigen::MatrixXd verticesCoordinates(8,3);
    Eigen::MatrixXi verticesDefiningEachFacet(12,3);
    verticesCoordinates <<
        0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00,
        20.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00,
        0.000000000000000000e+00, 10.000000000000000000e+00, 0.000000000000000000e+00,
        20.000000000000000000e+00, 10.000000000000000000e+00, 0.000000000000000000e+00,
        0.000000000000000000e+00, 0.000000000000000000e+00, 10.000000000000000000e+00,
        20.000000000000000000e+00, 0.000000000000000000e+00, 10.000000000000000000e+00,
        0.000000000000000000e+00, 10.000000000000000000e+00, 10.000000000000000000e+00,
        20.000000000000000000e+00, 10.000000000000000000e+00, 10.000000000000000000e+00;
    verticesCoordinates = verticesCoordinates * 1e3;
    verticesDefiningEachFacet <<
        2, 1, 0,
        1, 2, 3,
        4, 2, 0,
        2, 4, 6,
        1, 4, 0,
        4, 1, 5,
        6, 5, 7,
        5, 6, 4,
        3, 6, 7,
        6, 3, 2,
        5, 3, 7,
        3, 5, 1;

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 150.0; // tudat::physical_constants::JULIAN_DAY;
    const double environmentTimeBuffer = 0.0;

    // Create bodies in simulation
    std::vector< std::string > bodiesToCreate = { "Sun", "Earth" };
    // Load Spice kernels.
//    spice_interface::loadStandardSpiceKernels( );
//    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - environmentTimeBuffer,
//                                    simulationEndEpoch + environmentTimeBuffer,  "Earth", "J2000" );
//    bodySettings.resetFrames("SSB", "J2000");
    BodyListSettings bodySettings = BodyListSettings();

    bodySettings.addSettings("Sun");
    bodySettings.get("Sun")->gravityFieldSettings = //centralGravitySettings(1.327e20);
            polyhedronGravitySettings(
            1.327e20,
            verticesCoordinates, verticesDefiningEachFacet, associatedReferenceFrame);
    Eigen::Vector6d sunState;
    sunState << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    bodySettings.get("Sun")->ephemerisSettings = constantEphemerisSettings(
            sunState, bodySettings.getFrameOrigin(), bodySettings.getFrameOrientation());
    Eigen::Matrix3d initialOrientation;
    initialOrientation.setIdentity();
    bodySettings.get("Sun")->rotationModelSettings = simpleRotationModelSettings(
            bodySettings.getFrameOrientation(), associatedReferenceFrame, initialOrientation, simulationStartEpoch, 1e-5);

    bodySettings.addSettings("Earth");
    bodySettings.get("Earth")->gravityFieldSettings = centralGravitySettings(3.986e14);
    Eigen::Vector6d earthState;
    earthState << 1.5e11, 0.0, 0.0, 0.0, 0.0, 0.0;
    bodySettings.get("Earth")->ephemerisSettings = keplerEphemerisSettings(
            earthState, simulationStartEpoch, 1.327e20, bodySettings.getFrameOrigin(), bodySettings.getFrameOrientation());

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Asterix" );
    bodies.at( "Asterix" )->setConstantBodyMass( 0.0 );

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ] = { pointMassGravityAcceleration( ) };
    accelerationsOfAsterix[ "Sun" ] = { polyhedronAcceleration( ) }; // { polyhedronAcceleration( ) }; //{ pointMassGravityAcceleration( ) };

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationSettingsList;
    accelerationSettingsList[ "Asterix" ] = accelerationsOfAsterix;

    std::vector< std::string > bodiesToPropagate = { "Asterix" };
    std::vector< std::string > centralBodies = { "Earth" };

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationSettingsList, bodiesToPropagate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = getBodyGravitationalParameter( bodies, "Earth" );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch );

    const double fixedStepSize = 30.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > > ( rungeKutta4, 0.0, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    dynamicsSimulator.integrateEquationsOfMotion(asterixInitialState);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd finalIntegratedState = (--integrationResult.end( ) )->second;
    // Print the position (in km) and the velocity (in km/s) at t = 0.
    std::cout << "The initial position vector of Asterix is [km]:" << std::endl <<
                 asterixInitialState.segment( 0, 3 ).transpose() / 1E3 << std::endl;
    // Print the position (in km) and the velocity (in km/s) at t = 86400.
    std::cout << "After " << simulationEndEpoch <<
                 " seconds, the position vector of Asterix is [km]:" << std::endl <<
                 finalIntegratedState.segment( 0, 3 ).transpose() / 1E3 << std::endl;

    std::cerr << "All fine! ###########################################" << std::endl;
}
