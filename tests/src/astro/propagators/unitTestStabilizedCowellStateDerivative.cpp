/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <thread>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_stabilized_cowell_propagator )

// Test stabilized Cowell propagator for point mass central body.
BOOST_AUTO_TEST_CASE( testStabilizedCowellPopagatorForPointMassCentralBodies )
{
    // Test simulation for different central body cases
    for( unsigned int simulationCase = 0; simulationCase < 1; simulationCase++ )
    {
        //Using declarations.
        using namespace tudat::interpolators;
        using namespace tudat::numerical_integrators;
        using namespace tudat::spice_interface;
        using namespace tudat::simulation_setup;
        using namespace tudat::basic_astrodynamics;
        using namespace tudat::orbital_element_conversions;
        using namespace tudat::propagators;


        //Load spice kernels.
        spice_interface::loadStandardSpiceKernels( );

        // Define bodies in simulation.
        unsigned int totalNumberOfBodies = 7;
        std::vector< std::string > bodyNames;
        bodyNames.resize( totalNumberOfBodies );
        bodyNames[ 0 ] = "Earth";
        bodyNames[ 1 ] = "Mars";
        bodyNames[ 2 ] = "Sun";
        bodyNames[ 3 ] = "Venus";
        bodyNames[ 4 ] = "Moon";
        bodyNames[ 5 ] = "Mercury";
        bodyNames[ 6 ] = "Jupiter";

        double initialEphemerisTime = 1.0E7;
        double finalEphemerisTime = 2.0E7;
        double maximumTimeStep = 3600.0;
        double buffer = 5.0 * maximumTimeStep;

        // Create bodies needed in simulation
        SystemOfBodies bodies = createSystemOfBodies(
                    getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer ) );

        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
        accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        accelerationsOfEarth[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationMap[ "Earth" ] = accelerationsOfEarth;

//        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
//        accelerationsOfMars[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        accelerationsOfMars[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        accelerationsOfMars[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        accelerationMap[ "Mars" ] = accelerationsOfMars;
//
//        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
//        accelerationsOfMoon[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        accelerationsOfMoon[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//        accelerationMap[ "Moon" ] = accelerationsOfMoon;

        // Propagate Earth, Mars and Moon
        std::vector< std::string > bodiesToPropagate;
        bodiesToPropagate.push_back( "Earth" );
//        bodiesToPropagate.push_back( "Mars" );
//        bodiesToPropagate.push_back( "Moon" );

        unsigned int numberOfNumericalBodies = bodiesToPropagate.size( );

        // Define central bodies: all Sun for simulationCase = 0, Earth and Mars: Sun, Moon: Earth for simulationCase = 1
        std::vector< std::string > centralBodies;
        std::map< std::string, std::string > centralBodyMap;
        centralBodies.resize( numberOfNumericalBodies );
        for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
        {
//            if( i == 2 && simulationCase == 1 )
//            {
//                centralBodies[ i ] = "Earth";
//            }
//            else
//            {
//                centralBodies[ i ] = "Sun";
//            }
            centralBodies[ i ] = "Sun";
            centralBodyMap[ bodiesToPropagate[ i ] ] = centralBodies[ i ];
        }


        // Get initial states for bodies.
        Eigen::VectorXd systemInitialState = Eigen::VectorXd( bodiesToPropagate.size( ) * 6 );
        for( unsigned int i = 0; i < numberOfNumericalBodies ; i++ )
        {
            systemInitialState.segment( i * 6 , 6 ) =
                    bodies.at( bodiesToPropagate[ i ] )->getStateInBaseFrameFromEphemeris( initialEphemerisTime ) -
                    bodies.at( centralBodies[ i ] )->getStateInBaseFrameFromEphemeris( initialEphemerisTime );
        }

        // Create acceleration models.
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, centralBodyMap );

        // Create integrator settings: for cowell and for stabilized cowell
        std::shared_ptr< IntegratorSettings< > > cowellIntegratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4,
                  initialEphemerisTime, 250.0 );

        std::shared_ptr< IntegratorSettings< > > stabilizedCowellIntegratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4,
                  initialEphemerisTime, 4.626e-10 );


        // Select dependent variables to save
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

//        dependentVariables.push_back(
//                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
//                        point_mass_gravity, "Earth", "Jupiter", 1 ) );
        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        point_mass_gravity, "Earth", "Sun", 1 ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        total_acceleration_norm_dependent_variable, "Earth" ) );
//        dependentVariables.push_back(
//                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
//                        point_mass_gravity, "Earth", "Jupiter" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        point_mass_gravity, "Earth", "Sun" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        total_acceleration_dependent_variable, "Earth" ) );


        // Create propagation settings (Cowell)
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime,
                  cowell, std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
//
//        // Propagate orbit with Cowell method
//        SingleArcDynamicsSimulator< double > dynamicsSimulator2(
//                bodies, cowellIntegratorSettings, propagatorSettings, true, false, true );
//
//        // Define ephemeris interrogation settings.
//        double initialTestTime = initialEphemerisTime + 10.0 * maximumTimeStep;
//        double finalTestTime = finalEphemerisTime - 10.0 * maximumTimeStep;
//        double testTimeStep = 1.0E4;
//
//        // Get resutls of Cowell integration at given times.
//        double currentTestTime = initialTestTime;
//        std::map< double, Eigen::Matrix< double, 18, 1 > > cowellIntegrationResults;
//        std::map< double, Eigen::VectorXd > cowellDependentVariables = dynamicsSimulator2.getDependentVariableHistory( );

//        while( currentTestTime < finalTestTime )
//        {
//            cowellIntegrationResults[ currentTestTime ].segment( 0, 6 ) =
//                    bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( currentTestTime );
//            cowellIntegrationResults[ currentTestTime ].segment( 6, 6 ) =
//                    bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( currentTestTime );
//            cowellIntegrationResults[ currentTestTime ].segment( 12, 6 ) =
//                    bodies.at( "Moon" )->getStateInBaseFrameFromEphemeris( currentTestTime );
//
//            currentTestTime += testTimeStep;
//        }

        // Create propagation settings (stabilized cowell)
        propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime,
                  stabilized_cowell, std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

        // Propagate orbit with stabilized cowell method
        SingleArcDynamicsSimulator< double > dynamicsSimulator(
                bodies, stabilizedCowellIntegratorSettings, propagatorSettings, true, false, true );

        // Get results of stabilized cowell integration at given times.
//        currentTestTime = initialTestTime;
        std::map< double, Eigen::Matrix< double, 18, 1 > > stabilizedCowellIntegrationResults;
        std::map< double, Eigen::VectorXd > stabilizedCowellProcessedSolution = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > stabilizedCowellRawSolution = dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );
        std::map< double, Eigen::VectorXd > stabilizedCowellDependentVariables = dynamicsSimulator.getDependentVariableHistory( );

        for (std::map< double, Eigen::VectorXd >::iterator it = stabilizedCowellRawSolution.begin(); it != stabilizedCowellRawSolution.end(); it++)
        {
            std::cout << it->first - 1.0e7 << ':' << it->second.transpose() << std::endl;
        }


//        while( currentTestTime < finalTestTime )
//        {
//            stabilizedCowellIntegrationResults[ currentTestTime ].segment(0, 6 ) =
//                    bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( currentTestTime );
//            stabilizedCowellIntegrationResults[ currentTestTime ].segment(6, 6 ) =
//                    bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( currentTestTime );
//            stabilizedCowellIntegrationResults[ currentTestTime ].segment(12, 6 ) =
//                    bodies.at( "Moon" )->getStateInBaseFrameFromEphemeris( currentTestTime );
//            currentTestTime += testTimeStep;
//        }
//
//        // Compare results of Cowell and stabilized Cowell propagations
//        std::map< double, Eigen::Matrix< double, 18, 1 > >::iterator stabilizedCowellIterator = stabilizedCowellIntegrationResults.begin( );
//        std::map< double, Eigen::Matrix< double, 18, 1 > >::iterator cowellIterator = cowellIntegrationResults.begin( );
//        std::map< double, Eigen::VectorXd >::iterator stabilizedCowellDependentIterator = stabilizedCowellDependentVariables.begin( );
//        std::map< double, Eigen::VectorXd >::iterator cowellDependentIterator = cowellDependentVariables.begin( );
//        for(unsigned int i = 0; i < stabilizedCowellIntegrationResults.size( ); i++ )
//        {
//            for( int j= 0; j< 3; j++ )
//            {
//                BOOST_CHECK_SMALL(( stabilizedCowellIterator->second - cowellIterator->second ).segment(j, 1 )(0 ), 0.01 );
//            }
//
//            for( int j = 6; j < 9; j++ )
//            {
//                BOOST_CHECK_SMALL(( stabilizedCowellIterator->second - cowellIterator->second ).segment(j, 1 )(0 ), 0.01 );
//            }
//
//            for( int j = 12; j < 15; j++ )
//            {
//                BOOST_CHECK_SMALL(( stabilizedCowellIterator->second - cowellIterator->second ).segment(j, 1 )(0 ), 0.1 );
//            }
//
//            for( int j = 3; j < 6; j++ )
//            {
//                BOOST_CHECK_SMALL(( stabilizedCowellIterator->second - cowellIterator->second ).segment(j, 1 )(0 ), 1.0E-8 );
//            }
//
//            for( int j = 9; j < 12; j++ )
//            {
//                BOOST_CHECK_SMALL(( stabilizedCowellIterator->second - cowellIterator->second ).segment(j, 1 )(0 ), 1.0E-8 );
//
//            }
//
//            for( int j = 15; j < 18; j++ )
//            {
//                BOOST_CHECK_SMALL(( stabilizedCowellIterator->second - cowellIterator->second ).segment(j, 1 )(0 ), 1.0E-6 );
//
//            }
//
//            for( int j = 0; j < 12; j++ )
//            {
//                BOOST_CHECK_SMALL(( stabilizedCowellDependentIterator->second - cowellDependentIterator->second )(j ), 1.0E-11 );
//            }
//            stabilizedCowellIterator++;
//            cowellIterator++;
//            stabilizedCowellDependentIterator++;
//            cowellDependentIterator++;
//        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )


}

}

