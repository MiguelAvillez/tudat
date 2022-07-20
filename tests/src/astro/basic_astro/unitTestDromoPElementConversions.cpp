/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Regularised methods for high-efficiency propagation, Geul et al. (2015), AAS/AIAA Astrodynamics Specialist
 *          Conference, Vail, CO, United States. Code archive.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/basic_astro/dromoElementConversions.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"

namespace tudat
{
namespace unit_tests
{

//! Show the functionality of the unit tests.
BOOST_AUTO_TEST_SUITE( test_dromoP_element_conversions )

//! Unit test for conversion of Cartesian to modified equinoctial elements.
BOOST_AUTO_TEST_CASE( testConvertCartesianAndKeplerianElementsToDromoPElements )
{
    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    double tolerance = 1.0E-13;

    Eigen::Vector8d expectedDromoPElements = Eigen::VectorXd::Zero( 8 );
    Eigen::Vector8d computedDromoPElements = Eigen::VectorXd::Zero( 8 );
    Eigen::Vector6d cartesianElements = Eigen::VectorXd::Zero( 6 );

    // Set default Keplerian elements [m,-,rad,rad,rad,rad].
    Eigen::Vector6d keplerianElements = Eigen::VectorXd::Zero( 6 );
    // Initialize perturbing potential [m2/s2]
    double perturbingPotential = 0.0;

    double time = 0.42;
    double initialIndependentVariable = 0.0;
    double currentIndependentVariable = 0.5;

    double gravitationalParameter = 398600.44e9; // Earth's, but any parameter would do.

    for ( int useEnergyElement = 0; useEnergyElement < 1; ++useEnergyElement )
    {
        for ( int timeTypeIndex = 0; timeTypeIndex < 1; ++timeTypeIndex )
        {
            propagators::TimeElementType timeType;
            if ( timeTypeIndex == 0 )
            {
                timeType = propagators::scaled_physical_time;
            }
            else if ( timeTypeIndex == 1 )
            {
                timeType = propagators::linear_time_element;
            }
            else
            {
                timeType = propagators::constant_time_element;
            }

            // Case 1: Elliptical prograde orbit.
            {
                // Default, so no modification necessary
                keplerianElements( semiMajorAxisIndex ) = 1.0e7;
                keplerianElements( eccentricityIndex ) = 0.1;
                keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
                keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
                keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
                keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

                // Overwrite perturbing potential
                perturbingPotential = -5251.8932307137101816;

                // Set expected DromoP elements [-,-,-,-,-,-,-,-].
                // (Results were calculated using Jacco Geul's code)
                if ( timeType == propagators::scaled_physical_time )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0002304240862538;
                }
                else if ( timeType == propagators::linear_time_element )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0327866983245745;
                }
                else // Constant time element
                {
                    expectedDromoPElements( dromoPTimeIndex ) = -0.4015315308083531;
                }
                expectedDromoPElements( dromoPZeta1Index ) = -0.0825372218319105;
                expectedDromoPElements( dromoPZeta2Index ) = -0.0659301474773062;
                if ( useEnergyElement == 0 )
                {
                    expectedDromoPElements( dromoPZeta3Index ) = 1.0533730849570428;
                }
                else
                {
                    expectedDromoPElements( dromoPEnergyIndex ) = -0.5492178393889040;
                }
                expectedDromoPElements( dromoPZeta4Index ) = -0.2228512303422792;
                expectedDromoPElements( dromoPZeta5Index ) = 0.3590870706272545;
                expectedDromoPElements( dromoPZeta6Index ) = -0.8675165654431490;
                expectedDromoPElements( dromoPZeta7Index ) = -0.2623143410585717;

                // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
                cartesianElements = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
                // std::cout << std::setprecision(18) << cartesianElements.transpose() << std::endl;

                double unitOfLength = computeDromoUnitOfLengthFromCartesianElements( cartesianElements );

                // Convert keplerian to Dromo(P) elements
                computedDromoPElements = convertKeplerianToDromoElements(
                        keplerianElements, gravitationalParameter, unitOfLength, perturbingPotential, time,
                        useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, tolerance );

                // Convert cartesian to Dromo(P) elements
                computedDromoPElements = convertCartesianToDromoElements(
                        cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                        time, useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, tolerance );
            }

            // Case 2: Hyperbolic retrograde orbit.
            {
                // Modify Keplerian elements [m,-,rad,rad,rad,rad], i.e. overwrite them.
                keplerianElements( semiMajorAxisIndex ) = -1.0e7;
                keplerianElements( eccentricityIndex ) = 2.0;
                keplerianElements( inclinationIndex ) = convertDegreesToRadians( 170.0 );
                keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
                keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
                keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

                // Overwrite perturbing potential
                perturbingPotential = -8495.2852155969176238;

                // Set expected Dromo(P) elements [-,-,-,-,-,-,-,-].
                // Results were calculated using Jacco Geul's code or by hand (constant time element)
                if ( timeType == propagators::scaled_physical_time )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0002611479545499;
                }
                else if ( timeType == propagators::linear_time_element )
                {
                    // Linear time element is not defined for hyperbolic orbits
                    expectedDromoPElements( dromoPTimeIndex ) = TUDAT_NAN;
                }
                else // Constant time element
                {
                    // Value computed by hand
                    expectedDromoPElements( dromoPTimeIndex ) = -0.0996766923044830;
                }
                expectedDromoPElements( dromoPZeta1Index ) = 1.0995168509937607;
                expectedDromoPElements( dromoPZeta2Index ) = 0.3710210777579507;
                if ( useEnergyElement == 0 )
                {
                    expectedDromoPElements( dromoPZeta3Index ) = 0.5803384950652833;
                }
                else
                {
                    expectedDromoPElements( dromoPEnergyIndex ) = 0.50490058845263460;
                }
                expectedDromoPElements( dromoPZeta4Index ) = -0.9247979668906791;
                expectedDromoPElements( dromoPZeta5Index ) = -0.3703411899059824;
                expectedDromoPElements( dromoPZeta6Index ) = 0.0103557525398912;
                expectedDromoPElements( dromoPZeta7Index ) = -0.0865383260944448;

                // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
                cartesianElements = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
                // std::cout << std::setprecision(18) << cartesianElements.transpose() << std::endl;

                if ( timeType != propagators::linear_time_element )
                {
                    double unitOfLength = computeDromoUnitOfLengthFromCartesianElements( cartesianElements );

                     // Convert keplerian to Dromo(P) elements
                    computedDromoPElements = convertKeplerianToDromoElements(
                            keplerianElements, gravitationalParameter, unitOfLength, perturbingPotential, time,
                            useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                    // Compare computed and expected values
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, tolerance );

                    // Convert cartesian to Dromo(P) elements
                    computedDromoPElements = convertCartesianToDromoElements(
                            cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                            time, useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                    // Compare computed and expected values
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, tolerance );
                }
                // Check if error is thrown for linear time element
                else
                {
                    bool isExceptionFound = false;
                    try
                    {
                        double unitOfLength = computeDromoUnitOfLengthFromCartesianElements( cartesianElements );
                        computedDromoPElements = convertCartesianToDromoElements(
                                cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                                time, useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);

                    }
                    catch( std::runtime_error const& )
                    {
                        isExceptionFound = true;
                    }
                     // Check if runtime error has occurred
                    BOOST_CHECK( isExceptionFound );
                }
            }

            // Case 3: Parabolic retrograde orbit.
            {
                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                keplerianElements( semiLatusRectumIndex ) = 1.0e7;
                keplerianElements( eccentricityIndex ) = 1.0;
                keplerianElements( inclinationIndex ) = convertDegreesToRadians( 170.0 );
                keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
                keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
                keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

                // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
                cartesianElements = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
                // std::cout << std::setprecision(18) << cartesianElements.transpose() << std::endl;

                // Overwrite perturbing potential
                perturbingPotential = -0.0303869376125564;

                // Set expected Dromo(P) elements [-,-,-,-,-,-,-,-].
                // (Results were calculated using Jacco Geul's code)
                if ( timeType == propagators::scaled_physical_time )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0000004965365196;
                }
                else if ( timeType == propagators::linear_time_element )
                {
                    // Linear time element is not defined for parabolic orbits
                    expectedDromoPElements( dromoPTimeIndex ) = TUDAT_NAN;
                }
                else // Constant time element
                {
                    // Constant time element is not defined for parabolic orbits
                    expectedDromoPElements( dromoPTimeIndex ) = TUDAT_NAN;
                }
                expectedDromoPElements( dromoPZeta1Index ) = -6.3363763345734938;
                expectedDromoPElements( dromoPZeta2Index ) = -5.0669337369568375;
                if ( useEnergyElement == 0 )
                {
                    expectedDromoPElements( dromoPZeta3Index ) = 8.1131672390265432;
                }
                else
                {
                    expectedDromoPElements( dromoPEnergyIndex ) = -0.0000000501795938;
                }
                expectedDromoPElements( dromoPZeta4Index ) = -0.5253043567398263;
                expectedDromoPElements( dromoPZeta5Index ) = 0.8464391350216914;
                expectedDromoPElements( dromoPZeta6Index ) = -0.0834253569135866;
                expectedDromoPElements( dromoPZeta7Index ) = -0.0252256480142080;

                double unitOfLength = computeDromoUnitOfLengthFromCartesianElements( cartesianElements );

                // Conversion from keplerian to Dromo(P) elements not defined for parabolic orbits - check if error is
                // thrown
                bool isExceptionFound = false;
                try
                {
                    computedDromoPElements = convertKeplerianToDromoElements(
                        keplerianElements, gravitationalParameter, unitOfLength, perturbingPotential, time,
                        useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                }
                catch( std::runtime_error const& )
                {
                    isExceptionFound = true;
                }
                 // Check if runtime error has occurred
                BOOST_CHECK( isExceptionFound );

                // Convert cartesian to Dromo(P) elements
                // Don't test the linear and constant time elements as those aren't defined for parabolic orbits
                if ( timeType != propagators::linear_time_element and timeType != propagators::constant_time_element )
                {
                    computedDromoPElements = convertCartesianToDromoElements(
                            cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                            time, useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                    // Compare computed and expected values
                    // Because some of the elements are near-zero, a close fraction/percentage check will fail. Therefore,
                    // 1.0 is added to these elements to avoid this.
                    Eigen::Vector8d vectorToAdd = ( Eigen::Vector8d( ) << 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
                    Eigen::Vector8d expectedDromoPElementsPlusOne =  expectedDromoPElements + vectorToAdd;
                    Eigen::Vector8d computedDromoPElementsPlusOne = computedDromoPElements + vectorToAdd;
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElementsPlusOne, computedDromoPElementsPlusOne, tolerance );
                }

            }


            // Case 4: Circular prograde orbit.
            {
                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                keplerianElements( semiMajorAxisIndex ) = 1.0e7;
                keplerianElements( eccentricityIndex ) = 0.0;
                keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
                keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 0.0 );
                keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
                keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

                // Overwrite perturbing potential
                perturbingPotential = -8293.7736712349324080;

                // Set expected DromoP elements [-,-,-,-,-,-,-,-].
                // (Results were calculated using Jacco Geul's code)
                if ( timeType == propagators::scaled_physical_time )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0002651662075303;
                }
                else if ( timeType == propagators::linear_time_element )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0002651662075304;
                }
                else // Constant time element
                {
                    expectedDromoPElements( dromoPTimeIndex ) = -0.4994228875173207;
                }
                expectedDromoPElements( dromoPZeta1Index ) = -0.0003652773723825;
                expectedDromoPElements( dromoPZeta2Index ) = -0.0001995519380162;
                if ( useEnergyElement == 0 )
                {
                    expectedDromoPElements( dromoPZeta3Index ) = 1.0002081373298761;
                }
                else
                {
                    expectedDromoPElements( dromoPEnergyIndex ) = -0.5002080723661827;
                }
                expectedDromoPElements( dromoPZeta4Index ) = -0.1907067137786016;
                expectedDromoPElements( dromoPZeta5Index ) = 0.3771434004148780;
                expectedDromoPElements( dromoPZeta6Index ) = -0.8870776042295483;
                expectedDromoPElements( dromoPZeta7Index ) = -0.1857071051888838;

                // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
                cartesianElements = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
                // std::cout << std::setprecision(18) << cartesianElements.transpose() << std::endl;

                double unitOfLength = computeDromoUnitOfLengthFromCartesianElements( cartesianElements );

                // Convert keplerian to Dromo(P) elements
                computedDromoPElements = convertKeplerianToDromoElements(
                        keplerianElements, gravitationalParameter, unitOfLength, perturbingPotential, time,
                        useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                // Zeta1 and zeta2 require slightly higher tolerance... not sure why
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, 1e-12 );

                // Convert cartesian to Dromo(P) elements
                computedDromoPElements = convertCartesianToDromoElements(
                        cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                        time, useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                // Zeta1 and zeta2 require slightly higher tolerance... not sure why
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, 1e-12 );
            }

            // Case 5: 0 inclination orbit.
            {
                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                keplerianElements( semiMajorAxisIndex ) = 1.0e7;
                keplerianElements( eccentricityIndex ) = 0.1;
                keplerianElements( inclinationIndex ) = convertDegreesToRadians( 0.0 );
                keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 260.0 );
                keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 0.0 );
                keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

                // Overwrite perturbing potential
                perturbingPotential = -6613.9431102345906766;

                // Set expected DromoP elements [-,-,-,-,-,-,-,-].
                // (Results were calculated using Jacco Geul's code)
                if ( timeType == propagators::scaled_physical_time )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0002304240862538;
                }
                else if ( timeType == propagators::linear_time_element )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0327839604446311;
                }
                else // Constant time element
                {
                    expectedDromoPElements( dromoPTimeIndex ) = -0.4014897612054770;
                }
                expectedDromoPElements( dromoPZeta1Index ) = -0.0826104036599911;
                expectedDromoPElements( dromoPZeta2Index ) = -0.0659701268921979;
                if ( useEnergyElement == 0 )
                {
                    expectedDromoPElements( dromoPZeta3Index ) = 1.0534169470365562;
                }
                else
                {
                    expectedDromoPElements( dromoPEnergyIndex ) = -0.5492553639343897;
                }
                expectedDromoPElements( dromoPZeta4Index ) = -0.0000000000000000;
                expectedDromoPElements( dromoPZeta5Index ) = -0.0000000000000000;
                expectedDromoPElements( dromoPZeta6Index ) = -0.3530838749925205;
                expectedDromoPElements( dromoPZeta7Index ) = -0.9355916722696211;

                // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
                cartesianElements = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
                // std::cout << std::setprecision(18) << cartesianElements.transpose() << std::endl;

                double unitOfLength = computeDromoUnitOfLengthFromCartesianElements( cartesianElements );

                // Convert keplerian to Dromo(P) elements
                computedDromoPElements = convertKeplerianToDromoElements(
                        keplerianElements, gravitationalParameter, unitOfLength, perturbingPotential, time,
                        useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, tolerance );

                // Convert cartesian to Dromo(P) elements
                computedDromoPElements = convertCartesianToDromoElements(
                        cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                        time, useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, tolerance );
            }

            // Case 6: 180 inclination orbit.
            {
                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                keplerianElements( semiMajorAxisIndex ) = 1.0e7;
                keplerianElements( eccentricityIndex ) = 0.1;
                keplerianElements( inclinationIndex ) = PI; // = 180 deg
                keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 12.0 );
                keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 0.0 );
                keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 190.0 );

                // Overwrite perturbing potential
                perturbingPotential = -6613.9431102345870386;

                // Set expected DromoP elements [-,-,-,-,-,-,-,-].
                // (Results were calculated using Jacco Geul's code)
                if ( timeType == propagators::scaled_physical_time )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0002304240862538;
                }
                else if ( timeType == propagators::linear_time_element )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = -0.0323231122721236;
                }
                else // Constant time element
                {
                    expectedDromoPElements( dromoPTimeIndex ) = -0.4665968339222320;
                }
                expectedDromoPElements( dromoPZeta1Index ) = -0.1001465392300703;
                expectedDromoPElements( dromoPZeta2Index ) = -0.0338704460448820;
                if ( useEnergyElement == 0 )
                {
                    expectedDromoPElements( dromoPZeta3Index ) = 1.0534169470365560;
                }
                else
                {
                    expectedDromoPElements( dromoPEnergyIndex ) = -0.5492553639343895;
                }
                expectedDromoPElements( dromoPZeta4Index ) = -0.0579812459023547;
                expectedDromoPElements( dromoPZeta5Index ) = 0.9983176724488106;
                expectedDromoPElements( dromoPZeta6Index ) = 0.0000000000000001;
                expectedDromoPElements( dromoPZeta7Index ) = 0.0000000000000000;

                // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
                cartesianElements = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
                // std::cout << std::setprecision(20) << cartesianElements.transpose() << std::endl;

                double unitOfLength = computeDromoUnitOfLengthFromCartesianElements( cartesianElements );

                // Convert keplerian to Dromo(P) elements
                computedDromoPElements = convertKeplerianToDromoElements(
                        keplerianElements, gravitationalParameter, unitOfLength, perturbingPotential, time,
                        useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                // Because some of the elements are near-zero, a close fraction/percentage check will fail. Therefore,
                // 1.0 is added to these elements to avoid this.
                Eigen::Vector8d vectorToAdd = ( Eigen::Vector8d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 ).finished( );
                Eigen::Vector8d expectedDromoPElementsPlusOne =  expectedDromoPElements + vectorToAdd;
                Eigen::Vector8d computedDromoPElementsPlusOne = computedDromoPElements + vectorToAdd;
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElementsPlusOne, computedDromoPElementsPlusOne, tolerance );

                // Convert cartesian to Dromo(P) elements
                computedDromoPElements = convertCartesianToDromoElements(
                        cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                        time, useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                // Because some of the elements are near-zero, a close fraction/percentage check will fail. Therefore,
                // 1.0 is added to these elements to avoid this.
                computedDromoPElementsPlusOne = computedDromoPElements + vectorToAdd;
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElementsPlusOne, computedDromoPElementsPlusOne, tolerance );
            }

            // Case 7: 0 eccentricity and inclination orbit.
            {
                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                keplerianElements( semiMajorAxisIndex ) = 1.0e7;
                keplerianElements( eccentricityIndex ) = 0.0;
                keplerianElements( inclinationIndex ) = 0.0; // = 180 deg
                keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 0.0 );
                keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 0.0 );
                keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

                // Overwrite perturbing potential
                perturbingPotential = -8758.7284433312270266;

                // Set expected DromoP elements [-,-,-,-,-,-,-,-].
                // (Results were calculated using Jacco Geul's code)
                if ( timeType == propagators::scaled_physical_time )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0002651662075303;
                }
                else if ( timeType == propagators::linear_time_element )
                {
                    expectedDromoPElements( dromoPTimeIndex ) = 0.0002651662075303;
                }
                else // Constant time element
                {
                    expectedDromoPElements( dromoPTimeIndex ) = -0.4994054091925029;
                }
                expectedDromoPElements( dromoPZeta1Index ) = -0.0003857595798718;
                expectedDromoPElements( dromoPZeta2Index ) = -0.0002107414189651;
                if ( useEnergyElement == 0 )
                {
                    expectedDromoPElements( dromoPZeta3Index ) = 1.0002198095021353;
                }
                else
                {
                    expectedDromoPElements( dromoPEnergyIndex ) = -0.5002197370490443;
                }
                expectedDromoPElements( dromoPZeta4Index ) = -0.0000000000000000;
                expectedDromoPElements( dromoPZeta5Index ) = -0.0000000000000000;
                expectedDromoPElements( dromoPZeta6Index ) = 0.1620162307246385;
                expectedDromoPElements( dromoPZeta7Index ) = -0.9867880932509171;

                // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
                cartesianElements = convertKeplerianToCartesianElements( keplerianElements, gravitationalParameter );
                // std::cout << std::setprecision(20) << cartesianElements.transpose() << std::endl;

                double unitOfLength = computeDromoUnitOfLengthFromCartesianElements( cartesianElements );

                // Convert keplerian to Dromo(P) elements
                computedDromoPElements = convertKeplerianToDromoElements(
                        keplerianElements, gravitationalParameter, unitOfLength, perturbingPotential, time,
                        useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                // Because some of the elements are near-zero, a close fraction/percentage check will fail. Therefore,
                // 1.0 is added to these elements to avoid this.
                Eigen::Vector8d vectorToAdd = ( Eigen::Vector8d( ) << 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
                Eigen::Vector8d expectedDromoPElementsPlusOne =  expectedDromoPElements + vectorToAdd;
                Eigen::Vector8d computedDromoPElementsPlusOne = computedDromoPElements + vectorToAdd;
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElementsPlusOne, computedDromoPElementsPlusOne, tolerance );

                // Convert cartesian to Dromo(P) elements
                computedDromoPElements = convertCartesianToDromoElements(
                        cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                        time, useEnergyElement, timeType, initialIndependentVariable, currentIndependentVariable);
                // Compare computed and expected values
                // Because some of the elements are near-zero, a close fraction/percentage check will fail. Therefore,
                // 1.0 is added to these elements to avoid this.
                computedDromoPElementsPlusOne = computedDromoPElements + vectorToAdd;
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElementsPlusOne, computedDromoPElementsPlusOne, tolerance );
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat