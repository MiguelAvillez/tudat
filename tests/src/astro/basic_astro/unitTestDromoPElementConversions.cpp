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
#include "tudat/math/basic/basicMathematicsFunctions.h"

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
    keplerianElements( semiMajorAxisIndex ) = 1.0e7;
    keplerianElements( eccentricityIndex ) = 0.1;
    keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );
    double perturbingPotential = -5251.8932307137101816;
    double time = 0.42;
    double initialIndependentVariable = 0.0;
    double currentIndependentVariable = 0.5;

    double gravitationalParameter = 398600.44e9; // Earth's, but any parameter would do.

    for ( int useTimeElement = 0; useTimeElement < 2; ++useTimeElement )
    {
        for ( int timeTypeIndex = 0; timeTypeIndex < 3; ++timeTypeIndex )
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
                if ( useTimeElement == 0 )
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
                cartesianElements = convertKeplerianToCartesianElements( keplerianElements,
                                                                         gravitationalParameter );
                // std::cout << std::setprecision(18) << cartesianElements.transpose() << std::endl;

        //        // Convert to modified equinoctial elements.
        //        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
        //                                                                     gravitationalParameter,
        //                                                                     avoidSingularity );
        //
        //        // Compare, because element 2 is quite small, tolerance is less stringent than usual.
        //        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, 1.0E-13 );

                // Convert cartesian to Dromo(P) elements
                double unitOfLength = computeDromoUnitOfLength( cartesianElements.segment(0, 3) );
                computedDromoPElements = convertCartesianToDromoElements(
                        cartesianElements, gravitationalParameter, unitOfLength, perturbingPotential,
                        time, useTimeElement, timeType, initialIndependentVariable, currentIndependentVariable);

                // std::cout << computedDromoPElements.transpose() << std::endl;

                // Compare computed and expected values
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedDromoPElements, computedDromoPElements, tolerance );
            }
        }
    }



//    // Case 2: Hyperbolic retrograde orbit.
//    {
//        // Set Keplerian elements [m,-,rad,rad,rad,rad].
//        testKepler( semiMajorAxisIndex ) = -1.0e7;
//        testKepler( eccentricityIndex ) = 2.0;
//        testKepler( inclinationIndex ) = convertDegreesToRadians( 170.0 );
//        avoidSingularity = true;
//        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );
//
//        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
//        // hand)
//        testMEE( semiLatusRectumIndex ) = 3.0e7;
//        testMEE( fElementIndex )
//                = 1.8126155740732999264851053135086;
//        testMEE( gElementIndex )
//                = -0.84523652348139887237395697929546;
//        testMEE( hElementIndex )
//                = 0.0845075596072044152327702959491;
//        testMEE( kElementIndex )
//                = 0.02264373235107538825570191377426;
//        testMEE( trueLongitudeIndex )
//                = 6.0213859193804370403867331512857;
//
//        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
//        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
//                                                                     gravitationalParameter );
//
//        // Convert to modified equinoctial elements.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter,
//                                                                     avoidSingularity );
//
//        // Compare.
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, tolerance );
//
//        // Convert to modified equinoctial elements using direct function.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter );
//
//        // Compare.
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, tolerance );
//    }
//
//    // Case 3: Parabolic retrograde orbit.
//    {
//        // Set Keplerian elements [m,-,rad,rad,rad,rad].
//        testKepler( semiMajorAxisIndex ) = 1.0e7;
//        testKepler( eccentricityIndex ) = 1.0;
//        avoidSingularity = true;
//        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );
//
//        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
//        // hand).
//        testMEE( semiLatusRectumIndex ) = 1.0e7;
//        testMEE( fElementIndex )
//                = 0.90630778703664996324255265675432;
//        testMEE( gElementIndex )
//                = -0.42261826174069943618697848964773;
//        testMEE( hElementIndex )
//                = 0.0845075596072044152327702959491;
//        testMEE( kElementIndex )
//                = 0.02264373235107538825570191377426;
//        testMEE( trueLongitudeIndex )
//                = 2.5307274153917778865393516143085;
//
//        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
//        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
//                                                                     gravitationalParameter );
//
//        // Convert to modified equinoctial elements.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter,
//                                                                     avoidSingularity );
//
//        // Compare.
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, tolerance );
//
//        // Convert to modified equinoctial elements using direct function.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter );
//
//        // Compare.
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, tolerance );
//    }
//
//    // Case 4: Circular prograde orbit.
//    {
//        // Set Keplerian elements [m,-,rad,rad,rad,rad].
//        testKepler ( eccentricityIndex ) = 0.0;
//        testKepler ( inclinationIndex ) = convertDegreesToRadians( 50.0 );
//        avoidSingularity = false;
//        testKepler ( argumentOfPeriapsisIndex ) = 0.0; // e = 0, so actually undefined
//
//        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
//        // hand).
//        testMEE ( semiLatusRectumIndex ) = 1.0e7;
//        testMEE ( fElementIndex ) = 0.0;
//        testMEE ( gElementIndex ) = 0.0;
//        testMEE ( hElementIndex )
//                = 0.45041861000828740764931177254188;
//        testMEE ( kElementIndex )
//                = 0.12068930280766941437578622043344;
//        testMEE ( trueLongitudeIndex )
//                = 3.2288591161895097173088279217039;
//
//        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
//        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
//                                                                     gravitationalParameter );
//
//        // Convert to modified equinoctial elements.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter,
//                                                                     avoidSingularity );
//
//        // Check if computed elements match the expected values.
//        // Because two elements are near-zero, a close fraction/percentage check will fail.
//        // Therefore, 1.0 is added to the elements to avoid this
//        Eigen::Vector6d vectorToAdd
//                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
//        Eigen::Vector6d computedMeePlusOne = computedMEE + vectorToAdd;
//        Eigen::Vector6d testMeePlusOne = testMEE + vectorToAdd;
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
//
//        // Convert to modified equinoctial elements using direct function
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter );
//
//        // Check if computed elements match the expected values.
//        // Because two elements are near-zero, a close fraction/percentage check will fail.
//        // Therefore, 1.0 is added to the elements to avoid this.
//        computedMeePlusOne = computedMEE + vectorToAdd;
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
//    }
//
//    // Case 5: 0 inclination orbit.
//    {
//        // Set Keplerian elements [m,-,rad,rad,rad,rad].
//        testKepler ( eccentricityIndex ) = 0.1;
//        testKepler ( inclinationIndex ) = 0.0;
//        avoidSingularity = false;
//        testKepler ( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 260.0 );
//        testKepler ( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 0.0 );
//
//        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
//        // hand).
//        testMEE( semiLatusRectumIndex ) = 9.9e6;
//        testMEE( fElementIndex )
//                = -0.01736481776669303488517166267693;
//        testMEE( gElementIndex )
//                = -0.09848077530122080593667430245895;
//        testMEE( hElementIndex ) = 0.0;
//        testMEE( kElementIndex ) = 0.0;
//        testMEE( trueLongitudeIndex )
//                = 1.221730476396030703846583537942;
//
//        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
//        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
//                                                                     gravitationalParameter );
//
//        // Convert to modified equinoctial elements.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter,
//                                                                     avoidSingularity );
//
//        // Check if computed elements match the expected values.
//        // Because two elements are near-zero, a close fraction/percentage check will fail.
//        // Therefore, 1.0 is added to the elements to avoid this.
//        Eigen::Vector6d vectorToAdd
//                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
//        Eigen::Vector6d computedMeePlusOne = computedMEE + vectorToAdd;
//        Eigen::Vector6d testMeePlusOne = testMEE + vectorToAdd;
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
//
//        // Convert to modified equinoctial elements using direct function.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter );
//
//        // Check if computed elements match the expected values.
//        // Because two elements are near-zero, a close fraction/percentage check will fail.
//        // Therefore, 1.0 is added to the elements to avoid this.
//        computedMeePlusOne = computedMEE + vectorToAdd;
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
//    }
//
//    // Case 6: 180 inclination orbit.
//    {
//        // Set Keplerian elements [m,-,rad,rad,rad,rad].
//        testKepler( eccentricityIndex ) = 0.1;
//        testKepler( inclinationIndex ) = PI; // = 180 deg
//        avoidSingularity = true;
//        testKepler( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 12.0 );
//        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 190.0 );
//
//        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
//        // hand).
//        testMEE ( semiLatusRectumIndex ) = 9.9e6;
//        testMEE ( fElementIndex )
//                = 0.09781476007338056379285667478696;
//        testMEE ( gElementIndex )
//                = 0.02079116908177593371017422844051;
//        testMEE ( hElementIndex )
//                = 0.0;
//        testMEE ( kElementIndex )
//                = 0.0;
//        testMEE ( trueLongitudeIndex )
//                = 3.525565089028545745385855352347;
//
//        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
//        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
//                                                                     gravitationalParameter );
//
//        // Convert to modified equinoctial elements.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter,
//                                                                     avoidSingularity );
//
//        // Check if computed elements match the expected values.
//        // Because two elements are near-zero, a close fraction/percentage check will fail.
//        // Therefore, 1.0 is added to the elements to avoid this.
//        Eigen::Vector6d vectorToAdd
//                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
//        Eigen::Vector6d computedMeePlusOne = computedMEE + vectorToAdd;
//        Eigen::Vector6d testMeePlusOne = testMEE + vectorToAdd;
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
//
//        // Convert to modified equinoctial elements using direct function.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter );
//
//        // Check if computed elements match the expected values.
//        // Because two elements are near-zero, a close fraction/percentage check will fail.
//        // Therefore, 1.0 is added to the elements to avoid this.
//        computedMeePlusOne = computedMEE + vectorToAdd;
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
//    }
//
//    // Case 7: 0 eccentricity and inclination orbit.
//    {
//        // Set Keplerian elements [m,-,rad,rad,rad,rad].
//        testKepler( eccentricityIndex ) = 0.0;
//        testKepler( inclinationIndex ) = 0.0;
//        avoidSingularity = false;
//
//        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
//        // hand).
//        testMEE ( semiLatusRectumIndex ) = 1.0e7; // Circular
//        testMEE ( fElementIndex ) = 0.0;
//        testMEE ( gElementIndex ) = 0.0;
//        testMEE ( hElementIndex ) = 0.0;
//        testMEE ( kElementIndex ) = 0.0;
//        testMEE ( trueLongitudeIndex )
//                = 3.525565089028545745385855352347;
//
//        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
//        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
//                                                                     gravitationalParameter );
//
//        // Convert to modified equinoctial elements.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter,
//                                                                     avoidSingularity );
//
//        // Check if computed elements match the expected values.
//        // Because two elements are near-zero, a close fraction/percentage check will fail.
//        // Therefore, 1.0 is added to the elements to avoid this.
//        Eigen::Vector6d vectorToAdd
//                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
//        Eigen::Vector6d computedMeePlusOne = computedMEE + vectorToAdd;
//        Eigen::Vector6d testMeePlusOne = testMEE + vectorToAdd;
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
//
//        // Convert to modified equinoctial elements using direct function.
//        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
//                                                                     gravitationalParameter );
//
//        // Check if computed elements match the expected values.
//        // Because two elements are near-zero, a close fraction/percentage check will fail.
//        // Therefore, 1.0 is added to the elements to avoid this.
//        computedMeePlusOne = computedMEE + vectorToAdd;
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
//    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat