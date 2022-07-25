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
 *      A new set of integrals of motion to propagate the perturbed two-body problem, Bau et al. (2013),
 *          Celestial Mechanics and Dynamics Astronomy, 116, pp. 3-78
 *      Regularised methods for high-efficiency propagation, Geul et al. (2015), AAS/AIAA Astrodynamics Specialist
 *          Conference, Vail, CO, United States. Code archive.
 *      TIME ELEMENTS FOR ENHANCED PERFORMANCE OF THE DROMO ORBIT PROPAGATOR, Bau and Bombardelli (2014),
 *          The Astronomical Journal, 148:43
 *
 */

#ifndef TUDAT_DROMO_ELEMENT_CONVERSIONS_H
#define TUDAT_DROMO_ELEMENT_CONVERSIONS_H

#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Computes the unit of length used in the Dromo propagator.
/*!
 * Computes the unit of length used to make variables dimensionless in the Dromo propagator, starting from the cartesian
 * elements. The unit of length is selected to be the inital radial position. Bau (2013), section 2.
 *
 * @param initialCartesianElements Initial cartesian state.
 * @return Unit of length.
 */
double computeDromoUnitOfLengthFromCartesianElements (const Eigen::Vector6d& initialCartesianElements )
{
    return initialCartesianElements.segment( 0, 3 ).norm();
}

//! Computes the unit of time used in the Dromo(P) propagator.
/*!
 * Computes the unit of time used to make variables dimensionless in the Dromo(P) propagator. The unit of time is
 * selected to be angular rate of a circular orbit with the unit of length as radius. Bau (2013), section 2.
 *
 * @param centralBodyGravitationalParameter Central body gravitational parameter.
 * @return Unit of time.
 */
double computeDromoUnitOfTime (const double centralBodyGravitationalParameter,
                               const double unitOfLength)
{
    return 1 / std::sqrt( centralBodyGravitationalParameter / std::pow(unitOfLength, 3) );
}

//! Computes the MPhi matrix used in conversions for the Dromo(P) propagator.
/*!
 * Computes the MPhi matrix used when converting the Dromo(P) elements to cartesian and vice-versa. Eq. 50 of Bau (2013).
 *
 * @param initialIndependentVariable Initial independent variable.
 * @param currentIndependentVariable Current value of the independent variable.
 * @return MPhi matrix
 */
Eigen::Matrix3d computeDromoMPhiMatrix (const double initialIndependentVariable,
                                        const double currentIndependentVariable )
{
    double dPhi = currentIndependentVariable - initialIndependentVariable;
    return (Eigen::Matrix3d() <<
        std::cos( dPhi ), -std::sin( dPhi ), 0,
        std::sin( dPhi ), std::cos( dPhi ), 0,
        0, 0, 1).finished();
}

//! Convert the time to the dimensionless time used in Dromo(P).
/*!
 * Convert the time to the dimensionless time used in Dromo(P).
 *
 * @param physicalTime Dimensional physical time.
 * @param centralBodyGravitationalParameter Central body gravitational parameter.
 * @param unitOfLength Dromo(P) unit of length.
 * @return Dromo(P) dimensionless time.
 */
double convertPhysicalTimeToScaledTime(const double physicalTime,
                                       const double centralBodyGravitationalParameter,
                                       const double unitOfLength)
{
    return physicalTime / computeDromoUnitOfTime( centralBodyGravitationalParameter, unitOfLength);
}

//! Convert the dimensionless time used in Dromo(P) to the physical time.
/*!
 * Convert the dimensionless time used in Dromo(P) to the physical time.
 *
 * @param scaledTime Dromo(P) dimensionless time.
 * @param centralBodyGravitationalParameter Central body gravitational parameter.
 * @param unitOfLength Dromo(P) unit of length.
 * @return Dimensional physical time.
 */
double convertScaledTimeToPhysicalTime(const double scaledTime,
                                       const double centralBodyGravitationalParameter,
                                       const double unitOfLength)
{
    return scaledTime * computeDromoUnitOfTime( centralBodyGravitationalParameter, unitOfLength );
}

//! Computes the energy and zeta3 elements for Dromo(P).
/*!
 * Computes the energy and zeta3 elements for Dromo(P). If the dromo elements use the energy, then zeta3 is computed;
 * if the dromo elements use zeta3, then the energy is computed. Both are returned by reference
 *
 * @param dromoElementsExceptTime Dromo(P) state such time element doesn't have to be set.
 * @param usingEnergyElement Flag indicating whether the dromo state uses energy (true) or zeta3 (false).
 * @param energy Energy element.
 * @param zeta3 Zeta3 element.
 */
void computeEnergyAndDromoZeta3(const Eigen::Vector8d& dromoElementsExceptTime,
                                const bool usingEnergyElement,
                                double& energy,
                                double& zeta3)
{
    if ( usingEnergyElement )
    {
        energy = dromoElementsExceptTime( dromoPEnergyIndex );
        zeta3 = std::sqrt( std::pow( dromoElementsExceptTime( dromoPZeta1Index), 2 ) +
                           std::pow( dromoElementsExceptTime( dromoPZeta2Index), 2 ) - 2 * energy );
    }
    else
    {
        zeta3 = dromoElementsExceptTime ( dromoPZeta3Index );
        energy = ( std::pow( dromoElementsExceptTime( dromoPZeta1Index), 2 ) +
                   std::pow( dromoElementsExceptTime( dromoPZeta2Index), 2 ) -
                   std::pow( dromoElementsExceptTime( dromoPZeta3Index), 2 ) ) / 2;
    }
}

//! Compute the Dromo(P) term s.
/*!
 * Compute the Dromo(P) term s, Eq. 40 of Bau (2013).
 *
 * @param dromoElements Dromo(P) state.
 * @param currentIndependentVariable Current independent variable.
 * @param zeta3 Value of the zeta3 Dromo(P) element.
 * @return s
 */
double computeDromoS(const Eigen::Vector8d& dromoElements,
                     const double currentIndependentVariable,
                     const double zeta3)
{
    return zeta3 + dromoElements( dromoPZeta1Index ) * std::cos( currentIndependentVariable) +
           dromoElements( dromoPZeta2Index ) * std::sin( currentIndependentVariable );
}

//! Compute the radial velocity.
/*!
 * Compute the radial velocity, Eq. 41 of Bau (2013).
 *
 * @param dromoElements Dromo(P) state.
 * @param currentIndependentVariable Current independent variable.
 * @return Radial velocity
 */
double computeDromoRadialVelocity (const Eigen::Vector8d& dromoElements,
                                   const double currentIndependentVariable)
{
    return dromoElements( dromoPZeta1Index ) * std::sin( currentIndependentVariable) -
           dromoElements( dromoPZeta2Index ) * std::cos( currentIndependentVariable );
}

//! Compute the transverse velocity.
/*!
 * Compute the transverse velocity., Eq. 42 of Bau (2013).
 *
 * @param s Dromo(P) s term.
 * @param perturbingPotential Potential of perturbing accelerations.
 * @return Transverse velocity
 */
double computeDromoTransverseVelocity (const double s,
                                       const double perturbingPotential)
{
    return std::sqrt( std::pow( s, 2 ) - 2 * perturbingPotential );
}

//! Computes the difference between the linear time element and dimensionless time.
/*!
 * Computes the difference between the linear time element and dimensionless time, Eq. 6 of Geul (2016).
 *
 * @param dromoElementsExceptTime Dromo(P) state such time element doesn't have to be set.
 * @param currentIndependentVariable Current value of the independent variable.
 * @param usingEnergyElement Flag indicating whether the dromo state uses energy (true) or zeta3 (false).
 * @return Difference.
 */
double computeDromoScaledTimeToLinearTimeElementBias(const Eigen::Vector8d& dromoElementsExceptTime,
                                                     const double currentIndependentVariable,
                                                     const bool usingEnergyElement)
{
    double energy;
    double zeta3;
    computeEnergyAndDromoZeta3( dromoElementsExceptTime, usingEnergyElement, energy, zeta3 );

    if ( energy >= 0 )
    {
        throw std::runtime_error( "Error in Dromo(P) when converting between time and linear time element, the latter "
                                  "is not defined for parabolic and hyperbolic orbits." );
    }

    double s = computeDromoS( dromoElementsExceptTime, currentIndependentVariable, zeta3 );
    double radialVelocity = computeDromoRadialVelocity( dromoElementsExceptTime, currentIndependentVariable );

    return - radialVelocity / ( 2 * energy * zeta3 * s ) - 1 / ( energy * std::sqrt( -2 * energy) ) *
                                                           std::atan2( radialVelocity, s + std::sqrt( -2 * energy) );
}

//! Converts the dimensional physical time to linear time element.
/*!
 * Converts the dimensional physical time to linear time element, Eq. 6 of Geul (2016).
 *
 * @param physicalTime Dimensional physical time.
 * @param currentIndependentVariable Current value of the independent variable.
 * @param dromoElementsExceptTime Dromo(P) state such time element doesn't have to be set.
 * @param centralBodyGravitationalParameter Central body gravitational parameter.
 * @param unitOfLength Dromo(P) unit of length.
 * @param usingEnergyElement Flag indicating whether the dromo state uses energy (true) or zeta3 (false).
 * @return Linear time elemnent.
 */
double convertPhysicalTimeToDromoLinearTimeElement(const double physicalTime,
                                                   const double currentIndependentVariable,
                                                   const Eigen::Vector8d& dromoElementsExceptTime,
                                                   const double centralBodyGravitationalParameter,
                                                   const double unitOfLength,
                                                   const bool usingEnergyElement)
{
    return convertPhysicalTimeToScaledTime( physicalTime, centralBodyGravitationalParameter, unitOfLength ) +
           computeDromoScaledTimeToLinearTimeElementBias( dromoElementsExceptTime, currentIndependentVariable,
                                                           usingEnergyElement );
}

//! Converts linear time element to the dimensional physical time.
/*!
 * Converts linear time element to the dimensional physical time, Eq. 6 of Geul (2016).
 *
 * @param linearTimeElement Linear time elemnent.
 * @param currentIndependentVariable Current value of the independent variable.
 * @param dromoElementsExceptTime Dromo(P) state such time element doesn't have to be set.
 * @param centralBodyGravitationalParameter Central body gravitational parameter.
 * @param unitOfLength Dromo(P) unit of length.
 * @param usingEnergyElement Flag indicating whether the dromo state uses energy (true) or zeta3 (false).
 * @return Dimensional physical time.
 */
double convertDromoLinearTimeElementToPhysicalTime(const double linearTimeElement,
                                                   const double currentIndependentVariable,
                                                   const Eigen::Vector8d& dromoElementsExceptTime,
                                                   const double centralBodyGravitationalParameter,
                                                   const double unitOfLength,
                                                   const bool usingEnergyElement)
{
    double energy;
    double zeta3;
    computeEnergyAndDromoZeta3( dromoElementsExceptTime, usingEnergyElement, energy, zeta3 );

    if ( energy >= 0 )
    {
        throw std::runtime_error( "Error in Dromo(P) when converting between linear time element and time, the former "
                                  "is not defined for parabolic and hyperbolic orbits." );
    }

    const double scaledTime = linearTimeElement - computeDromoScaledTimeToLinearTimeElementBias(
            dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement );
    return convertScaledTimeToPhysicalTime( scaledTime, centralBodyGravitationalParameter, unitOfLength );
}

//! Computes the difference between the constant time element and linear time element.
/*!
 * Computes the difference between the constant time element and linear time element, Eq. 7 of Geul (2016).
 *
 * @param dromoElementsExceptTime Dromo(P) state such time element doesn't have to be set.
 * @param currentIndependentVariable Current value of the independent variable.
 * @param usingEnergyElement Flag indicating whether the dromo state uses energy (true) or zeta3 (false).
 * @return Difference.
 */
double computeDromoLinearToConstantTimeElementBias(const Eigen::Vector8d& dromoElementsExceptTime,
                                                   const double currentIndependentVariable,
                                                   const bool usingEnergyElement)
{
    double energy;
    double zeta3;
    computeEnergyAndDromoZeta3( dromoElementsExceptTime, usingEnergyElement, energy, zeta3 );

    double semiMajorAxis = -1 / ( 2 * energy );
    return - std::pow( semiMajorAxis, 3.0/2.0 ) * currentIndependentVariable;
}

//! Computes the difference between the constant time element and dimensionless time.
/*!
 * Computes the difference between the constant time element and dimensionless time, Eq. B1 of Bau (2014).
 *
 * @param dromoElementsExceptTime Dromo(P) state such time element doesn't have to be set.
 * @param currentIndependentVariable Current value of the independent variable.
 * @param usingEnergyElement Flag indicating whether the dromo state uses energy (true) or zeta3 (false).
 * @return Difference.
 */
double computeDromoScaledTimeToHyperbolicConstantTimeElementBias(const Eigen::Vector8d& dromoElementsExceptTime,
                                                                 const double currentIndependentVariable,
                                                                 const bool usingEnergyElement)
{
    double energy;
    double zeta3;
    computeEnergyAndDromoZeta3( dromoElementsExceptTime, usingEnergyElement, energy, zeta3 );

    double s = computeDromoS( dromoElementsExceptTime, currentIndependentVariable, zeta3 );
    double radialVelocity = computeDromoRadialVelocity( dromoElementsExceptTime, currentIndependentVariable );

    double positiveSemiMajorAxis = 1 / ( 2 * energy );

    double generalizedEccentricity = std::sqrt(1 + 2.0 * energy / std::pow(zeta3, 2) );

    double k3 =
        dromoElementsExceptTime( orbital_element_conversions::dromoPZeta1Index ) * std::cos( currentIndependentVariable ) +
        dromoElementsExceptTime( orbital_element_conversions::dromoPZeta2Index ) * std::sin( currentIndependentVariable );

    return - radialVelocity * positiveSemiMajorAxis / ( zeta3 * s ) + 2 * std::pow( positiveSemiMajorAxis, 3.0 / 2.0 ) *
        std::atanh( radialVelocity / ( zeta3 * generalizedEccentricity + k3 ) *
        std::sqrt( ( generalizedEccentricity - 1 ) / ( generalizedEccentricity + 1 ) ) );
}

//! Converts the dimensional physical time to constant time element.
/*!
 * Converts the dimensional physical time to linear time element, Eq. 7 of Geul (2016).
 *
 * @param physicalTime Dimensional physical time.
 * @param currentIndependentVariable Current value of the independent variable.
 * @param dromoElementsExceptTime Dromo(P) state such time element doesn't have to be set.
 * @param centralBodyGravitationalParameter Central body gravitational parameter.
 * @param unitOfLength Dromo(P) unit of length.
 * @param usingEnergyElement Flag indicating whether the dromo state uses energy (true) or zeta3 (false).
 * @return Constant time elemnent.
 */
double convertPhysicalTimeToDromoConstantTimeElement(const double physicalTime,
                                                     const double currentIndependentVariable,
                                                     const Eigen::Vector8d& dromoElementsExceptTime,
                                                     const double centralBodyGravitationalParameter,
                                                     const double unitOfLength,
                                                     const bool usingEnergyElement)
{
    double energy;
    double zeta3;
    computeEnergyAndDromoZeta3( dromoElementsExceptTime, usingEnergyElement, energy, zeta3 );

    double constantTimeElement;

    const double tolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    if ( std::abs(energy) < tolerance ) // Parabolic orbit
    {
        throw std::runtime_error( "Error in Dromo(P) when converting between time and constant time element, the latter "
                                  "is not defined for parabolic orbits." );
    }
    else if ( energy < - tolerance ) // Eccentric orbit
    {
        double linearTimeElement = convertPhysicalTimeToDromoLinearTimeElement(
            physicalTime, currentIndependentVariable, dromoElementsExceptTime, centralBodyGravitationalParameter, unitOfLength,
            usingEnergyElement );

        constantTimeElement = linearTimeElement + computeDromoLinearToConstantTimeElementBias(
            dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement);
    }
    else
    {
        constantTimeElement = convertPhysicalTimeToScaledTime( physicalTime, centralBodyGravitationalParameter, unitOfLength ) +
                computeDromoScaledTimeToHyperbolicConstantTimeElementBias( dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement );
    }

    return constantTimeElement;
}

//! Converts constant time element to the dimensional physical time.
/*!
 * Converts the constant time element to dimensional physical time, Eq. 7 of Geul (2016).
 *
 * @param constantTimeElement Linear time elemnent.
 * @param currentIndependentVariable Current value of the independent variable.
 * @param dromoElementsExceptTime Dromo(P) state such time element doesn't have to be set.
 * @param centralBodyGravitationalParameter Central body gravitational parameter.
 * @param unitOfLength Dromo(P) unit of length.
 * @param usingEnergyElement Flag indicating whether the dromo state uses energy (true) or zeta3 (false).
 * @return Dimensional physical time.
 */
double convertDromoConstantTimeElementToPhysicalTime(const double constantTimeElement,
                                                     const double currentIndependentVariable,
                                                     const Eigen::Vector8d& dromoElementsExceptTime,
                                                     const double centralBodyGravitationalParameter,
                                                     const double unitOfLength,
                                                     const bool usingEnergyElement)
{
    double time;

    double energy;
    double zeta3;
    computeEnergyAndDromoZeta3( dromoElementsExceptTime, usingEnergyElement, energy, zeta3 );

    const double tolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    if ( std::abs(energy) < tolerance ) // Parabolic orbit
    {
        throw std::runtime_error( "Error in Dromo(P) when converting between constant time element and time, the former "
                                  "is not defined for parabolic orbits." );
    }
    else if ( energy < - tolerance ) // Eccentric orbit
    {
        double linearTimeElement = constantTimeElement - computeDromoLinearToConstantTimeElementBias(
            dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement);

        time = convertDromoLinearTimeElementToPhysicalTime(
            linearTimeElement, currentIndependentVariable, dromoElementsExceptTime, centralBodyGravitationalParameter,
            unitOfLength, usingEnergyElement );
    }
    else
    {
        double scaledTime = constantTimeElement - computeDromoScaledTimeToHyperbolicConstantTimeElementBias(
                dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement );

        time = convertScaledTimeToPhysicalTime( scaledTime, centralBodyGravitationalParameter, unitOfLength );
    }

    return time;
}

//! Converts Cartesian to Dromo(P) elements, following Bau (2013) and Bau (2014).
/*!
 * Converts Cartesian to Dromo(P) elements, following Bau (2013) and Bau (2014). Zeta7 is always taken to be negative.
 * @param cartesianElements Converted state in Cartesian elements. The order of elements is fixed!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 * @param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * @param unitOfLength Unit of length [m] used to make variables dimensionless. Usually it is taken to be the the position
 *          vector norm at the initial state.
 * @param perturbingPotential Potential of perturbing accelerations [m^2/s^2].
 * @param timeFromPropagationStart Time since the start of the propagation [s].
 * @param useEnergyElement Flag indicating whether the dromo state should use energy (true) or zeta3 (false).
 * @param timeType Enum indicating the type of time element which should be used.
 * @param initialIndependentVariable Initial value of the independent variable. Default is TUDAT_NAN, in which case
 *          the initial independent variable is taken to be the true anomaly.
 * @param currentIndependentVariable Current value of the independent variable. Default is TUDAT_NAN, in which case
 *          the current independent variable is taken to be the same as the initial independent variable.
 * @return dromoElements State in Dromo(P) elements.
 *         The order of elements is fixed
 *         dromoElements( 0 ) = Time element,                 [-]
 *         dromoElements( 1 ) = Zeta1 element,                [-]
 *         dromoElements( 2 ) = Zeta2 element,                [-]
 *         dromoElements( 3 ) = Energy or zeta3 element,      [-]
 *         dromoElements( 4 ) = Zeta4 element,                [-]
 *         dromoElements( 5 ) = Zeta5 element,                [-]
 *         dromoElements( 6 ) = Zeta6 element.                [-]
 *         dromoElements( 7 ) = Zeta7 element.                [-]
 */
Eigen::Vector8d convertCartesianToDromoElements(const Eigen::Vector6d& cartesianElements,
                                                const double centralBodyGravitationalParameter,
                                                const double unitOfLength,
                                                const double perturbingPotential,
                                                const double timeFromPropagationStart,
                                                const bool useEnergyElement,
                                                const propagators::TimeElementType timeType,
                                                double initialIndependentVariable = TUDAT_NAN,
                                                double currentIndependentVariable = TUDAT_NAN)
{
    // Declaring eventual output vector.
    Eigen::Vector8d dromoElements;

    // Check if provided unit of length is valid
    if ( unitOfLength <= 0.0 )
    {
        throw std::runtime_error( "Error when converting Cartesian to Dromo(P) elements, unit of length (" +
                std::to_string( unitOfLength ) + ") is <= zero." );
    }

    // If current independent variable isn't defined set values to default ones
    if ( isnan(currentIndependentVariable) )
    {
        if ( isnan(initialIndependentVariable) )
        {
            Eigen::Vector6d keplerElements = convertCartesianToKeplerianElements( cartesianElements,
                                                                                  centralBodyGravitationalParameter );
            initialIndependentVariable = keplerElements( trueAnomalyIndex );
        }
        currentIndependentVariable = initialIndependentVariable;
    }
    else if ( currentIndependentVariable < initialIndependentVariable )
    {
        throw std::runtime_error( "Error when converting Cartesian to Dromo(P) elements, initial independent variable (" +
                std::to_string( initialIndependentVariable ) + ") is larger than current independent variable (" +
                std::to_string( currentIndependentVariable ) + ")." );
    }

    // Compute dimensionless cartesian position and velocity
    double unitOfTime = computeDromoUnitOfTime(centralBodyGravitationalParameter, unitOfLength);
    Eigen::Vector3d cartesianPosition = cartesianElements.segment(0, 3 ) / unitOfLength;
    double cartesianPositionNorm = cartesianPosition.norm();
    Eigen::Vector3d cartesianVelocity = cartesianElements.segment(3, 3 ) / unitOfLength * unitOfTime;
    // Compute dimensionless potential
    double dimensionlessPerturbingPotential = perturbingPotential * std::pow( unitOfTime, 2 ) / std::pow( unitOfLength, 2 );

    Eigen::Vector3d angularMomentum = cartesianPosition.cross( cartesianVelocity );
    double angularMomentumNorm = angularMomentum.norm();

    double radialVelocityNorm = cartesianPosition.dot( cartesianVelocity ) / cartesianPositionNorm;
    double transverseVelocityNorm = angularMomentumNorm / cartesianPositionNorm;

    if ( perturbingPotential > - std::pow( transverseVelocityNorm, 2 ) / 2 )
    {
        throw std::runtime_error(
            "Value of perturbing potential (" + std::to_string( perturbingPotential ) + ") is too large (maximum value is " +
            std::to_string( - std::pow( transverseVelocityNorm, 2 ) / 2 ) + "), which will lead to singularities. For "
            "details see Section 6 of A New Set of Integrals of Motion to Propagate the Perturbed Two-Body Problem, Bau et al. (2013)." );
    }

    // Auxiliary variables
    double s = std::sqrt( std::pow( transverseVelocityNorm, 2 ) + 2 * dimensionlessPerturbingPotential );
    double zeta3 = 1 / s;

    Eigen::Matrix3d QRI;
    QRI << cartesianPosition / cartesianPositionNorm,
        angularMomentum.cross( cartesianPosition ) / angularMomentumNorm / cartesianPositionNorm,
        angularMomentum / angularMomentumNorm;
    Eigen::Matrix3d Q0 = QRI * computeDromoMPhiMatrix( initialIndependentVariable, currentIndependentVariable ).transpose();

    // Compute dromo elements
    dromoElements( dromoPZeta1Index ) = ( s - zeta3 ) * std::cos( currentIndependentVariable ) +
            radialVelocityNorm * std::sin( currentIndependentVariable );
    dromoElements( dromoPZeta2Index ) = ( s - zeta3 ) * std::sin( currentIndependentVariable ) -
            radialVelocityNorm * std::cos( currentIndependentVariable );

    int sign = -1; // Can either take the value -1 or 1
    dromoElements( dromoPZeta7Index ) = sign * 0.5 * std::sqrt( 1 + Q0( 0, 0) + Q0( 1, 1) + Q0( 2, 2) );

    // Define the tolerance of a singularity
    const double tolerance = 20.0 * std::numeric_limits< double >::epsilon( );
    if ( std::abs( dromoElements( dromoPZeta7Index ) ) > tolerance )
    {
        dromoElements( dromoPZeta4Index ) = ( Q0( 2, 1) - Q0( 1, 2) ) / ( 4 * dromoElements( dromoPZeta7Index ) );
        dromoElements( dromoPZeta5Index ) = ( Q0( 0, 2) - Q0( 2, 0) ) / ( 4 * dromoElements( dromoPZeta7Index ) );
        dromoElements( dromoPZeta6Index ) = ( Q0( 1, 0) - Q0( 0, 1) ) / ( 4 * dromoElements( dromoPZeta7Index ) );
    }
    else
    {
        dromoElements( dromoPZeta6Index ) = sign * std::sqrt( ( Q0( 2, 2) + 1 ) / 2 );
        if ( std::abs( dromoElements( dromoPZeta6Index ) ) > tolerance )
        {
            dromoElements( dromoPZeta4Index ) = Q0( 0, 2) / ( 2 * dromoElements( dromoPZeta6Index ) );
            dromoElements( dromoPZeta5Index ) = Q0( 1, 2) / ( 2 * dromoElements( dromoPZeta6Index ) );
        }
        else
        {
            dromoElements( dromoPZeta4Index ) = sign * std::sqrt( ( 1 - Q0( 1, 1) ) / 2 );
            if ( std::abs( dromoElements( dromoPZeta4Index ) ) > tolerance )
            {
                dromoElements( dromoPZeta5Index ) = Q0( 0, 1) / ( 2 * dromoElements( dromoPZeta4Index ) );
            }
            else
            {
                dromoElements( dromoPZeta5Index ) = sign * 1;
            }
        }
    }

    // Select the dromo element 3 to use based on the flag
    if ( useEnergyElement )
    {
        // Computing energy: the second term is -1/r instead of -mu/r because the position and velocity are dimensionless
        dromoElements( dromoPZeta3Index ) = std::pow( cartesianVelocity.norm(), 2 ) / 2 -
                1 / cartesianPositionNorm + dimensionlessPerturbingPotential;
    }
    else
    {
        dromoElements( dromoPZeta3Index ) = zeta3;
    }

    // Compute the time element
    if ( timeType == propagators::scaled_physical_time )
    {
        dromoElements( dromoPTimeIndex ) = convertPhysicalTimeToScaledTime( timeFromPropagationStart,
                                                                            centralBodyGravitationalParameter,
                                                                            unitOfLength );
    }
    else if ( timeType == propagators::linear_time_element )
    {
        dromoElements( dromoPTimeIndex ) = convertPhysicalTimeToDromoLinearTimeElement(
                timeFromPropagationStart, currentIndependentVariable, dromoElements, centralBodyGravitationalParameter,
                unitOfLength, useEnergyElement );
    }
    else if ( timeType == propagators::constant_time_element )
    {
        dromoElements( dromoPTimeIndex ) = convertPhysicalTimeToDromoConstantTimeElement(
                timeFromPropagationStart, currentIndependentVariable, dromoElements, centralBodyGravitationalParameter,
                unitOfLength, useEnergyElement );
    }
    else
    {
        throw std::runtime_error( "Error when converting Cartesian to Dromo(P) elements, specified time type (" +
                std::to_string( timeType ) + ") is not valid." );
    }

    return dromoElements;
}

//! Converts Keplerian to Dromo(P) elements, following Bau (2013) and Bau (2014).
/*!
 * Converts Keplerian to Dromo(P) elements, following Bau (2013) and Bau (2014). Zeta7 is always taken to be negative.
 * \param keplerianElements Vector containing Keplerian elements. Order of elements is important!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination (in range [0,PI]),                            [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 * @param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * @param unitOfLength Unit of length [m] used to make variables dimensionless. Usually it is taken to be the the position
 *          vector norm at the initial state.
 * @param perturbingPotential Potential of perturbing accelerations [m^2/s^2].
 * @param timeFromPropagationStart Time since the start of the propagation [s].
 * @param useEnergyElement Flag indicating whether the dromo state should use energy (true) or zeta3 (false).
 * @param timeType Enum indicating the type of time element which should be used.
 * @param initialIndependentVariable Initial value of the independent variable. Default is TUDAT_NAN, in which case
 *          the initial independent variable is taken to be the true anomaly.
 * @param currentIndependentVariable Current value of the independent variable. Default is TUDAT_NAN, in which case
 *          the current independent variable is taken to be the same as the initial independent variable.
 * @return dromoElements State in Dromo(P) elements.
 *         The order of elements is fixed
 *         dromoElements( 0 ) = Time element,                 [-]
 *         dromoElements( 1 ) = Zeta1 element,                [-]
 *         dromoElements( 2 ) = Zeta2 element,                [-]
 *         dromoElements( 3 ) = Energy or zeta3 element,      [-]
 *         dromoElements( 4 ) = Zeta4 element,                [-]
 *         dromoElements( 5 ) = Zeta5 element,                [-]
 *         dromoElements( 6 ) = Zeta6 element.                [-]
 *         dromoElements( 7 ) = Zeta7 element.                [-]
 */
Eigen::Vector8d convertKeplerianToDromoElements(const Eigen::Vector6d& keplerianElements,
                                                const double centralBodyGravitationalParameter,
                                                const double unitOfLength,
                                                const double perturbingPotential,
                                                const double timeFromPropagationStart,
                                                const bool useEnergyElement,
                                                const propagators::TimeElementType timeType,
                                                double initialIndependentVariable = TUDAT_NAN,
                                                double currentIndependentVariable = TUDAT_NAN)
{
    // Declaring eventual output vector.
    Eigen::Vector8d dromoElements;

    // Check if provided keplerian elements are valid

    using mathematical_constants::PI;

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // If eccentricity is outside range [0,inf)
    if ( keplerianElements( eccentricityIndex ) < 0.0 )
    {
        //Define the error message
        std::stringstream errorMessage;
        errorMessage << "Eccentricity is expected in range [0,inf)\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << std::endl;

        // Throw exception
        throw std::runtime_error(errorMessage.str( ) );
    }

    // If inclination is outside range [0,PI]
    else if ( ( keplerianElements( inclinationIndex ) < 0.0 ) || ( keplerianElements( inclinationIndex ) > PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Inclination is expected in range [0," << PI << "]\n"
                     << "Specified inclination: " << keplerianElements( inclinationIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // If argument of pericenter is outside range [0,2.0 * PI]
    else if ( ( keplerianElements( argumentOfPeriapsisIndex ) < 0.0 ) ||
        ( keplerianElements( argumentOfPeriapsisIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Argument of periapsis is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( argumentOfPeriapsisIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // If right ascension of ascending node is outside range [0,2.0 * PI]
    else if ( ( keplerianElements( longitudeOfAscendingNodeIndex ) < 0.0 ) ||
         ( keplerianElements( longitudeOfAscendingNodeIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "RAAN is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( longitudeOfAscendingNodeIndex ) << " rad."
                     << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // If true anomaly is outside range [0,2.0 * PI]
    else if ( ( keplerianElements( trueAnomalyIndex ) < 0.0 ) ||
        ( keplerianElements( trueAnomalyIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "True anomaly is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( trueAnomalyIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // If inclination is zero and the right ascension of ascending node is non-zero
    else if ( ( std::fabs( keplerianElements( inclinationIndex ) ) < singularityTolerance ) &&
         ( std::fabs( keplerianElements( longitudeOfAscendingNodeIndex ) ) > singularityTolerance ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the inclination is zero, the right ascending node should be zero by definition\n"
                     << "Specified right ascension of ascending node: " <<
                        keplerianElements( longitudeOfAscendingNodeIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // If eccentricity is zero and the argument of pericenter is non-zero
    else if ( ( std::fabs( keplerianElements( eccentricityIndex ) ) < singularityTolerance ) &&
         ( std::fabs( keplerianElements( argumentOfPeriapsisIndex ) ) > singularityTolerance ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the eccentricity is zero, the argument of pericenter should be zero by definition\n"
                     << "Specified argument of pericenter: " <<
                        keplerianElements( argumentOfPeriapsisIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // If semi-major axis is negative and the eccentricity is smaller or equal to one
    else if ( ( keplerianElements( semiMajorAxisIndex ) < 0.0 ) && ( keplerianElements( eccentricityIndex ) <= 1.0 ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the semi-major axis is negative, the eccentricity should be larger than one\n"
                     << "Specified semi-major axis: " << keplerianElements( semiMajorAxisIndex ) << " m.\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }

    // If semi-major axis is positive and the eccentricity is larger than one
    else if ( ( keplerianElements( semiMajorAxisIndex ) > 0.0 ) && ( keplerianElements( eccentricityIndex ) > 1.0 ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the semi-major axis is positive, the eccentricity should be smaller than or equal to one\n"
                     << "Specified semi-major axis: " << keplerianElements( semiMajorAxisIndex ) << " m.\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( errorMessage.str( ) );
    }
    //Else, nothing wrong and continue

    // Check if provided unit of length is valid
    if ( unitOfLength <= 0.0 )
    {
        throw std::runtime_error( "Error when setting converting Keplerian to Dromo(P) elements, unit of length (" +
                std::to_string( unitOfLength ) + ") is <= zero." );
    }

    // If current independent variable isn't defined set values to default ones
    if ( isnan(currentIndependentVariable) )
    {
        if ( isnan(initialIndependentVariable) )
        {
            initialIndependentVariable = keplerianElements( trueAnomalyIndex );
        }
        currentIndependentVariable = initialIndependentVariable;
    }
    else if ( currentIndependentVariable < initialIndependentVariable )
    {
        throw std::runtime_error( "Error when converting Keplerian to Dromo(P) elements, initial independent variable (" +
                std::to_string( initialIndependentVariable ) + ") is larger than current independent variable (" +
                std::to_string( currentIndependentVariable ) + ")." );
    }

    // Extract kepler elements and make them dimensionless (just the semi-major axis)
    double ecc = keplerianElements( eccentricityIndex );
    double inc = keplerianElements( inclinationIndex );
    double raan = keplerianElements( longitudeOfAscendingNodeIndex );
    double aop = keplerianElements( argumentOfPeriapsisIndex );
    double trAnom = keplerianElements( trueAnomalyIndex );

    // Compute dimensionless potential
    double unitOfTime = computeDromoUnitOfTime(centralBodyGravitationalParameter, unitOfLength);
    double dimensionlessPerturbingPotential = perturbingPotential * std::pow( unitOfTime, 2 ) / std::pow( unitOfLength, 2 );

    // Compute dimensionless angular momentum depending on whether orbit is or not parabolic
    double angularMomentumNorm;
    if ( ecc == 1.0 )
    {
        // Parabolic orbit: compute using dimensionless semi latus rectum
        angularMomentumNorm = std::sqrt( keplerianElements( semiLatusRectumIndex ) / unitOfLength );
    }
    else
    {
        // Non-parabolic orbit: compute using dimensionless semi major axis and eccentricity
        angularMomentumNorm = std::sqrt( keplerianElements( semiMajorAxisIndex ) / unitOfLength *
                ( 1 - std::pow( ecc, 2) ) );
    }

    double dimensionlessCartesianPositionNorm =
            convertKeplerianToCartesianElements(keplerianElements, centralBodyGravitationalParameter).segment(0, 3).norm() / unitOfLength;
    double transverseVelocityNorm = angularMomentumNorm / dimensionlessCartesianPositionNorm;
    if ( perturbingPotential > - std::pow( transverseVelocityNorm, 2 ) / 2 )
    {
        throw std::runtime_error(
            "Value of perturbing potential (" + std::to_string( perturbingPotential ) + ") is too large (maximum value is " +
            std::to_string( - std::pow( transverseVelocityNorm, 2 ) / 2 ) + "), which will lead to singularities. For "
            "details see Section 6 of A New Set of Integrals of Motion to Propagate the Perturbed Two-Body Problem, Bau et al. (2013)." );
    }

    // Compute dromo elements
    double eccCosTrAnom = ecc * std::cos( trAnom );
    double term0 = angularMomentumNorm * std::sqrt( std::pow( ( 1 + eccCosTrAnom ), 2 ) +
            2 * dimensionlessPerturbingPotential * std::pow( angularMomentumNorm, 2 ) );
    double term1 = ( eccCosTrAnom * ( 1 + eccCosTrAnom ) + 2 * dimensionlessPerturbingPotential
            * std::pow( angularMomentumNorm, 2) ) / term0;

    dromoElements( dromoPZeta1Index ) = term1 * std::cos( currentIndependentVariable ) + ecc / angularMomentumNorm *
            std::sin( trAnom ) * std::sin( currentIndependentVariable );
    dromoElements( dromoPZeta2Index ) = term1 * std::sin( currentIndependentVariable ) - ecc / angularMomentumNorm *
            std::sin( trAnom ) * std::cos( currentIndependentVariable );

    double zeta3 = ( 1 + eccCosTrAnom ) / term0;
    // Select the dromo element 3 to use based on the flag
    if ( useEnergyElement )
    {
        // Computing energy
        dromoElements( dromoPZeta3Index ) = ( std::pow( dromoElements( dromoPZeta1Index ), 2 ) +
                std::pow( dromoElements( dromoPZeta2Index ), 2 ) - std::pow( zeta3, 2 ) ) / 2;
    }
    else
    {
        dromoElements( dromoPZeta3Index ) = zeta3;
    }

    double dPhi = currentIndependentVariable - initialIndependentVariable;

    // Note on the selection of the sign of [z4,z5,z6,z7], z=zeta:
    // The Dromo(P) equations of motion and element conversions have an ambiguity when it comes to the sign of [z4,z5,z6,z7],
    // i.e. its the same using [z4,z5,z6,z7] or -[z4,z5,z6,z7]. When converting cartesian to Dromo(P) elements, that sign
    // is explicitly selected when computing z7. However, following the equations in Bau(2013), that is not done when
    // converting between Kelplerian and Dromo(P) elements. To ensure consistency in the conversions cartesian->Dromo(P)
    // and keplerian->Dromo(P) zeta7 is always taken to be negative.

    // Compute zeta7 and select it to be negative
    dromoElements( dromoPZeta7Index ) = std::cos( inc / 2 ) * std::cos( ( raan + aop + trAnom - dPhi ) / 2 );
    int signFlip;
    if ( dromoElements( dromoPZeta7Index ) >= 0)
    {
        signFlip = -1;
    }
    else
    {
        signFlip = 1;
    }
    dromoElements( dromoPZeta7Index ) *= signFlip;

    dromoElements( dromoPZeta4Index ) = signFlip * std::sin( inc / 2 ) * std::cos( ( raan - aop - trAnom + dPhi ) / 2 );
    dromoElements( dromoPZeta5Index ) = signFlip * std::sin( inc / 2 ) * std::sin( ( raan - aop - trAnom + dPhi ) / 2 );
    dromoElements( dromoPZeta6Index ) = signFlip * std::cos( inc / 2 ) * std::sin( ( raan + aop + trAnom - dPhi ) / 2 );


    // Compute the time element
    if ( timeType == propagators::scaled_physical_time )
    {
        dromoElements( dromoPTimeIndex ) = convertPhysicalTimeToScaledTime( timeFromPropagationStart,
                                                                            centralBodyGravitationalParameter,
                                                                            unitOfLength );
    }
    else if ( timeType == propagators::linear_time_element )
    {
        dromoElements( dromoPTimeIndex ) = convertPhysicalTimeToDromoLinearTimeElement(
                timeFromPropagationStart, currentIndependentVariable, dromoElements, centralBodyGravitationalParameter,
                unitOfLength, useEnergyElement );
    }
    else if ( timeType == propagators::constant_time_element )
    {
        dromoElements( dromoPTimeIndex ) = convertPhysicalTimeToDromoConstantTimeElement(
                timeFromPropagationStart, currentIndependentVariable, dromoElements, centralBodyGravitationalParameter,
                unitOfLength, useEnergyElement );
    }
    else
    {
        throw std::runtime_error( "Error when converting Keplerian to Dromo(P) elements, specified time type (" +
                std::to_string( timeType ) + ") is not valid." );
    }

    return dromoElements;
}

Eigen::Vector6d convertDromoToCartesianElements(const Eigen::Vector8d dromoElements,
                                                const double centralBodyGravitationalParameter,
                                                const double initialIndependentVariable,
                                                const double currentIndependentVariable,
                                                const double unitOfLength,
                                                const double perturbingPotential,
                                                const bool usingEnergyElement,
                                                const bool forceQuaternionNormalization )
{
    // Declaring eventual output vector.
    Eigen::Vector6d cartesianElements;

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Extract quaternion elements
    Eigen::Vector8d normalizedDromoElements = dromoElements;
    const double normOfQuaternionElements = normalizedDromoElements.segment(dromoPZeta4Index, 4).norm( );
    if ( std::fabs( normOfQuaternionElements - 1.0 ) > singularityTolerance )
    {
        if ( forceQuaternionNormalization )
        {
            normalizedDromoElements.segment(dromoPZeta4Index, 4) /= normOfQuaternionElements;
        }
        else
        {
            // Define the error message.
            std::stringstream errorMessage;
            errorMessage << "Error when converting Dromo(P) to cartesian elements: the norm of the quaternion should be equal to one.\n"
                         << "Norm of the specified quaternion is: " << normOfQuaternionElements - 1.0 << " + 1." << std::endl;

            // Throw exception.
            throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
        }
    }
    // Else, nothing wrong and continue

    double unitOfTime = computeDromoUnitOfTime(centralBodyGravitationalParameter, unitOfLength);
    // Compute dimensionless potential
    double dimensionlessPerturbingPotential = perturbingPotential * std::pow( unitOfTime, 2 ) / std::pow( unitOfLength, 2 );

    double energy;
    double zeta3;
    computeEnergyAndDromoZeta3( normalizedDromoElements, usingEnergyElement, energy, zeta3 );

    const double positionNorm = 1 / ( zeta3 * computeDromoS( normalizedDromoElements, currentIndependentVariable, zeta3 ) );

    const double zeta4sq = std::pow( normalizedDromoElements( dromoPZeta4Index ), 2 );
    const double zeta5sq = std::pow( normalizedDromoElements( dromoPZeta5Index ), 2 );
    const double zeta6sq = std::pow( normalizedDromoElements( dromoPZeta6Index ), 2 );
    const double zeta45 = normalizedDromoElements( dromoPZeta4Index ) * normalizedDromoElements( dromoPZeta5Index );
    const double zeta46 = normalizedDromoElements( dromoPZeta4Index ) * normalizedDromoElements( dromoPZeta6Index );
    const double zeta47 = normalizedDromoElements( dromoPZeta4Index ) * normalizedDromoElements( dromoPZeta7Index );
    const double zeta56 = normalizedDromoElements( dromoPZeta5Index ) * normalizedDromoElements( dromoPZeta6Index );
    const double zeta57 = normalizedDromoElements( dromoPZeta5Index ) * normalizedDromoElements( dromoPZeta7Index );
    const double zeta67 = normalizedDromoElements( dromoPZeta6Index ) * normalizedDromoElements( dromoPZeta7Index );

    Eigen::Matrix3d Q0;
    Q0 << 1 - 2*zeta5sq - 2*zeta6sq, 2*zeta45 - 2*zeta67, 2*zeta46 + 2*zeta57,
          2*zeta45 + 2*zeta67, 1 - 2*zeta4sq - 2*zeta6sq, 2*zeta56 - 2*zeta47,
          2*zeta46 - 2*zeta57, 2*zeta56 + 2 * zeta47, 1 - 2*zeta4sq - 2*zeta5sq;

    Eigen::Matrix3d QRI = Q0 * computeDromoMPhiMatrix( initialIndependentVariable, currentIndependentVariable );

    // Compute dimensionless position and then make it dimensional
    cartesianElements.segment(0, 3) = QRI * (Eigen::Vector3d() << positionNorm, 0, 0).finished() * unitOfLength;

    // Compute dimensionless velocity and then make it dimensionless
    cartesianElements.segment(3, 3) = QRI * (Eigen::Vector3d() <<
            computeDromoRadialVelocity( normalizedDromoElements, currentIndependentVariable ),
            computeDromoTransverseVelocity( computeDromoS( normalizedDromoElements, currentIndependentVariable, zeta3 ),
                                            dimensionlessPerturbingPotential ),
            0).finished() * unitOfLength / unitOfTime;

    return cartesianElements;
}

//Eigen::Vector6d convertDromoToKeplerianElements(const Eigen::Vector8d dromoElements,
//                                                const double centralBodyGravitationalParameter,
//                                                const double initialIndependentVariable,
//                                                const double currentIndependentVariable,
//                                                const double unitOfLength,
//                                                const double perturbingPotential,
//                                                const bool usingEnergyElement)
//{
//    // Declaring eventual output vector.
//    Eigen::Vector6d keplerianElements;
//
//    double unitOfTime = computeDromoUnitOfTime(centralBodyGravitationalParameter, unitOfLength);
//    // Compute dimensionless potential
//    double dimensionlessPerturbingPotential = perturbingPotential * std::pow( unitOfTime, 2 ) / std::pow( unitOfLength, 2 );
//
//    double energy;
//    double zeta3;
//    computeEnergyAndDromoZeta3( dromoElements, usingEnergyElement, energy, zeta3 );
//    double s = computeDromoS( dromoElements, currentIndependentVariable, zeta3 );
//    double radialVelocityNorm = computeDromoRadialVelocity( dromoElements, currentIndependentVariable );
//    double transverseVelocityNorm = computeDromoTransverseVelocity( s, dimensionlessPerturbingPotential );
//
//    keplerianElements( semiMajorAxisIndex ) = - 1.0 / ( 2.0 * energy) * unitOfLength;
//    keplerianElements( eccentricityIndex ) = std::sqrt( std::pow( transverseVelocityNorm, 2.0 ) *
//            ( std::pow( dromoElements(dromoPZeta1Index), 2.0) * std::pow(dromoPZeta2Index, 2.0) - 2.0 * dimensionlessPerturbingPotential ) +
//            2.0 * dimensionlessPerturbingPotential * std::pow( zeta3, 2.0 ) ) / ( zeta3 * s );
//    keplerianElements( inclinationIndex ) = std::acos( 1 - 2.0 * std::pow( dromoElements(dromoPZeta4Index), 2 ) -
//            2.0 * std::pow( dromoElements(dromoPZeta5Index), 2) );
//    keplerianElements( longitudeOfAscendingNodeIndex ) = std::atan2()
//
//}

double convertDromoTimeToPhysicalTime(const Eigen::Vector8d dromoElements,
                                      const double centralBodyGravitationalParameter,
                                      const double currentIndependentVariable,
                                      const double unitOfLength,
                                      const bool usingEnergyElement,
                                      const propagators::TimeElementType timeType)
{
    double timeFromPropagationStart;

    if ( timeType == propagators::scaled_physical_time )
    {
        timeFromPropagationStart = convertScaledTimeToPhysicalTime(
                dromoElements( dromoPTimeIndex ), centralBodyGravitationalParameter, unitOfLength );
    }
    else if ( timeType == propagators::linear_time_element )
    {
        timeFromPropagationStart = convertDromoLinearTimeElementToPhysicalTime(
                dromoElements( dromoPTimeIndex ), currentIndependentVariable, dromoElements, centralBodyGravitationalParameter,
                unitOfLength, usingEnergyElement);
    }
    else if ( timeType == propagators::constant_time_element )
    {
        timeFromPropagationStart = convertDromoConstantTimeElementToPhysicalTime(
                dromoElements( dromoPTimeIndex ), currentIndependentVariable, dromoElements, centralBodyGravitationalParameter,
                unitOfLength, usingEnergyElement);
    }
    else
    {
        throw std::runtime_error( "Error when setting converting Dromo(P) time to physical time, specified time type (" +
                std::to_string( timeType ) + ") is not valid." );
    }

    return timeFromPropagationStart;
}

} // namespace orbital_element_conversions

} // close tudat

#endif //TUDAT_DROMO_ELEMENT_CONVERSIONS_H
