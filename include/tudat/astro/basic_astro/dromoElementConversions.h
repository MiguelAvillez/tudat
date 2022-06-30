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
 */

#ifndef TUDATBUNDLE_DROMOELEMENTCONVERSIONS_H
#define TUDATBUNDLE_DROMOELEMENTCONVERSIONS_H

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
 * Computes the unit of length used to make variables dimensionless in the Dromo propagator. The unit of length is
 * selected to be the inital radial position. Bau (2013), section 2.
 *
 * @param initialCartesianPosition Initial cartesian position.
 * @return Unit of length.
 */
double computeDromoUnitOfLength ( const Eigen::Vector3d& initialCartesianPosition )
{
    return initialCartesianPosition.norm();
}

//! Computes the unit of time used in the Dromo propagator.
/*!
 * Computes the unit of time used to make variables dimensionless in the Dromo propagator. The unit of time is
 * selected to be angular rate of a circular orbit with the unit of length as radius. Bau (2013), section 2.
 *
 * @param initialCartesianPosition Initial cartesian position.
 * @return Unit of time.
 */
double computeDromoUnitOfTime (const double centralBodyGravitationalParameter,
                               const double unitOfLength)
{
    return std::sqrt( centralBodyGravitationalParameter / std::pow(unitOfLength, 3) );
}

//! Computes the MPhi matrix used in conversions for the Dromo propagator.
/*!
 * Computes the MPhi matrix used when converting the Dromo elements to cartesian and vice-versa. Eq. 50 of Bau (2013).
 *
 * @param initialIndependentVariable Initial independent variable.
 * @param currentIndependentVariable Independent variable for which the matrix is to be computed.
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

//! Convert the time to the time used in Dromo.
/*!
 * Convert the time to the time used in Dromo, i.e. makes the time dimensionless.
 *
 * @param time Dimensional time.
 * @param centralBodyGravitationalParameter Central body gravitational parameter.
 * @param unitOfLength Dromo unit of length.
 * @return Dromo time.
 */
double convertTimeToDromoTime(const double time,
                              const double centralBodyGravitationalParameter,
                              const double unitOfLength)
{
    return time / computeDromoUnitOfTime(centralBodyGravitationalParameter, unitOfLength);
}

double convertDromoTimeToTime(const double dimensionlessDromoTime,
                              const double centralBodyGravitationalParameter,
                              const double unitOfLength)
{
    return dimensionlessDromoTime * computeDromoUnitOfTime(centralBodyGravitationalParameter, unitOfLength);
}

void computeEnergyAndDromoElement3(const Eigen::Vector8d& dromoElementsExceptTime,
                                   const bool usingEnergyElement,
                                   double& energy,
                                   double& dromoElement3)
{
    if ( usingEnergyElement )
    {
        energy = dromoElementsExceptTime( dromoEnergyIndex );
        dromoElement3 = std::sqrt( std::pow( dromoElementsExceptTime( dromoElement1Index), 2 ) +
                std::pow( dromoElementsExceptTime( dromoElement2Index), 2 ) - 2 * energy );
    }
    else
    {
        dromoElement3 = dromoElementsExceptTime ( dromoElement3Index );
        energy = ( std::pow( dromoElementsExceptTime( dromoElement1Index), 2 ) +
                std::pow( dromoElementsExceptTime( dromoElement2Index), 2 ) -
                std::pow( dromoElementsExceptTime( dromoElement3Index), 2 ) ) / 2;
    }
}

//! Compute the Dromo term s.
/*!
 * Compute the Dromo term s, Eq. 40 of Bau (2013).
 *
 * @param dromoElements Dromo state
 * @param currentIndependentVariable Current independent variable.
 * @param dromoElement3 Value of the zeta3 Dromo element.
 * @return s
 */
double computeDromoS(const Eigen::Vector8d& dromoElements,
                     const double currentIndependentVariable,
                     const double dromoElement3)
{
    return dromoElement3 + dromoElements( dromoElement1Index ) * std::cos( currentIndependentVariable) +
           dromoElements( dromoElement2Index ) * std::sin( currentIndependentVariable );
}

//! Compute the radial velocity.
/*!
 * Compute the radial velocity, Eq. 41 of Bau (2013).
 *
 * @param dromoElements Dromo state
 * @param currentIndependentVariable Current independent variable.
 * @return Radial velocity
 */
double computeDromoRadialVelocity (const Eigen::Vector8d& dromoElements,
                                   const double currentIndependentVariable)
{
    return dromoElements( dromoElement1Index ) * std::sin( currentIndependentVariable) -
            dromoElements( dromoElement2Index ) * std::cos( currentIndependentVariable );
}

//! Compute the transverse velocity.
/*!
 * Compute the transverse velocity., Eq. 42 of Bau (2013).
 *
 * @param s Dromo s term.
 * @param perturbingPotential Potential of perturbing accelerations.
 * @return Transverse velocity
 */
double computeDromoTransverseVelocity (const double s,
                                       const double perturbingPotential)
{
    return std::sqrt( std::pow( s, 2 ) - 2 * perturbingPotential );
}

double computeDromoTimeToTimeElementBias(const Eigen::Vector8d& dromoElementsExceptTime,
                                         const double currentIndependentVariable,
                                         const bool usingEnergyElement)
{
    double energy;
    double dromoElement3;
    computeEnergyAndDromoElement3(dromoElementsExceptTime, usingEnergyElement, energy, dromoElement3);

    double s = computeDromoS( dromoElementsExceptTime, currentIndependentVariable, dromoElement3 );
    double radialVelocity = computeDromoRadialVelocity( dromoElementsExceptTime, currentIndependentVariable );

    return - radialVelocity / ( 2 * energy * dromoElement3 * s ) - 1 / ( energy * std::sqrt( -2 * energy) ) *
        std::atan2( radialVelocity, s + std::sqrt( -2 * energy) );
}

double convertTimeToDromoLinearTimeElement(const double time,
                                           const double currentIndependentVariable,
                                           const Eigen::Vector8d& dromoElementsExceptTime,
                                           const double centralBodyGravitationalParameter,
                                           const double unitOfLength,
                                           const bool usingEnergyElement)
{
    return convertTimeToDromoTime( time, centralBodyGravitationalParameter, unitOfLength ) +
            computeDromoTimeToTimeElementBias( dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement);
}

double convertDromoLinearTimeElementToTime(const double linearTimeElement,
                                           const double currentIndependentVariable,
                                           const Eigen::Vector8d& dromoElementsExceptTime,
                                           const double centralBodyGravitationalParameter,
                                           const double unitOfLength,
                                           const bool usingEnergyElement)
{
    const double dimensionlessDromoTime = linearTimeElement - computeDromoTimeToTimeElementBias(
            dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement);
    return convertDromoTimeToTime( dimensionlessDromoTime, centralBodyGravitationalParameter, unitOfLength );
}

double computeDromoLinearToConstantTimeElementBias(const Eigen::Vector8d& dromoElementsExceptTime,
                                                   const double currentIndependentVariable,
                                                   const bool usingEnergyElement)
{
    double energy;
    double dromoElement3;
    computeEnergyAndDromoElement3(dromoElementsExceptTime, usingEnergyElement, energy, dromoElement3);

    double semiMajorAxis = -1 / ( 2 * energy );

    return - std::pow( semiMajorAxis, 3/2 ) * currentIndependentVariable;
}

double convertTimeToDromoConstantTimeElement(const double time,
                                             const double currentIndependentVariable,
                                             const Eigen::Vector8d& dromoElementsExceptTime,
                                             const double centralBodyGravitationalParameter,
                                             const double unitOfLength,
                                             const bool usingEnergyElement)
{
    double linearTimeElement = convertTimeToDromoLinearTimeElement(
            time, currentIndependentVariable, dromoElementsExceptTime, centralBodyGravitationalParameter, unitOfLength,
            usingEnergyElement);

    return linearTimeElement + computeDromoLinearToConstantTimeElementBias(
            dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement);
}

double convertDromoConstantTimeElementToTime(const double constantTimeElement,
                                             const double currentIndependentVariable,
                                             const Eigen::Vector8d& dromoElementsExceptTime,
                                             const double centralBodyGravitationalParameter,
                                             const double unitOfLength,
                                             const bool usingEnergyElement)
{
    const double linearTimeElement = constantTimeElement - computeDromoLinearToConstantTimeElementBias(
            dromoElementsExceptTime, currentIndependentVariable, usingEnergyElement);
    return convertDromoLinearTimeElementToTime(
            linearTimeElement, currentIndependentVariable, dromoElementsExceptTime, centralBodyGravitationalParameter,
            unitOfLength, usingEnergyElement);
}

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

    // Auxiliary variables
    double s = std::sqrt( std::pow( transverseVelocityNorm, 2 ) + 2 * dimensionlessPerturbingPotential );
    double dromoElement3 = 1 / s;

    Eigen::Matrix3d QRI;
    QRI << cartesianPosition.transpose() / cartesianPositionNorm,
        angularMomentum.cross( cartesianPosition ).transpose() / angularMomentumNorm / cartesianPositionNorm,
        angularMomentum.transpose() / angularMomentumNorm;
    Eigen::Matrix3d Q0 = QRI * computeDromoMPhiMatrix( initialIndependentVariable, currentIndependentVariable ).transpose();

    // Compute dromo elements
    dromoElements( dromoElement1Index ) = ( s - dromoElement3 ) * std::cos( currentIndependentVariable ) +
            radialVelocityNorm * std::sin( currentIndependentVariable );
    dromoElements( dromoElement2Index ) = ( s - dromoElement3 ) * std::sin( currentIndependentVariable ) -
            radialVelocityNorm * std::sin( currentIndependentVariable );

    dromoElements( dromoElement7Index ) = -0.5 * std::sqrt( 1 + Q0(1,1) + Q0(2,2) + Q0(3,3) );

    const double tolerance = 1e-10;
    if ( std::abs( dromoElements( dromoElement7Index ) ) > tolerance )
    {
        dromoElements( dromoElement4Index ) = ( Q0(3,2) - Q0(2,3) ) / ( 4 * dromoElements( dromoElement7Index ) );
        dromoElements( dromoElement5Index ) = ( Q0(1,3) - Q0(3,1) ) / ( 4 * dromoElements( dromoElement7Index ) );
        dromoElements( dromoElement6Index ) = ( Q0(2,1) - Q0(1,2) ) / ( 4 * dromoElements( dromoElement7Index ) );
    }
    else
    {
        dromoElements( dromoElement6Index ) = std::sqrt( ( Q0(3,3) + 1 ) / 2 );
        if ( std::abs( dromoElements( dromoElement6Index ) ) > tolerance )
        {
            dromoElements( dromoElement4Index ) = Q0(1,3) / ( 2 * dromoElements( dromoElement6Index ) );
            dromoElements( dromoElement5Index ) = Q0(2,3) / ( 2 * dromoElements( dromoElement6Index ) );
        }
        else
        {
            dromoElements( dromoElement4Index ) = std::sqrt( ( 1 - Q0(2,2) ) / 2 );
            if ( std::abs( dromoElements( dromoElement4Index ) ) > tolerance )
            {
                dromoElements( dromoElement5Index ) = Q0(1,2) / ( 2 * dromoElements( dromoElement4Index ) );
            }
            else
            {
                dromoElements( dromoElement5Index ) = 1;
            }
        }
    }

    // Select the dromo element 3 to use based on the flag
    if ( useEnergyElement )
    {
        // Computing energy: the second term is -1/r instead of -mu/r because the position and velocity are dimensionless
        dromoElements( dromoElement3Index ) = std::pow( cartesianVelocity.norm(), 2 ) / 2 -
                1 / cartesianPositionNorm + dimensionlessPerturbingPotential;
    }
    else
    {
        dromoElements( dromoElement3Index ) = dromoElement3;
    }

    if ( timeType == propagators::physical_time)
    {
        dromoElements( dromoTimeIndex ) = convertTimeToDromoTime( timeFromPropagationStart, centralBodyGravitationalParameter,
                                                                  unitOfLength );
    }
    else if ( timeType == propagators::linear_time_element )
    {
        dromoElements( dromoTimeIndex ) = convertTimeToDromoLinearTimeElement(
                timeFromPropagationStart, currentIndependentVariable, dromoElements, centralBodyGravitationalParameter,
                unitOfLength, useEnergyElement);
    }
    else if ( timeType == propagators::constant_time_element )
    {
        dromoElements( dromoTimeIndex ) = convertTimeToDromoConstantTimeElement(
                timeFromPropagationStart, currentIndependentVariable, dromoElements, centralBodyGravitationalParameter,
                unitOfLength, useEnergyElement);
    }

    return dromoElements;
}

Eigen::Vector6d convertDromoToCartesianElements(const Eigen::Vector8d dromoElements,
                                                const double centralBodyGravitationalParameter,
                                                const double initialIndependentVariable,
                                                const double currentIndependentVariable,
                                                const double unitOfLength,
                                                const double perturbingPotential,
                                                const bool usingEnergyElement)
{
    // Declaring eventual output vector.
    Eigen::Vector6d cartesianElements;

    double energy;
    double dromoElement3;
    computeEnergyAndDromoElement3( dromoElements, usingEnergyElement, energy, dromoElement3 );

    const double positionNorm = 1 / ( dromoElement3 * computeDromoS( dromoElements, currentIndependentVariable, dromoElement3 ) );

    const double zeta4sq = std::pow( dromoElements( dromoElement4Index ), 2 );
    const double zeta5sq = std::pow( dromoElements( dromoElement5Index ), 2 );
    const double zeta6sq = std::pow( dromoElements( dromoElement6Index ), 2 );
    const double zeta45 = dromoElements( dromoElement4Index ) * dromoElements( dromoElement5Index );
    const double zeta46 = dromoElements( dromoElement4Index ) * dromoElements( dromoElement6Index );
    const double zeta47 = dromoElements( dromoElement4Index ) * dromoElements( dromoElement7Index );
    const double zeta56 = dromoElements( dromoElement5Index ) * dromoElements( dromoElement6Index );
    const double zeta57 = dromoElements( dromoElement5Index ) * dromoElements( dromoElement7Index );
    const double zeta67 = dromoElements( dromoElement6Index ) * dromoElements( dromoElement7Index );

    Eigen::Matrix3d Q0;
    Q0 << 1 - 2*zeta5sq - 2*zeta6sq, 2*zeta45 - 2*zeta67, 2*zeta46 + 2*zeta57,
          2*zeta45 + 2*zeta67, 1 - 2*zeta4sq - 2*zeta6sq, 2*zeta56 - 2*zeta47,
          2*zeta46 - 2*zeta57, 2*zeta56 + 2 * zeta47, 1 - 2*zeta4sq - 2*zeta5sq;

    Eigen::Matrix3d QRI = Q0 * computeDromoMPhiMatrix( initialIndependentVariable, currentIndependentVariable );

    // Compute dimensionless position and then make it dimensional
    cartesianElements.segment(0, 3) = QRI * (Eigen::Vector3d() << positionNorm, 0, 0).finished() * unitOfLength;

    // Compute dimensionless velocity and then make it dimensionless
    cartesianElements.segment(3, 3) = QRI * (Eigen::Vector3d() <<
            computeDromoRadialVelocity( dromoElements, currentIndependentVariable ),
            computeDromoTransverseVelocity( computeDromoS( dromoElements, currentIndependentVariable, dromoElement3 ), perturbingPotential ),
            0).finished() * unitOfLength / computeDromoUnitOfTime( centralBodyGravitationalParameter, unitOfLength );

    return cartesianElements;
}

} // namespace orbital_element_conversions

} // close tudat

#endif //TUDATBUNDLE_DROMOELEMENTCONVERSIONS_H
