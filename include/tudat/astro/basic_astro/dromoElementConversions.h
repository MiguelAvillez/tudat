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

double computeDromoTimeToTimeElementBias(const Eigen::Vector8d& dromoElementsExceptTime,
                                         const double currentIndependentVariable,
                                         const bool usingEnergyElement)
{
    double energy;
    double dromoElement3;

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

    double s = dromoElement3 + dromoElementsExceptTime( dromoElement1Index ) * std::cos( currentIndependentVariable) +
            dromoElementsExceptTime( dromoElement2Index ) * std::sin( currentIndependentVariable );
    double radialVelocity = dromoElementsExceptTime( dromoElement1Index ) * std::sin( currentIndependentVariable) -
            dromoElementsExceptTime( dromoElement2Index ) * std::cos( currentIndependentVariable );

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


Eigen::Vector8d convertCartesianToDromoElementsExceptTime(const Eigen::Vector6d& cartesianElements,
                                                          const double centralBodyGravitationalParameter,
                                                          const double unitOfLength,
                                                          const double initialIndependentVariable,
                                                          const double currentIndependentVariable,
                                                          const double perturbingPotential,
                                                          const bool useEnergyElement)
{
    // Declaring eventual output vector.
    Eigen::Vector8d dromoElements;
    dromoElements( 0 ) = TUDAT_NAN;

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

    dromoElements( dromoElement4Index ) = ( Q0(3,2) - Q0(2,3) ) / ( 4 * dromoElements( dromoElement7Index ) );
    dromoElements( dromoElement5Index ) = ( Q0(1,3) - Q0(3,1) ) / ( 4 * dromoElements( dromoElement7Index ) );
    dromoElements( dromoElement6Index ) = ( Q0(2,1) - Q0(1,2) ) / ( 4 * dromoElements( dromoElement7Index ) );

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

    return dromoElements;
}

} // namespace orbital_element_conversions

} // close tudat

#endif //TUDATBUNDLE_DROMOELEMENTCONVERSIONS_H
