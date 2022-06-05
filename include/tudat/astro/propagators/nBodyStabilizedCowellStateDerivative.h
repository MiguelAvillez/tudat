/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATBUNDLE_NBODYSTABILIZEDCOWELLSTATEDERIVATIVE_H
#define TUDATBUNDLE_NBODYSTABILIZEDCOWELLSTATEDERIVATIVE_H

#include "tudat/astro/propagators/nBodyStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to compute the derivative of the physical time with respect to the fictitious time (Janin, 1974, Eq. 3.10).
//! @param currentStabilizedCowellState Current state in stabilized Cowell formulation of the body for which the equations of
//! motions are to be evaluated
//! @param sundmanConstant Multiplying constant used in the Sundman transformation
//! @return Derivative of the physical time with respect to the fictitious time
double computePhysicalTimeDerivativeForStabilizedCowell(
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant );

//! Function to compute the derivative of the linear time element with respect to the fictitious time (Janin, 1974, Eq. 3.4 and 3.11).
//! \param currentStabilizedCowellState Current state in stabilized Cowell formulation of the body for which the equations of
//! motions are to be evaluated
//! \param sundmanConstant Multiplying constant used in the Sundman transformation
//! \param energyDerivative Derivative of the enrgy with respect to the fictitious time
//! \param centralBodyGravitationalParameter Gravitational parameter of sum of central body and body for which orbit is propagated.
//! \param nonConservativeAccelerationsInInertialFrame Vector with non-conservative perturbing accelerations written in inertial frame
//! \return Derivative of the linear time element with respect to the fictitious time
double computeLinearTimeElementDerivativeForStabilizedCowell(
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant,
        const double energyDerivative,
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& nonConservativeAccelerationsInInertialFrame);

//! Function to convert the linear time element to physical time
//! \param currentStabilizedCowellState Current state in stabilized Cowell formulation of the body for which the equations of
//////! motions are to be evaluated
//! \param sundmanConstant Multiplying constant used in the Sundman transformation
//! \return Physical time
double convertLinearTimeElementToPhysicalTime (
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant);

//! Function to compute the derivative of the stabilized cowell elements, except time, with respect to the fictitious time (Janin, 1974, Eq. 3.14).
//! \param currentStabilizedCowellState Current state in stabilized Cowell formulation of the body for which the equations of
//! motions are to be evaluated
//! \param sundmanConstant Multiplying constant used in the Sundman transformation
//! \param conservativeAccelerationsPotential Total potential of conservative accelerations
//! \param nonConservativeAccelerationsInInertialFrame Vector with non-conservative perturbing accelerations written in inertial frame
//! \param conservativeAccelerationsInInertialFrame Vector with conservative perturbing accelerations (i.e. the accelerations
//! associated with the specified conservativeAccelerationsPotential) written in inertial frame
//! \return Derivative of the stabilized cowell elements, except time, with respect to the fictitious time
Eigen::Vector7d computeStateDerivativeExceptTimeForStabilizedCowell(
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant,
        const double conservativeAccelerationsPotential,
        const Eigen::Vector3d& nonConservativeAccelerationsInInertialFrame,
        const Eigen::Vector3d& conservativeAccelerationsInInertialFrame);

//! Class for computing the state derivative of translational motion of N bodies, using a stabilized Cowell propagator.
template< typename StateScalarType = double, typename TimeType = double >
class NBodyStabilizedCowellStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on.  This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     */
    NBodyStabilizedCowellStateDerivative(
            const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
            const std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
            const std::vector< std::string >& bodiesToIntegrate,
            const RegularizedPropagatorTimeType timeType,
            const std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& initialKeplerElements ):
        NBodyStateDerivative< StateScalarType, TimeType >(
                accelerationModelsPerBody, centralBodyData, stabilized_cowell, bodiesToIntegrate ),
        timeType_( timeType )
    {
        // Check if specified timeType is valid
        if ( timeType_ != physical_time and timeType_ != linear_time_element )
        {
            throw std::runtime_error( "Error when setting N Body Stabilized Cowell propagator, specified time type (" +
                std::to_string( timeType_ ) + ") is not implemented." );
        }

        // Check if a single body is to be integrated
        if ( bodiesToIntegrate.size() != 1 )
        {
            throw std::runtime_error( "Error when setting N Body Stabilized Cowell propagator, only possible to integrate a single body ("
                + std::to_string( bodiesToIntegrate.size() ) + " bodies specified)." );
        }

        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameter_ = removeCentralGravityAccelerations(
                centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                this->accelerationModelsPerBody_, this->removedCentralAccelerations_ ).at( 0 );
        this->createAccelerationModelList( );

        // Compute Sundman transformation multiplying term
        double semiMajorAxis = initialKeplerElements.at( 0 )( orbital_element_conversions::semiMajorAxisIndex );
        if ( semiMajorAxis <= 0 )
        {
            throw std::runtime_error( "Error when setting N Body Stabilized Cowell propagator, only possible to integrate orbit with semi-major axis (" +
                std::to_string( semiMajorAxis ) + ") larger than zero." );
        }
        sundmanConstant_ = std::sqrt(semiMajorAxis / std::sqrt(centralBodyGravitationalParameter_( ) ) );
    }

    //! Destructor
    ~NBodyStabilizedCowellStateDerivative( ){ }

private:

    //! Gravitational parameter of the central body used to convert Cartesian to Keplerian orbits, and vice versa
    std::function< double( ) > centralBodyGravitationalParameter_;

    //! Type of time that is to be integrated
    RegularizedPropagatorTimeType timeType_;

    //! Multiplying constant used in the Sundman transformation
    double sundmanConstant_;

};


extern template class NBodyStabilizedCowellStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class NBodyStabilizedCowellStateDerivative< long double, double >;
extern template class NBodyStabilizedCowellStateDerivative< double, Time >;
extern template class NBodyStabilizedCowellStateDerivative< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat

#endif //TUDATBUNDLE_NBODYSTABILIZEDCOWELLSTATEDERIVATIVE_H
