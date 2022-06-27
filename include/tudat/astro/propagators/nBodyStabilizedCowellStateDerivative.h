/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
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
//! @param currentCartesianPosition Current cartesian position
//! @param sundmanConstant Multiplying constant used in the Sundman transformation
//! @return Derivative of the physical time with respect to the fictitious time
double computePhysicalTimeDerivativeForStabilizedCowell(
            const Eigen::Vector3d& currentCartesianPosition,
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

//! Function to compute the difference between the linear time element and physical time
//! \param currentStabilizedCowellState Current state in stabilized Cowell formulation of the body for which the equations of
//! motions are to be evaluated
//! \param sundmanConstant Multiplying constant used in the Sundman transformation
//! \return Difference between the linear time element and physical time
double computeLinearTimeElementToPhysicalTimeBias (
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant);

//! Function to convert the linear time element to physical time
//! \param currentStabilizedCowellState Current state in stabilized Cowell formulation of the body for which the equations of
//! motions are to be evaluated
//! \param sundmanConstant Multiplying constant used in the Sundman transformation
//! \return Physical time
double convertLinearTimeElementToPhysicalTime (
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant);

//! Function to convert the physical time to a linear time element
//! \param currentStabilizedCowellState Current state in stabilized Cowell formulation of the body for which the equations of
//! motions are to be evaluated
//! \param sundmanConstant Multiplying constant used in the Sundman transformation
//! \return Linear time element
double convertPhysicalTimeToLinearTimeElement (
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

//! Computes the energy for the current cartesian state.
//! \param centralBodyGravitationalParameter Gravitational parameter of sum of central body and body for which orbit is propagated.
//! \param currentCartesianState State in 'conventional form'
//! \param conservativeAccelerationsPotential Total potential of conservative accelerations
//! \return Energy
double computeEnergy(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector6d& currentCartesianState,
        const double conservativeAccelerationsPotential);

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
        centralBodyGravitationalParameterFunction_ = removeCentralGravityAccelerations(
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
        sundmanConstant_ = std::sqrt(semiMajorAxis / std::sqrt(centralBodyGravitationalParameterFunction_( ) ) );
    }

    //! Destructor
    ~NBodyStabilizedCowellStateDerivative( ){ }

    //! Calculates the state derivative of the translational motion of the system, using the equations of motion for
    //! stabilized Cowell.
    /*!
     *  Calculates the state derivative of the translational motion of the system, using the equations of motion for
     *  stabilized Cowell. The input is the current state in stabilized cowell elements. The state derivative
     *  of this set is computed.
     *  \param time Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 8 * bodiesToBeIntegratedNumerically_.size( ), containing stabilized cowell
     *  elements of the bodies being integrated.
     *  The order of the values is defined by the order of bodies in bodiesToBeIntegratedNumerically_
     *  \param stateDerivative Current derivative of the stabilized cowell elements of the
     *  system of bodies integrated numerically (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        std::cerr << time << std::endl;
        // Get total inertial accelerations acting on bodies
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > totalAccelerationInInertialFrame;
        totalAccelerationInInertialFrame.resizeLike(currentCartesianLocalSolution_ );
        this->sumStateDerivativeContributions(stateOfSystemToBeIntegrated, totalAccelerationInInertialFrame, false );

        // TODO: Computation of conservative accelerations
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > conservativeAccelerationInInertialFrame;
        conservativeAccelerationInInertialFrame.resizeLike(currentCartesianLocalSolution_ );
        conservativeAccelerationInInertialFrame.setZero();

        Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > nonConservativeAccelerationInInertialFrame =
                totalAccelerationInInertialFrame - conservativeAccelerationInInertialFrame;

        // TODO: computation of potential of conservative accelerations
        double conservativeAccelerationsPotential = 0.0;

        // Evaluate equations of motion
        stateDerivative.setZero( );
        stateDerivative.block( 0, 0, getPropagatedStateSize() - 1, 1 ) = computeStateDerivativeExceptTimeForStabilizedCowell(
                stateOfSystemToBeIntegrated, sundmanConstant_, conservativeAccelerationsPotential,
                nonConservativeAccelerationInInertialFrame, conservativeAccelerationInInertialFrame);
        // Evaluate time derivative
        if ( timeType_ == physical_time )
        {
            stateDerivative(orbital_element_conversions::stabilizedCowellTimeIndex, 0 ) = computePhysicalTimeDerivativeForStabilizedCowell(
                    stateOfSystemToBeIntegrated.block(0, 0, 3, 1), sundmanConstant_);
        }
        else // Linear time element
        {
            stateDerivative(orbital_element_conversions::stabilizedCowellTimeIndex, 0 ) = computeLinearTimeElementDerivativeForStabilizedCowell(
                    stateOfSystemToBeIntegrated, sundmanConstant_,
                    stateDerivative(orbital_element_conversions::stabilizedCowellEnergyIndex, 0 ), centralBodyGravitationalParameterFunction_(),
                    nonConservativeAccelerationInInertialFrame);
        }
//        std::cerr << "State: " << stateOfSystemToBeIntegrated.transpose() << std::endl;
//        std::cerr << "State derivative: " << stateDerivative.transpose() << std::endl;

    }

    //! Function to convert the stabilized-cowell states of the bodies to the conventional form.
    /*!
     * Function to convert the stabilized-cowell elements state to the conventional form. For the stabilized-cowell
     * propagator, this transforms stabilized-cowell elements w.r.t. the central bodies to the Cartesian states w.r.t. these
     * same central bodies: In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in stabilized cowell elements (i.e. form that is used in
     * numerical integration)
     * \param independentVariable Current independent variable at which the state is valid
     * \param currentCartesianLocalSolution State (internalSolution, which is stabilized-cowell formulation),
     * converted to the 'conventional form' (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& independentVariable,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSolution )
    {
        currentCartesianLocalSolution.block( 0, 0, 3, 1 ) = internalSolution.block( 0, 0, 3, 1 );
        currentCartesianLocalSolution.block( 3, 0, 3, 1 ) = internalSolution.block( 3, 0, 3, 1 ) /
                 computePhysicalTimeDerivativeForStabilizedCowell( internalSolution.block(0, 0, 3, 1), sundmanConstant_);
        currentCartesianLocalSolution_ = currentCartesianLocalSolution;
    }

    //! Function to convert the state in the conventional form to the stabilized-cowell form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For the stabilized-cowell propagator,
     * this transforms the Cartesian state w.r.t. the central body (conventional form) to the stabilized-cowell elements
     * \param cartesianSolution State in 'conventional form'
     * \param time Current time at which the state is valid.
     * \return State (outputSolution), converted to the stabilized-cowell elements
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& physicalTime )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( getPropagatedStateSize( ) );

        // Position
        currentState.block( 0, 0, 3, 1 ) = cartesianSolution.block( 0, 0, 3, 1 );

        // Velocity
        currentState.block( 3, 0, 3, 1 ) = cartesianSolution.block( 3, 0, 3, 1 ) *
                computePhysicalTimeDerivativeForStabilizedCowell( cartesianSolution.block(0, 0, 3, 1), sundmanConstant_);

        // Energy
        // TODO: computation of potential of conservative accelerations
        double conservativeAccelerationsPotential = 0.0;
        currentState(orbital_element_conversions::stabilizedCowellEnergyIndex, 0 ) =
                computeEnergy( centralBodyGravitationalParameterFunction_(), cartesianSolution, conservativeAccelerationsPotential);

        if ( timeType_ == physical_time )
        {
            currentState(orbital_element_conversions::stabilizedCowellTimeIndex, 0 ) = physicalTime;
        }
        else // Linear time element
        {
            currentState(orbital_element_conversions::stabilizedCowellTimeIndex, 0 ) = convertPhysicalTimeToLinearTimeElement(
                    currentState, sundmanConstant_);
        }

        std::cerr << currentState.transpose() << std::endl;

        return currentState;
    }

    //! Function to return the size of the state handled by the object.
    /*!
     * Function to return the size of the state handled by the object.
     * \return Size of the state under consideration (8 times the number if integrated bodies, which is 1).
     */
    int getPropagatedStateSize( )
    {
        return 8 * this->bodiesToBeIntegratedNumerically_.size( );
    }

    //! Function to return the physical time.
    //! \param independentVariable Independent variable.
    //! \return Physical time.
    virtual TimeType getPhysicalTime(
            const TimeType& independentVariable,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated )
    {
        if ( timeType_ == physical_time )
        {
            return stateOfSystemToBeIntegrated( orbital_element_conversions::stabilizedCowellTimeIndex );
        }
        else // Linear time element
        {
            return convertLinearTimeElementToPhysicalTime( stateOfSystemToBeIntegrated, sundmanConstant_ );
        }
    }

    //! Function returns whether the time is part of the state (i.e. whether it is a dependent variable).
    //! \return Boolean informing whether state derivative includes time.
    bool timeIsADependentVariable( )
    {
        return true;
    }

   // Function to compute the initial independent variable.
   /*
    * Function to compute the initial independent variable to use in the propagation. For stabilized cowell, the initial
    * independent variable is a fictitious angle, thus its initial value is considered to be 0.
    * \param initialOutputSolution Initial state in 'conventional form'
    * \param initialTime Initial propagation time, at which the state is valid.
    * \return Independent variable
    */
    virtual TimeType computeInitialIndependentVariable (
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialOutputSolution,
            const TimeType& initialTime)
    {
        return 0.0;
    }

private:

    //! Gravitational parameter of the central body used to convert Cartesian to Keplerian orbits, and vice versa
    std::function< double( ) > centralBodyGravitationalParameterFunction_;

    //! Type of time that is to be integrated
    RegularizedPropagatorTimeType timeType_;

    //! Multiplying constant used in the Sundman transformation
    double sundmanConstant_;

    //! Current full Cartesian state of the propagated bodies, w.r.t. the central bodies
    /*!
     *  Current full Cartesian state of the propagated bodies, w.r.t. the central bodies. These variables are set when calling
     *  the convertToOutputSolution function.
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentCartesianLocalSolution_;
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
