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

#ifndef TUDATBUNDLE_NBODYDROMOSTATEDERIVATIVE_H
#define TUDATBUNDLE_NBODYDROMOSTATEDERIVATIVE_H

#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/basic_astro/dromoElementConversions.h"

namespace tudat
{

namespace propagators
{

Eigen::Vector8d computeStateDerivativeForDromoP(
        const Eigen::Vector8d& currentDromoState,
        const double perturbingPotential,
        const double perturbingPotentialDerivativeWrtDromoTime,
        const double perturbingPotentialDerivativeWrtDromoElement3,
        const Eigen::Vector3d& perturbingNonPotentialAccelerationsInRswFrame,
        const Eigen::Vector3d& perturbingPotentialAccelerationsInRswFrame,
        const double initialIndependentVariable,
        const double currentIndependentVariable,
        const bool usingEnergyElement,
        const TimeElementType timeType )
{
    Eigen::Vector8d stateDerivative;

    double energy;
    double zeta3;
    orbital_element_conversions::computeEnergyAndDromoElement3( currentDromoState, usingEnergyElement, energy, zeta3 );

    double s = orbital_element_conversions::computeDromoS( currentDromoState, currentIndependentVariable, zeta3 );
    double radialVelocity = orbital_element_conversions::computeDromoRadialVelocity( currentDromoState, currentIndependentVariable );
    double transverseVelocity = orbital_element_conversions::computeDromoTransverseVelocity( s, perturbingPotential );

    double energyDerivative = 1 / ( zeta3 * std::pow( s, 2 ) ) * ( radialVelocity * perturbingNonPotentialAccelerationsInRswFrame( 0 ) +
            perturbingNonPotentialAccelerationsInRswFrame( 1 ) * transverseVelocity + perturbingPotentialDerivativeWrtDromoTime );
    double zeta3Derivative = zeta3 / std::pow( s, 2 ) * ( radialVelocity / s *
            ( zeta3 * s / ( s + zeta3 ) * perturbingPotentialDerivativeWrtDromoElement3 +
            perturbingNonPotentialAccelerationsInRswFrame( 2 ) / ( zeta3 * s ) - 2 * perturbingPotential) - energyDerivative );

    if ( usingEnergyElement )
    {
        stateDerivative( orbital_element_conversions::dromoElement3Index ) = energyDerivative;
    }
    else
    {
        stateDerivative( orbital_element_conversions::dromoElement3Index ) = zeta3Derivative;
    }

    double sinPhi = std::sin( currentIndependentVariable );
    double cosPhi = std::cos( currentIndependentVariable );
    double sinDeltaPhi = std::sin( currentIndependentVariable - initialIndependentVariable );
    double cosDeltaPhi = std::cos( currentIndependentVariable - initialIndependentVariable );
    double zeta4 = currentDromoState( orbital_element_conversions::dromoElement4Index );
    double zeta5 = currentDromoState( orbital_element_conversions::dromoElement5Index );
    double zeta6 = currentDromoState( orbital_element_conversions::dromoElement6Index );
    double zeta7 = currentDromoState( orbital_element_conversions::dromoElement7Index );

    Eigen::Vector3d perturbingAccelerationInRswFrame = perturbingNonPotentialAccelerationsInRswFrame + perturbingPotentialAccelerationsInRswFrame;
    double term1 = perturbingAccelerationInRswFrame( 0 ) / ( zeta3 * s - 2 * perturbingPotential );
    double term2 = s / zeta3 - 1;
    double term3 = perturbingAccelerationInRswFrame( 1 ) / ( zeta3 * s * transverseVelocity );

    stateDerivative( orbital_element_conversions::dromoElement1Index ) = sinPhi / s * term1 - term2 * cosPhi * zeta3Derivative;
    stateDerivative( orbital_element_conversions::dromoElement2Index ) = - cosPhi / s * term1 - term2 * sinPhi * zeta3Derivative;
    stateDerivative( orbital_element_conversions::dromoElement4Index ) =
            1 / ( 2 * s ) * ( term3 * ( zeta7 * cosDeltaPhi - zeta6 * sinDeltaPhi ) + zeta5 * ( transverseVelocity - s ) );
    stateDerivative( orbital_element_conversions::dromoElement5Index ) =
            1 / ( 2 * s ) * ( term3 * ( zeta6 * cosDeltaPhi + zeta7 * sinDeltaPhi ) - zeta4 * ( transverseVelocity - s ) );
    stateDerivative( orbital_element_conversions::dromoElement4Index ) =
            - 1 / ( 2 * s ) * ( term3 * ( zeta5 * cosDeltaPhi - zeta4 * sinDeltaPhi ) - zeta7 * ( transverseVelocity - s ) );
    stateDerivative( orbital_element_conversions::dromoElement4Index ) =
            - 1 / ( 2 * s ) * ( term3 * ( zeta4 * cosDeltaPhi + zeta5 * sinDeltaPhi ) + zeta6 * ( transverseVelocity - s ) );

    // Compute time derivative
    if ( timeType == scaled_physical_time )
    {
        stateDerivative( orbital_element_conversions::dromoTimeIndex ) = 1 / ( zeta3 * std::pow( s, 2 ) );
    }
    else if ( timeType == linear_time_element or timeType == constant_time_element )
    {
        double semiMajorAxis = - 1 / ( 2 * energy );
        double k4 = zeta3 + std::sqrt( - 2 * energy );
        double k3 = currentDromoState( orbital_element_conversions::dromoElement1Index ) * cosPhi +
                currentDromoState( orbital_element_conversions::dromoElement2Index ) * sinPhi;
        double k1 = std::sqrt( semiMajorAxis ) * radialVelocity / std::pow( s, 2 ) * ( ( zeta3 + s ) / k4 + 2 * k3 / zeta3 + 1 );
        double k2 = 1 / std::pow( s, 2 ) * ( k4 / zeta3 + k3 / k4 + std::pow( radialVelocity, 2 ) / ( k4 * s ) );

        double termT1 = 6 * semiMajorAxis * std::pow( radialVelocity, k4 + k3 );
        double termT2 = ( perturbingAccelerationInRswFrame( 0 ) / ( zeta3 * s ) - 2 * perturbingPotential ) * k2;

        if ( timeType == linear_time_element )
        {
            stateDerivative( orbital_element_conversions::dromoTimeIndex ) = std::pow( semiMajorAxis, 3/2 ) *
                    ( 1 + energyDerivative * ( termT1 + k1 ) + termT2 );
        }
        else // Constant time element
        {
            stateDerivative( orbital_element_conversions::dromoTimeIndex ) = std::pow( semiMajorAxis, 3/2 ) *
                    ( energyDerivative * ( termT1 - 3 * semiMajorAxis * currentIndependentVariable + k1 ) ) + termT2;
        }
    }
    else
    {
        throw std::runtime_error( "Error when computing Dromo(P) state derivative: selected time element (" +  std::to_string( timeType ) +
                     + ") is not implemented." );
    }

    return stateDerivative;
}

//! Class for computing the state derivative of translational motion of N bodies, using Dromo propagator.
template< typename StateScalarType = double, typename TimeType = double >
class NBodyDromoPStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType > {
public:

    NBodyDromoPStateDerivative (
            const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
            const std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
            const std::vector< std::string >& bodiesToIntegrate,
            const TimeElementType timeType,
            const std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& initialCartesianElements,
            const double initialPhysicalTime,
            const bool useEnergyElement) :
            NBodyStateDerivative< StateScalarType, TimeType >(
                    accelerationModelsPerBody, centralBodyData, dromo_p, bodiesToIntegrate ),
            timeType_( timeType ),
            initialPhysicalTime_( initialPhysicalTime ),
            useEnergyElement_( useEnergyElement ) {
        // Check if specified timeType is valid
        if ( not( timeType_ == scaled_physical_time or timeType_ == linear_time_element or
                  timeType_ == constant_time_element ) ) {
            throw std::runtime_error( "Error when setting N Body Dromo(P) propagator, specified time type (" +
                                      std::to_string( timeType_ ) + ") is not implemented." );
        }

        // Check if a single body is to be integrated
        if ( bodiesToIntegrate.size( ) != 1 ) {
            throw std::runtime_error(
                    "Error when setting N Body Dromo(P) propagator, only possible to integrate a single body ("
                    + std::to_string( bodiesToIntegrate.size( ) ) + " bodies specified)." );
        }

        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameterFunction_ = removeCentralGravityAccelerations(
                centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                this->accelerationModelsPerBody_, this->removedCentralAccelerations_ ).at( 0 );
        this->createAccelerationModelList( );

        // Select unit of length
        unitOfLength_ = orbital_element_conversions::computeDromoUnitOfLength(
                initialCartesianElements.block( 0, 0, 3, 1 ) );

        // Select initial independent variable
        Eigen::Vector6d keplerElements = convertCartesianToKeplerianElements( initialCartesianElements, centralBodyGravitationalParameterFunction_() );
        initialIndependentVariable_ = keplerElements( orbital_element_conversions::trueAnomalyIndex );
    }

    void convertToOutputSolution (
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution,
            const TimeType& currentIndependentVariable,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSolution)
    {
        currentCartesianLocalSolution.block( 0, 0, 6, 1 ) = orbital_element_conversions::convertDromoToCartesianElements(
                internalSolution, centralBodyGravitationalParameterFunction_( ), initialIndependentVariable_,
                currentIndependentVariable, unitOfLength_, fixedPerturbingPotential_, useEnergyElement_ );
    }

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution (
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& physicalTime)
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState = orbital_element_conversions::convertCartesianToDromoElements(
                cartesianSolution.block( 0, 0, 6, 1 ), centralBodyGravitationalParameterFunction_(), unitOfLength_,
                fixedPerturbingPotential_, physicalTime - initialPhysicalTime_ , useEnergyElement_, timeType_);
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
        double timeFromPropagationStart = orbital_element_conversions::convertDromoTimeToPhysicalTime(
                stateOfSystemToBeIntegrated, centralBodyGravitationalParameterFunction_(), independentVariable,
                unitOfLength_, useEnergyElement_, timeType_);

        return timeFromPropagationStart + initialPhysicalTime_;
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
        return initialIndependentVariable_;
    }

private:

    //! Gravitational parameter of the central body used to convert Cartesian to Keplerian orbits, and vice versa
    std::function< double( ) > centralBodyGravitationalParameterFunction_;

    //! Type of time that is to be integrated
    TimeElementType timeType_;

    //! Initial value of the independent variable
    double initialIndependentVariable_;

    //! Initial value of the physical time
    double initialPhysicalTime_;

    //! Flag indicating whether energy element is or not being used
    bool useEnergyElement_;

    //! Unit of length used to make coordinates dimensionless
    double unitOfLength_;

    //!
    double fixedPerturbingPotential_ = 0.0;
};

extern template class NBodyDromoPStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class NBodyDromoPStateDerivative< long double, double >;
extern template class NBodyDromoPStateDerivative< double, Time >;
extern template class NBodyDromoPStateDerivative< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat


#endif //TUDATBUNDLE_NBODYDROMOSTATEDERIVATIVE_H
