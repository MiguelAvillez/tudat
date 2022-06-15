#include "tudat/astro/propagators/nBodyStabilizedCowellStateDerivative.h"

namespace tudat
{

namespace propagators
{

double computePhysicalTimeDerivativeForStabilizedCowell(
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant )
{
    const Eigen::Vector3d position( currentStabilizedCowellState.segment( 0, 3 ) );
    return sundmanConstant * position.squaredNorm();
}

double computeLinearTimeElementDerivativeForStabilizedCowell(
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant,
        const double energyDerivative,
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& nonConservativeAccelerationsInInertialFrame)
{
    const double energy = currentStabilizedCowellState( orbital_element_conversions::stabilizedCowellEnergy );
    const Eigen::Vector3d position( currentStabilizedCowellState.segment( 0, 3 ) );
    const double positionNorm = position.norm();
    const Eigen::Vector3d velocity( currentStabilizedCowellState.segment( 3, 3 ) );

    double linearTimeElementDerivative = sundmanConstant / ( 2.0 * energy) * ( centralBodyGravitationalParameter + positionNorm *
            position.dot(nonConservativeAccelerationsInInertialFrame - velocity * energyDerivative /
            ( std::pow( sundmanConstant, 2) * std::pow( positionNorm, 2) * energy ) ) );

    return linearTimeElementDerivative;
}

double computeLinearTimeElementToPhysicalTimeBias (
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant)
{
    const double energy = currentStabilizedCowellState( orbital_element_conversions::stabilizedCowellEnergy );
    const Eigen::Vector3d position( currentStabilizedCowellState.segment( 0, 3 ) );
    const Eigen::Vector3d velocity( currentStabilizedCowellState.segment( 3, 3 ) );

    return position.dot( velocity ) / ( 2.0 * energy * sundmanConstant * position.norm() );
}

double convertLinearTimeElementToPhysicalTime (
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant)
{
    const double linearTimeElement = currentStabilizedCowellState( orbital_element_conversions::stabilizedCowellTime );

    double physicalTime = linearTimeElement -
            computeLinearTimeElementToPhysicalTimeBias(currentStabilizedCowellState, sundmanConstant);

    return physicalTime;
}

double convertPhysicalTimeToLinearTimeElement (
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant)
{
    const double physicalTime = currentStabilizedCowellState( orbital_element_conversions::stabilizedCowellTime );

    double linearTimeElement = physicalTime +
            computeLinearTimeElementToPhysicalTimeBias(currentStabilizedCowellState, sundmanConstant);

    return linearTimeElement;
}

Eigen::Vector7d computeStateDerivativeExceptTimeForStabilizedCowell(
        const Eigen::Vector8d& currentStabilizedCowellState,
        const double sundmanConstant,
        const double conservativeAccelerationsPotential,
        const Eigen::Vector3d& nonConservativeAccelerationsInInertialFrame,
        const Eigen::Vector3d& conservativeAccelerationsInInertialFrame)
{
    const Eigen::Vector3d position( currentStabilizedCowellState.segment( 0, 3 ) );
    const double positionNorm = position.norm();
    const Eigen::Vector3d velocity( currentStabilizedCowellState.segment( 3, 3 ) );
    const double velocityNorm = velocity.norm();
    const double energy = currentStabilizedCowellState( orbital_element_conversions::stabilizedCowellEnergy );

    Eigen::Vector7d stateDerivative;

    stateDerivative( orbital_element_conversions::stabilizedCowellEnergy ) =
            - ( nonConservativeAccelerationsInInertialFrame - conservativeAccelerationsInInertialFrame).dot( velocity );

    stateDerivative( orbital_element_conversions::stabilizedCowellXCartesianPositionIndex ) =
            currentStabilizedCowellState( orbital_element_conversions::stabilizedCowellXCartesianVelocityIndex );
    stateDerivative( orbital_element_conversions::stabilizedCowellYCartesianPositionIndex ) =
            currentStabilizedCowellState( orbital_element_conversions::stabilizedCowellYCartesianVelocityIndex );
    stateDerivative( orbital_element_conversions::stabilizedCowellZCartesianPositionIndex ) =
            currentStabilizedCowellState( orbital_element_conversions::stabilizedCowellZCartesianVelocityIndex );

    Eigen::Vector3d acceleration = velocity * position.dot(velocity) / std::pow( positionNorm, 2) -
            std::pow( sundmanConstant, 2) *
            (
                ( energy + 0.5 * std::pow( velocityNorm / positionNorm, 2) + conservativeAccelerationsPotential ) * position +
                std::pow( positionNorm, 2) * ( nonConservativeAccelerationsInInertialFrame + conservativeAccelerationsInInertialFrame )
            );
    stateDerivative( orbital_element_conversions::stabilizedCowellXCartesianVelocityIndex ) = acceleration(0);
    stateDerivative( orbital_element_conversions::stabilizedCowellYCartesianVelocityIndex ) = acceleration(1);
    stateDerivative( orbital_element_conversions::stabilizedCowellZCartesianVelocityIndex ) = acceleration(2);

    return stateDerivative;
}

double computeEnergy(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector6d& currentCartesianState,
        const double conservativeAccelerationsPotential)
{
    const Eigen::Vector3d position( currentCartesianState.segment( 0, 3 ) );
    const Eigen::Vector3d velocity( currentCartesianState.segment( 3, 3 ) );

    return centralBodyGravitationalParameter / position.norm() - std::pow( velocity.norm(), 2 ) / 2.0 - conservativeAccelerationsPotential;
}

template class NBodyStabilizedCowellStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class NBodyStabilizedCowellStateDerivative< long double, double >;
template class NBodyStabilizedCowellStateDerivative< double, Time >;
template class NBodyStabilizedCowellStateDerivative< long double, Time >;
#endif


} // namespace propagators

} // namespace tudat
