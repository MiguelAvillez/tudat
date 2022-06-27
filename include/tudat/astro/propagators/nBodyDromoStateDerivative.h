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

namespace tudat
{

namespace propagators
{

//! Class for computing the state derivative of translational motion of N bodies, using Dromo propagator.
template< typename StateScalarType = double, typename TimeType = double >
class NBodyDromoStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

};

extern template class NBodyDromoStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class NBodyDromoStateDerivative< long double, double >;
extern template class NBodyDromoStateDerivative< double, Time >;
extern template class NBodyDromoStateDerivative< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat


#endif //TUDATBUNDLE_NBODYDROMOSTATEDERIVATIVE_H
