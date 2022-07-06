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

#include "tudat/astro/propagators/nBodyDromoStateDerivative.h"

namespace tudat
{

namespace propagators
{


template class NBodyDromoPStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class NBodyDromoPStateDerivative< long double, double >;
template class NBodyDromoPStateDerivative< double, Time >;
template class NBodyDromoPStateDerivative< long double, Time >;
#endif


} // namespace propagators

} // namespace tudat