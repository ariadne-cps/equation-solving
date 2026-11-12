/***************************************************************************
 *            equation-solving/solver_attributes.hpp
 *
 *  Copyright  2011-26  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Equation Solving.
 *
 *  Ariadne Equation Solving is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Equation Solving is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Equation Solving. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file equation-solving/solver_attributes.hpp
 *  @brief Typed configuration attributes for equation solvers.
 */

#ifndef ARIADNE_SOLVER_ATTRIBUTES_HPP
#define ARIADNE_SOLVER_ATTRIBUTES_HPP

#include "numeric/builtin.hpp"
#include "utility/attribute.hpp"

namespace Ariadne {

/** @brief Specify a maximum permissible equation residual. */
struct MaximumError : Attribute<ApproximateDouble> {
    /** @brief Construct a maximum error attribute from an approximate double. */
    MaximumError(ApproximateDouble x) : Attribute<ApproximateDouble>(x) {}
    /** @brief Construct a maximum error attribute from a native double. */
    MaximumError(double x) : Attribute<ApproximateDouble>(x) {}
};
/** @brief Configure the numerical sweep threshold. */
struct SweepThreshold : Attribute<ApproximateDouble> {
    using Attribute<ApproximateDouble>::Attribute;
};
/** @brief Limit the number of iterative solver steps. */
struct MaximumNumericTypeOfSteps : Attribute<Nat> {
    using Attribute<Nat>::Attribute;
};

/** @brief Construct the maximum error attribute using the generator. */
static const Generator<MaximumError> maximum_error = Generator<MaximumError>();
/** @brief Construct the sweep threshold attribute using the generator. */
static const Generator<SweepThreshold> sweep_threshold = Generator<SweepThreshold>();
/** @brief Construct the maximum step count attribute using the generator. */
static const Generator<MaximumNumericTypeOfSteps> maximum_number_of_steps = Generator<MaximumNumericTypeOfSteps>();

} // namespace Ariadne

#endif /* ARIADNE_SOLVER_ATTRIBUTES_HPP */
