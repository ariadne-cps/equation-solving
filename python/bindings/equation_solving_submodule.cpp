/***************************************************************************
 *            solver_submodule.cpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11.hpp"
#include "function-utilities.hpp"

#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "equation-solving/equation_solver_interface.hpp"
#include "equation-solving/equation_solver.hpp"
#include "equation-solving/interval_newton_method.hpp"
#include "equation-solving/krawczyk_method.hpp"
#include "function/taylor_function.hpp"


using namespace Ariadne;

namespace Ariadne {

template<class T> Nat __hash__(const T&);
template<> Nat __hash__<FloatDPBoundsVector>(const FloatDPBoundsVector&) {
    return 0;
}

class EquationSolverWrapper
  : public pybind11::wrapper< EquationSolverInterface >
{
    typedef EquationSolverInterface::ValidatedNumericType ValidatedNumericType;
    typedef EquationSolverInterface::ApproximateNumericType ApproximateNumericType;
  public:
    EquationSolverInterface* clone() const { return this->get_override("clone")(); }
    Void set_maximum_error(ApproximateDouble me) { this->get_override("set_maximum_error")(me); }
    ExactDouble maximum_error() const { return this->get_override("maximum_error")(); }
    Void set_maximum_number_of_steps(Nat ns) { this->get_override("set_maximum_number_of_steps")(ns); }
    Nat maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    Vector<ValidatedNumericType> zero(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("zero")(f,bx); }
    Vector<ValidatedNumericType> fixed_point(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("fixed_point")(f,bx); }
    Vector<ValidatedNumericType> solve(const ValidatedVectorMultivariateFunction& f, const Vector<ValidatedNumericType>& pt) const {
        return this->get_override("solve")(f,pt); }
    Vector<ValidatedNumericType> solve(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        return this->get_override("solve")(f,bx); }
    ValidatedVectorMultivariateFunctionPatch implicit(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& pd, const ExactBoxType& bx) const {
        return this->get_override("implicit")(f,pd,bx); }
    ValidatedScalarMultivariateFunctionPatch implicit(const ValidatedScalarMultivariateFunction& f, const ExactBoxType& pd, const ExactIntervalType& ivl) const {
        return this->get_override("implicit")(f,pd,ivl); }
    ValidatedVectorMultivariateFunctionPatch continuation(const ValidatedVectorMultivariateFunction& f, const Vector< ApproximateNumericType>& a, const ExactBoxType& X,  const ExactBoxType& A) const {
        return this->get_override("continuation")(f,a,X,A); }
    List< Vector<ValidatedNumericType> > list_solve_all(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx) const {
        Set< Vector<ValidatedNumericType> > set = this->get_override("solve_all")(f,bx);
        return make_list(set);
    }
    Void _write(OutputStream& os) const { this->get_override("_write")(os); }
};




} // namespace Ariadne


Void export_equation_solvers(pybind11::module& module)
{
    typedef EquationSolverInterface::ValidatedNumericType ValidatedNumericType;
    pybind11::class_<EquationSolverInterface,EquationSolverWrapper> equation_solver_interface_class(module,"EquationSolverInterface");
    equation_solver_interface_class.def("maximum_error", &EquationSolverInterface::maximum_error);
    equation_solver_interface_class.def("set_maximum_error", &EquationSolverInterface::set_maximum_error);
    equation_solver_interface_class.def("maximum_number_of_steps", &EquationSolverInterface::maximum_number_of_steps);
    equation_solver_interface_class.def("set_maximum_number_of_steps", &EquationSolverInterface::set_maximum_number_of_steps);
    equation_solver_interface_class.def("zero", &EquationSolverInterface::zero);
    equation_solver_interface_class.def("fixed_point", &EquationSolverInterface::fixed_point);
    equation_solver_interface_class.def("solve", (Vector<ValidatedNumericType>(EquationSolverInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&)const) &EquationSolverInterface::solve);
    equation_solver_interface_class.def("implicit",(ValidatedVectorMultivariateFunctionPatch(EquationSolverInterface::*)(const ValidatedVectorMultivariateFunction&,const ExactBoxType&,const ExactBoxType&)const) &EquationSolverInterface::implicit);
    equation_solver_interface_class.def("implicit",(ValidatedScalarMultivariateFunctionPatch(EquationSolverInterface::*)(const ValidatedScalarMultivariateFunction&,const ExactBoxType&,const ExactIntervalType&)const) &EquationSolverInterface::implicit);
    equation_solver_interface_class.def("solve_all",[](EquationSolverInterface const& s, const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bx){ return static_cast< List<Vector<ValidatedNumericType>> >( s.solve_all(f,bx) );});
    equation_solver_interface_class.def("__str__",&__cstr__<EquationSolverInterface>);

    pybind11::class_<IntervalNewtonSolver, EquationSolverInterface> interval_newton_solver_class(module,"IntervalNewtonSolver");
    interval_newton_solver_class.def(pybind11::init<ApproximateDouble,Nat>());

    pybind11::class_<KrawczykSolver, EquationSolverInterface> krawczyk_solver_class(module,"KrawczykSolver");
    krawczyk_solver_class.def(pybind11::init<ApproximateDouble,Nat>());

    pybind11::class_<FactoredKrawczykSolver, EquationSolverInterface> factored_krawczyk_solver_class(module,"FactoredKrawczykSolver");
    factored_krawczyk_solver_class.def(pybind11::init<ApproximateDouble,Nat>());

}





