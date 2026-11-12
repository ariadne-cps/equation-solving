/***************************************************************************
 *            tests/test_equation_solvers.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file tests/test_equation_solvers.cpp
 *  @brief Test root finding, contraction and implicit equation solving.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "equation-solving/equation_solver.hpp"
#include "equation-solving/interval_newton_method.hpp"
#include "equation-solving/krawczyk_method.hpp"
#include "function/function.hpp"
#include "function/function_patch.hpp"
#include "function/taylor_function.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "algebra/differential.hpp"
#include "function/formula.hpp"
#include "logging/logging.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace std;

typedef Vector<DyadicInterval> DyadicIntervalVector;

class TestEquationSolver : public Ariadne::Testing::ContextFixture
{
  private:
    std::unique_ptr<EquationSolverInterface> solver;
    StringType solver_class_name;
  public:
    TestEquationSolver(Ariadne::Testing::Context& context,
                       const EquationSolverInterface& s, const char* n)
        : ContextFixture(context), solver(s.clone()), solver_class_name(n) { }

    Void test_solve() {
        ARIADNE_TESTING_NOTE("*solver = " << *solver);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(1,0);
        ExactBoxType d({{0.0_x,1.0_x}});
        EffectiveVectorMultivariateFunction f({(x*x+1)*x-1});
        FloatDPBoundsVector p=solver->solve(f,d);
        ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(contains,ExactIntervalType(0.6823_pr,0.6824_pr),p[0]);
    }

    Void test_implicit() {
        //TaylorModelAccuracy::set_default_sweep_threshold(1e-12);

        // FIXME: Should be able to use numbers yielding exact results for p,r
        EffectiveScalarMultivariateFunction aa=EffectiveScalarMultivariateFunction::coordinate(1,0);
        EffectiveScalarMultivariateFunction a=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction bb;
        ExactBoxType p,r;
        EffectiveVectorMultivariateFunction f;
        ValidatedVectorMultivariateFunctionPatch h;
        EffectiveVectorMultivariateFunction e;

        // Test solution of x-a=0. This should be very easy to solve.
        p=ExactBoxType({{-0.25_x,0.25_x}});
        r=ExactBoxType({{-2.0_x,2.0_x}});
        f=EffectiveVectorMultivariateFunction({x-a});
        ARIADNE_TESTING_NOTE("f = " << (f));
        h=solver->implicit(f,p,r);
        ARIADNE_TESTING_NOTE("h = " << (h));
        e=EffectiveVectorMultivariateFunction(1u,aa);
        ARIADNE_TESTING_NOTE("e = " << (e));
        ARIADNE_TESTING_EXPECT_LESS(norm(h-e), 1e-8_pr);

        // Test solution of 4x^2+x-4-a=0 on [0.875_x,1.125_x]. There is a unique solution with positive derivative.
        p=ExactBoxType({ExactIntervalType(0.875_x,1.125_x)});
        r=ExactBoxType({ExactIntervalType(0.25_x,1.25_x)});
        f=EffectiveVectorMultivariateFunction({(x*x+1)*x-a});
        ARIADNE_TESTING_NOTE("f = " << (f));
        h=solver->implicit(f,p,r);
        ARIADNE_TESTING_NOTE("h = " << (h));
        bb=EffectiveScalarMultivariateFunction(aa-DyadicInterval(p[0]).midpoint())/DyadicInterval(p[0]).radius();
        Decimal a0(0.682328_dec), a1(0.0521547_dec), a2(-0.0023232_dec), a3(0.000147778_dec);
        e=EffectiveVectorMultivariateFunction( { a0+bb*(a1+bb*(a2+bb*a3)) } );
        ARIADNE_TESTING_NOTE("e = " << (e));
        ARIADNE_TESTING_EXPECT_LESS(norm(h-e), 1e-4_pr);

        // Test solution of 4x^2+x-4-a=0 on [-0.25_x,0.25_x]. There is a unique solution with positive derivative.
        p=ExactBoxType({ExactIntervalType(-0.25_x,0.25_x)});
        r=ExactBoxType({ExactIntervalType(0.25_x,2.0_x)});
        f=EffectiveVectorMultivariateFunction({4*x+x*x-a-4});
        ARIADNE_TESTING_NOTE("f = " << (f));
        h=solver->implicit(f,p,r);
        ARIADNE_TESTING_NOTE("h = " << (h));
        bb=EffectiveScalarMultivariateFunction(aa-DyadicInterval(p[0]).midpoint())/DyadicInterval(p[0]).radius();
        Decimal c0(0.828427_dec), c1(0.0441942_dec), c2(-0.000345267_dec), c3(0.00000539468_dec);
        e=EffectiveVectorMultivariateFunction( { c0+bb*(c1+bb*(c2+bb*c3)) } );
        ARIADNE_TESTING_NOTE("e = " << (e));
        ARIADNE_TESTING_EXPECT_LESS(norm(h-e), 1e-4_pr);

        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1].
        // There is at most one solution, but this lies partially outside the range.
        // Should obtain PartialSolutionException
        p=ExactBoxType({ExactIntervalType(-1,1)});
        r=ExactBoxType({ExactIntervalType(-1,1)});
        f=EffectiveVectorMultivariateFunction({x-2*a});
        ARIADNE_TESTING_NOTE("f = " << (f));
        try {
            h=solver->implicit(f,p,r);
            ARIADNE_TESTING_NOTE(solver_class_name<<" silently returns partially defined vector implicit function.");
        }
        catch(const EquationSolverException&) {
            ARIADNE_TESTING_EXPECT_THROWS(EquationSolverException, solver->implicit(f,p,r));
            ARIADNE_TESTING_NOTE(solver_class_name<<" throws error on partially defined vector implicit function.");
        }

    }

    Void test_scalar_implicit() {
        //TaylorModelAccuracy::set_default_sweep_threshold(1e-12);

        EffectiveScalarMultivariateFunction aa=EffectiveScalarMultivariateFunction::coordinate(1,0);
        EffectiveScalarMultivariateFunction a=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,1);
        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1]. There is at most one solution.
        // Uses scalar implicit
        ExactBoxType p; ExactIntervalType r;
        EffectiveScalarMultivariateFunction e,f,s; // s is unscaling functions
        ValidatedScalarMultivariateFunctionPatch h;

        ARIADNE_TESTING_NOTE("*solver = " << (*solver));

        // Test solution of 4x^2+x-4-a=0 on [0.875_x,1.125_x]. There is a unique solution with positive derivative.
        p=ExactBoxType({ExactIntervalType(0.875_x,1.125_x)});
        r=ExactIntervalType(0.25_x,1.25_x);
        f=EffectiveScalarMultivariateFunction((x*x+1)*x-a);
        ARIADNE_TESTING_NOTE("f = " << (f));
        h=solver->implicit(f,p,r);
        ARIADNE_TESTING_NOTE("h = " << (h));
        s=EffectiveScalarMultivariateFunction(aa-DyadicInterval(p[0]).midpoint())/DyadicInterval(p[0]).radius();
        Decimal a0(0.682328_dec), a1(0.0521547_dec), a2(-0.0023232_dec), a3(0.000147778_dec);
        e=EffectiveScalarMultivariateFunction( a0+s*(a1+s*(a2+s*a3)) );
        ARIADNE_TESTING_NOTE("e = " << (e));
        ARIADNE_TESTING_EXPECT_LESS(mag((h-e).range()), 1e-4_pr);

        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1].
        // There is at most one solution, but this lies partially outside the range.
        // Should obtain PartialSolutionException
        p=ExactBoxType({ExactIntervalType(-1,1)});
        r=ExactIntervalType(-1,1);
        f=EffectiveScalarMultivariateFunction(x-2*a);
        ARIADNE_TESTING_NOTE("f = " << (f));
        ValidatedScalarMultivariateFunctionPatch g=ValidatedScalarMultivariateTaylorFunctionModelDP(product(p,r),f,ThresholdSweeper<FloatDP>(dp,1e-12));
        ARIADNE_TESTING_NOTE("g = " << (g));
        try {
            h=solver->implicit(g,p,r);
            ARIADNE_TESTING_NOTE(solver_class_name<<" silently returns partially defined scalar implicit function.");
        }
        catch(const EquationSolverException&) {
            ARIADNE_TESTING_EXPECT_THROWS(EquationSolverException, solver->implicit(g,p,r));
            ARIADNE_TESTING_NOTE(solver_class_name<<" throws error on partially defined scalar implicit function.");
        }
        catch(const DomainException&) {
            ARIADNE_TESTING_EXPECT_THROWS(DomainException, solver->implicit(g,p,r));
            ARIADNE_TESTING_NOTE(solver_class_name<<" throws DomainException on partially defined scalar implicit function.");
        }
    }

};


class TestEquationSolverBase : public EquationSolverBase
{
  public:
    enum class StepMode { Identity, Disjoint, Wide, HighError };

    TestEquationSolverBase(StepMode mode, ApproximateDouble max_error=1e-5, Nat max_steps=0)
        : EquationSolverBase(max_error,max_steps), _mode(mode) { }

    virtual TestEquationSolverBase* clone() const { return new TestEquationSolverBase(*this); }
    virtual Void _write(OutputStream& os) const { os << "TestEquationSolverBase"; }

  protected:
    virtual Vector<ValidatedNumericType>
    step(const ValidatedVectorMultivariateFunction&,
         const Vector<ValidatedNumericType>& x) const
    {
        if(_mode==StepMode::Identity || _mode==StepMode::HighError) { return x; }
        Vector<ValidatedNumericType> r(x.size(),dp);
        for(SizeType i=0; i!=x.size(); ++i) {
            if(_mode==StepMode::Disjoint) {
                r[i]=ValidatedNumericType(10,11,dp);
            } else {
                r[i]=ValidatedNumericType(-10,10,dp);
            }
        }
        return r;
    }

    virtual ValidatedVectorMultivariateFunctionPatch
    implicit_step(const ValidatedVectorMultivariateFunction&,
                  const ValidatedVectorMultivariateFunctionPatch&,
                  const ValidatedVectorMultivariateFunctionPatch& x) const
    {
        if(_mode==StepMode::HighError) {
            ValidatedScalarMultivariateTaylorFunctionModelDP high_error=
                ValidatedScalarMultivariateTaylorFunctionModelDP::zero(
                    x.domain(),ThresholdSweeper<FloatDP>(dp,1e-12));
            high_error.set_error(FloatDPError(2u,dp));
            return ValidatedVectorMultivariateFunctionPatch(
                ValidatedVectorMultivariateTaylorFunctionModelDP(x.result_size(),high_error));
        }
        return x;
    }

  private:
    StepMode _mode;
};

class TestEquationSolverContracts : public Ariadne::Testing::ContextFixture
{
  public:
    using ContextFixture::ContextFixture;

    Void test_equation_solver_base_paths()
{
    EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(1,0);
    EffectiveVectorMultivariateFunction root({x});
    EffectiveVectorMultivariateFunction no_root({x+2});
    ExactBoxType unit({{-1.0_x,1.0_x}});

    TestEquationSolverBase identity(TestEquationSolverBase::StepMode::Identity);
    TestEquationSolverBase disjoint(TestEquationSolverBase::StepMode::Disjoint,1e-5,1);
    TestEquationSolverBase wide(TestEquationSolverBase::StepMode::Wide);
    TestEquationSolverBase wide_loose(TestEquationSolverBase::StepMode::Wide,3.0,0);
    TestEquationSolverBase high_error(TestEquationSolverBase::StepMode::HighError,1e-8,1);

    identity.set_function_factory(identity.function_factory());

    FloatDPBoundsVector interval_point({FloatDPBounds(-1,1,dp)});
    ARIADNE_TESTING_NOTE("identity.solve(root,interval_point) = " << (identity.solve(root,interval_point)));

    ARIADNE_TESTING_EXPECT_THROWS(NoSolutionException, identity.zero(no_root,unit));
    ARIADNE_TESTING_NOTE("identity.zero(root,unit) = " << (identity.zero(root,unit)));
    ARIADNE_TESTING_EXPECT_THROWS(EquationSolverException, wide.zero(root,unit));
    ARIADNE_TESTING_EXPECT_THROWS(NoSolutionException, disjoint.zero(root,unit));
    ARIADNE_TESTING_EXPECT_EQUAL(disjoint.solve_all(root,unit).size(),0u);

    ARIADNE_TESTING_EXPECT_EQUAL(wide_loose.solve_all(root,unit).size(),0u);

    EffectiveVectorMultivariateFunction fixed_map({(x+1)/2});
    ExactBoxType fixed_domain({{0.0_x,2.0_x}});
    IntervalNewtonSolver newton(1e-8,12);
    FloatDPBoundsVector fp=newton.fixed_point(fixed_map,fixed_domain);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(contains,ExactIntervalType(1.0_x,1.0_x),fp[0]);

    EffectiveVectorMultivariateFunction impossible({x*x+1});
    ARIADNE_TESTING_EXPECT_EQUAL(newton.solve_all(impossible,unit).size(),0u);
    ARIADNE_TESTING_EXPECT_THROWS(NoSolutionException, newton.solve(impossible,unit));

    IntervalNewtonSolver split_newton(1e-3,12);
    EffectiveVectorMultivariateFunction two_roots({x*x-1});
    ExactBoxType two_root_domain({{-2.0_x,2.0_x}});
    auto roots=split_newton.solve_all(two_roots,two_root_domain);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(roots.size(), 2u);
    ARIADNE_TESTING_EXPECT_THROWS(EquationSolverException, split_newton.solve(two_roots,two_root_domain));

    IntervalNewtonSolver loose_newton(3.0,4);
    EffectiveVectorMultivariateFunction singular({x*x});
    ExactBoxType small({{-0.125_x,0.125_x}});
    ARIADNE_TESTING_EXPECT_EQUAL(loose_newton.solve_all(singular,small).size(),0u);

    FloatDPApproximationVector parameter(1u,dp);
    parameter[0]=FloatDPApproximation(0.0_x,dp);
    ARIADNE_TESTING_EXPECT_EXCEPTION(identity.continuation(root,parameter,unit,unit));

    EffectiveScalarMultivariateFunction a=EffectiveScalarMultivariateFunction::coordinate(2,0);
    EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
    EffectiveVectorMultivariateFunction implicit_function({y-a});
    ExactBoxType par({{-1.0_x,1.0_x}});
    ExactBoxType implicit_range({{-1.0_x,1.0_x}});

    IntervalNewtonSolver zero_step_newton(1e-8,0);
    ARIADNE_TESTING_EXPECT_THROWS(NoSolutionException, zero_step_newton.implicit(implicit_function,par,implicit_range));

    IntervalNewtonSolver one_step_newton(1e-20,1);
    ARIADNE_TESTING_NOTE("one_step_newton.implicit(implicit_function,par,implicit_range) = " << (one_step_newton.implicit(implicit_function,par,implicit_range)));

    EffectiveVectorMultivariateFunction disjoint_implicit_function({y-10});
    ARIADNE_TESTING_EXPECT_THROWS(NoSolutionException, newton.implicit(disjoint_implicit_function,par,implicit_range));
    ARIADNE_TESTING_EXPECT_THROWS(UnknownSolutionException, high_error.implicit(implicit_function,par,implicit_range));

    auto& logging_configuration=Logger::instance().configuration();
    unsigned int const old_verbosity=logging_configuration.verbosity();
    logging_configuration.set_verbosity(10u);
    TestEquationSolverBase logging_identity(TestEquationSolverBase::StepMode::Identity,1e-5,1);
    ARIADNE_TESTING_NOTE("logging_identity.solve_all(root,unit) = " << (logging_identity.solve_all(root,unit)));
    ARIADNE_TESTING_NOTE("logging_identity.zero(root,unit) = " << (logging_identity.zero(root,unit)));
    ARIADNE_TESTING_NOTE("one_step_newton.implicit(implicit_function,par,implicit_range) = " << (one_step_newton.implicit(implicit_function,par,implicit_range)));
    EffectiveScalarMultivariateFunction scalar_implicit_function(y-a);
    ARIADNE_TESTING_NOTE("one_step_newton.implicit(scalar_implicit_function,par,implicit_range[0]) = " << (one_step_newton.implicit(scalar_implicit_function,par,implicit_range[0])));
    logging_configuration.set_verbosity(old_verbosity);

    ARIADNE_TESTING_EXPECT_EXCEPTION(newton.implicit(root,par,implicit_range));

    ExactBoxType bad_range;
    ARIADNE_TESTING_EXPECT_EXCEPTION(newton.implicit(implicit_function,par,bad_range));
}


Void test_interval_newton_paths()
{
    auto vars=EffectiveScalarMultivariateFunction::coordinates(3);
    EffectiveScalarMultivariateFunction a=vars[0];
    EffectiveScalarMultivariateFunction u=vars[1];
    EffectiveScalarMultivariateFunction v=vars[2];
    ExactBoxType parameter_domain({{-0.25_x,0.25_x}});
    ExactBoxType solution_domain({{-1.0_x,1.0_x},{-1.0_x,1.0_x}});

    IntervalNewtonSolver solver(1e-8,12);
    EffectiveVectorMultivariateFunction coupled({u-a,v+a});
    auto solutions=solver.implicit(coupled,parameter_domain,solution_domain);
    EffectiveScalarMultivariateFunction p=EffectiveScalarMultivariateFunction::coordinate(1,0);
    EffectiveVectorMultivariateFunction expected({p,-p});
    ARIADNE_TESTING_EXPECT_EQUAL(solutions.result_size(),2u);
    ARIADNE_TESTING_EXPECT_LESS(norm(solutions-expected), 1e-4_pr);

    auto scalar_vars=EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction singular({scalar_vars[1]*scalar_vars[1]-scalar_vars[0]});
    ExactBoxType scalar_domain({{-1.0_x,1.0_x}});
    ARIADNE_TESTING_EXPECT_THROWS(SingularJacobianException, solver.implicit(singular,parameter_domain,scalar_domain));

    auto& config=Logger::instance().configuration();
    unsigned int const old_verbosity=config.verbosity();
    config.set_verbosity(10u);
    EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(1,0);
    EffectiveVectorMultivariateFunction linear({x-0.5_x});
    ExactBoxType root_domain({{0.0_x,1.0_x}});
    FloatDPBoundsVector root=solver.solve(linear,root_domain);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(
        contains,ExactIntervalType(0.25_x,0.75_x),root[0]);
    config.set_verbosity(old_verbosity);
}


Void test_krawczyk_logging_paths()
{
    EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(1,0);
    EffectiveVectorMultivariateFunction linear({x-0.5_x});
    FloatDPBoundsVector search_box({FloatDPBounds(0.0_x,1.0_x,dp)});

    KrawczykSolver krawczyk(1e-8,12);
    FactoredKrawczykSolver factored(1e-8,12);
    auto& configuration=Logger::instance().configuration();
    unsigned int const previous_verbosity=configuration.verbosity();
    configuration.set_verbosity(10u);

    auto krawczyk_step=krawczyk.step(linear,search_box);
    auto factored_step=factored.step(linear,search_box);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(
        contains,ExactIntervalType(0.0_x,1.0_x),krawczyk_step[0]);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(
        contains,ExactIntervalType(0.0_x,1.0_x),factored_step[0]);

    auto coordinates=EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction implicit_equation({coordinates[1]-coordinates[0]});
    ExactBoxType parameter_box({{-0.25_x,0.25_x}});
    ExactBoxType solution_box({{-1.0_x,1.0_x}});
    auto implicit_solution=krawczyk.implicit(implicit_equation,parameter_box,solution_box);
    EffectiveScalarMultivariateFunction p=EffectiveScalarMultivariateFunction::coordinate(1,0);
    EffectiveVectorMultivariateFunction expected({p});
    ARIADNE_TESTING_EXPECT_LESS(norm(implicit_solution-expected), 1e-8_pr);

    configuration.set_verbosity(previous_verbosity);
}


// Exercise the public inline accessors and the exception overloads directly.
Void test_equation_solver_header_contracts()
{
    TestEquationSolverBase solver(TestEquationSolverBase::StepMode::Identity);
    solver.set_maximum_error(ApproximateDouble(0.125));
    ARIADNE_TESTING_EXPECT_EQUAL(solver.maximum_error(),cast_exact(ApproximateDouble(0.125)));
    solver.set_maximum_number_of_steps(17u);
    ARIADNE_TESTING_EXPECT_EQUAL(solver.maximum_number_of_steps(),17u);

    ARIADNE_TESTING_EXPECT_EQUAL(std::string(EquationSolverException("equation").what()),std::string("equation"));
    ARIADNE_TESTING_EXPECT_EQUAL(std::string(UnknownSolutionException("unknown").what()),std::string("unknown"));
    ARIADNE_TESTING_EXPECT_EQUAL(std::string(NoSolutionException("no solution").what()),std::string("no solution"));
    ARIADNE_TESTING_EXPECT_EQUAL(std::string(SingularJacobianException("singular").what()),std::string("singular"));

    // A generic equation solver exception must also catch the specialised cases.
    try {
        throw SingularJacobianException("singular Jacobian");
    } catch(const EquationSolverException& error) {
        ARIADNE_TESTING_EXPECT_EQUAL(std::string(error.what()),std::string("singular Jacobian"));
    }

    // Lexicographic ordering must first check dimensions, then compare
    // the first differing component; equal vectors are not less-than.
    FloatDPApproximationVector shorter(1u,dp), longer(2u,dp);
    ARIADNE_TESTING_PROVE(shorter<longer);
    ARIADNE_TESTING_PROVE(!(longer<shorter));
    FloatDPApproximationVector low(2u,dp), high(2u,dp);
    low[0]=FloatDPApproximation(1.0_x,dp);
    high[0]=FloatDPApproximation(1.0_x,dp);
    low[1]=FloatDPApproximation(2.0_x,dp);
    high[1]=FloatDPApproximation(3.0_x,dp);
    ARIADNE_TESTING_PROVE(low<high);
    ARIADNE_TESTING_PROVE(!(high<low));
    ARIADNE_TESTING_PROVE(!(low<low));
    // The Bounds instantiation used by the solver must cover unequal sizes.
    FloatDPBoundsVector shorter_bounds({FloatDPBounds(0,1,dp)});
    FloatDPBoundsVector longer_bounds({FloatDPBounds(0,1,dp),FloatDPBounds(1,2,dp)});
    ARIADNE_TESTING_PROVE(shorter_bounds<longer_bounds);
    ARIADNE_TESTING_PROVE(!(longer_bounds<shorter_bounds));
}

    Void test_maximum_error_constructors() {
        // Both constructors must preserve the exact stored double, not an
        // indeterminate comparison between approximate numeric objects.
        MaximumError maximum_error_from_approximate(ApproximateDouble(1e-5));
        MaximumError maximum_error_from_double(1e-5);
        ARIADNE_TESTING_EXPECT_EQUAL(maximum_error_from_approximate.value().get_d(), 1e-5);
        ARIADNE_TESTING_EXPECT_EQUAL(maximum_error_from_double.value().get_d(), 1e-5);
    }

    Void test_krawczyk_function_factory() {
        KrawczykSolver solver(maximum_error=1e-5_pr,maximum_number_of_steps=12);
        ARIADNE_TESTING_NOTE("solver.function_factory() = " << solver.function_factory());
        ARIADNE_TESTING_NOTE("INCOMPLETE");
    }
};

int main(int argc, char** argv)
{
    // Preserve logging coverage without flooding the test reporter.
#ifdef _WIN32
    Logger::instance().redirect_to_file("NUL");
#else
    Logger::instance().redirect_to_file("/dev/null");
#endif
    ARIADNE_TESTING_SUITE("test_equation_solvers");
    ARIADNE_TESTING_ADD_TEST(TestEquationSolverContracts, test_equation_solver_header_contracts);
    ARIADNE_TESTING_ADD_TEST(TestEquationSolverContracts, test_equation_solver_base_paths);
    ARIADNE_TESTING_ADD_TEST(TestEquationSolverContracts, test_interval_newton_paths);
    ARIADNE_TESTING_ADD_TEST(TestEquationSolverContracts, test_krawczyk_logging_paths);
    ARIADNE_TESTING_ADD_TEST(TestEquationSolverContracts, test_maximum_error_constructors);
    ARIADNE_TESTING_ADD_TEST(TestEquationSolverContracts, test_krawczyk_function_factory);

    auto interval_newton = Ariadne::Testing::fixture_variant<TestEquationSolver>(
        ariadne_test_suite, "IntervalNewtonSolver",
        [](Ariadne::Testing::Context& context) {
            return TestEquationSolver(context,
                IntervalNewtonSolver(maximum_error=1e-5_pr,maximum_number_of_steps=12),
                "IntervalNewtonSolver");
        });
    ARIADNE_TESTING_ADD_VARIANT_TEST(interval_newton, test_solve);
    ARIADNE_TESTING_ADD_VARIANT_TEST(interval_newton, test_implicit);
    ARIADNE_TESTING_ADD_VARIANT_TEST(interval_newton, test_scalar_implicit);

    auto krawczyk = Ariadne::Testing::fixture_variant<TestEquationSolver>(
        ariadne_test_suite, "KrawczykSolver",
        [](Ariadne::Testing::Context& context) {
            return TestEquationSolver(context,
                KrawczykSolver(maximum_error=1e-5_pr,maximum_number_of_steps=12),
                "KrawczykSolver");
        });
    ARIADNE_TESTING_ADD_VARIANT_TEST(krawczyk, test_solve);
    ARIADNE_TESTING_ADD_VARIANT_TEST(krawczyk, test_implicit);
    ARIADNE_TESTING_ADD_VARIANT_TEST(krawczyk, test_scalar_implicit);

    auto factored_krawczyk = Ariadne::Testing::fixture_variant<TestEquationSolver>(
        ariadne_test_suite, "FactoredKrawczykSolver",
        [](Ariadne::Testing::Context& context) {
            return TestEquationSolver(context,
                FactoredKrawczykSolver(maximum_error=1e-5_pr,maximum_number_of_steps=12),
                "FactoredKrawczykSolver");
        });
    ARIADNE_TESTING_ADD_VARIANT_TEST(factored_krawczyk, test_solve);
    ARIADNE_TESTING_ADD_VARIANT_TEST(factored_krawczyk, test_implicit);
    ARIADNE_TESTING_ADD_VARIANT_TEST(factored_krawczyk, test_scalar_implicit);

    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
