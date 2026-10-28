#include <iostream>
#include <string>
#include "ariadne-equation-solving.hpp"
using namespace Ariadne;
using std::cout;
namespace {
void section(const std::string& title) { cout << "\n== " << title << " ==\n"; }
template<class T> void show(const std::string& label, const T& value) { cout << label << ": " << value << "\n"; }
}
int main() {
    section("Algebraic equations with Interval Newton");

    // EquationSolverBase exposes a common interface for algebraic solvers.  We use
    // IntervalNewtonSolver here as the representative implementation: the
    // tutorial is about the solving workflow, not comparing algorithms.
    IntervalNewtonSolver solver(1e-10, 16);
    show("maximum error", solver.maximum_error());
    show("maximum number of steps", solver.maximum_number_of_steps());

    auto x = EffectiveScalarMultivariateFunction::coordinate(1, 0);
    EffectiveVectorMultivariateFunction root_function({
        (x*x + 1) * x - 1
    });
    ExactBoxType root_domain({
        {0.0_x, 1.0_x}
    });

    // solve() searches for a unique transverse zero in the supplied box and
    // returns a validated enclosure of the root.
    FloatDPBoundsVector root = solver.solve(root_function, root_domain);
    show("root function", root_function);
    show("search domain", root_domain);
    show("validated root", root);

    // zero() is the explicit root-solving entry point.  It is useful when
    // writing code against EquationSolverInterface and makes the intended equation
    // f(x)=0 clear at the call site.
    FloatDPBoundsVector root_via_zero = solver.zero(root_function, root_domain);
    show("root via zero()", root_via_zero);

    section("Fixed points");

    // fixed_point() solves f(x)=x.  For f(x)=(x+1)/2 the unique fixed point
    // in [0,2] is x=1.
    EffectiveVectorMultivariateFunction contraction({
        (x + 1) / 2
    });
    ExactBoxType fixed_point_domain({
        {0.0_x, 2.0_x}
    });
    FloatDPBoundsVector fixed_point =
        solver.fixed_point(contraction, fixed_point_domain);

    show("fixed-point map", contraction);
    show("fixed-point enclosure", fixed_point);

    section("Implicit equations");

    // The implicit solver turns an equation F(a,y)=0 into a validated
    // function y=h(a).  Here F(a,y)=y^2-a and the selected branch is near
    // y=1, so the result is a validated square-root function on a parameter
    // interval around one.
    auto z = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction implicit_equation({
        z[1] * z[1] - z[0]
    });
    ExactBoxType parameter_domain({
        {0.75_x, 1.25_x}
    });
    ExactBoxType solution_range({
        {0.75_x, 1.25_x}
    });

    ValidatedVectorMultivariateFunctionPatch implicit_solution =
        solver.implicit(implicit_equation, parameter_domain, solution_range);

    show("implicit equation", implicit_equation);
    show("parameter domain", parameter_domain);
    show("solution search range", solution_range);
    show("implicit solution patch", implicit_solution);
    show("implicit solution range", implicit_solution.range());

    section("Krawczyk method");
    KrawczykSolver krawczyk(1e-10,16);
    show("Krawczyk validated root",krawczyk.solve(root_function,root_domain));
    FactoredKrawczykSolver factored(1e-10,16);
    show("Factored Krawczyk validated root",factored.solve(root_function,root_domain));
    return 0;
}
