#!/usr/bin/env python3
from pyariadne import *

def section(title):
    print(f"\n== {title} ==")

def show(label, value):
    print(f"{label}: {value}")

def exact_box(bounds):
    return FloatDPExactBox([{exact(lower): exact(upper)} for lower, upper in bounds])

def tutorial_equation_solving():
    section("Algebraic equations with Interval Newton")

    solver = IntervalNewtonSolver(1e-10, 16)
    show("solver", solver)
    show("maximum error", solver.maximum_error())
    show("maximum number of steps", solver.maximum_number_of_steps())

    x = EffectiveScalarMultivariateFunction.coordinate(1, 0)
    root_function = EffectiveVectorMultivariateFunction([
        (x * x + 1) * x - 1,
    ])
    root_domain = exact_box([
        (0, 1),
    ])

    root = solver.solve(root_function, root_domain)
    show("root function", root_function)
    show("search domain", root_domain)
    show("validated root", root)

    root_via_zero = solver.zero(root_function, root_domain)
    show("root via zero()", root_via_zero)

    section("Fixed points")

    contraction = EffectiveVectorMultivariateFunction([
        (x + 1) / 2,
    ])
    fixed_point_domain = exact_box([
        (0, 2),
    ])
    fixed_point = solver.fixed_point(contraction, fixed_point_domain)

    show("fixed-point map", contraction)
    show("fixed-point enclosure", fixed_point)

    section("Implicit equations")

    a = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    y = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    implicit_equation = EffectiveVectorMultivariateFunction([
        y * y - a,
    ])
    parameter_domain = exact_box([
        (0.75, 1.25),
    ])
    solution_range = exact_box([
        (0.75, 1.25),
    ])

    implicit_solution = solver.implicit(
        implicit_equation,
        parameter_domain,
        solution_range,
    )

    show("implicit equation", implicit_equation)
    show("parameter domain", parameter_domain)
    show("solution search range", solution_range)
    show("implicit solution patch", implicit_solution)
    show("implicit solution range", implicit_solution.range())

    section("Krawczyk methods")
    show("Krawczyk root", KrawczykSolver(1e-10, 16).solve(root_function, root_domain))
    show("Factored Krawczyk root", FactoredKrawczykSolver(1e-10, 16).solve(root_function, root_domain))

if __name__ == "__main__":
    tutorial_equation_solving()
