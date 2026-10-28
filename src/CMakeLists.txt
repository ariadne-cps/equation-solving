add_library(EQUATION_SOLVING_SRC OBJECT
    equation_solver.cpp
    interval_newton_method.cpp
    krawczyk_method.cpp
)
target_include_directories(EQUATION_SOLVING_SRC PRIVATE "${PROJECT_SOURCE_DIR}/include/equation-solving")
