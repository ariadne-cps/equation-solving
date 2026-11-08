cmake_minimum_required(VERSION 3.15)

project(EquationSolving VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT TARGET ariadne-build)
    find_package(Git REQUIRED)
    # Keep the initial Build checkout output in English.
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env "LC_ALL=C"
                "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}"
                submodule update --init -- ss/bui
        RESULT_VARIABLE _result)
    if(_result)
        message(FATAL_ERROR "Cannot initialize build submodule")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/ss/bui/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Pinned child commits guarantee their internal consistency.
    # Testing is development-only and excluded from dependency consistency.
    require_same_dependency_commit(Build CHAINS
        "build" "function|build" "logging|build")
endif()

setup_project()

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(ss)

add_library(equation_solving INTERFACE)
target_include_directories(equation_solving INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-equation-solving>
)
target_link_libraries(equation_solving INTERFACE function logging)

add_subdirectory(src)
target_link_libraries(EQUATION_SOLVING_SRC PRIVATE equation_solving)

ariadne_register_public_headers(
    TARGET equation_solving
    DIRECTORIES "${PROJECT_SOURCE_DIR}/include/equation-solving"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-equation-solving.hpp"
)

setup_project_library(
    TARGET ariadne-equation-solving
    OBJECTS EQUATION_SOLVING_SRC
    LINK_LIBRARIES equation_solving
)

setup_standalone_project(
    TARGET ariadne-equation-solving
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    ariadne_install_dependency_bundle(
        TARGET equation_solving
        DESTINATION include/ariadne-equation-solving
    )
endif()
