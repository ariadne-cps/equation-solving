cmake_minimum_required(VERSION 3.15)

project(EquationSolving VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT TARGET ariadne-configuration)
    find_package(Git REQUIRED)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}"
                    submodule update --init -- ss/configuration
                    RESULT_VARIABLE _result)
    if(_result)
        message(FATAL_ERROR "Cannot initialize configuration submodule")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/ss/configuration/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        CHAINS
            "configuration"
            "function|configuration"
    )
endif()

setup_project()

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(ss)

add_library(equation_solving INTERFACE)
target_include_directories(equation_solving INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-equation-solving>
)
target_link_libraries(equation_solving INTERFACE function)

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

setup_standalone_project_tests(
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
