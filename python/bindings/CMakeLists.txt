include_directories(${Python3_INCLUDE_DIRS})

get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-equation_solving-submodule-obj OBJECT
    equation_solving_submodule.cpp
)
target_include_directories(pyariadne-equation_solving-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    $<TARGET_PROPERTY:pyariadne-function,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-equation_solving-submodule-obj PRIVATE
    pyariadne-python-common
    pybind11::module
    equation_solving
)
set_target_properties(pyariadne-equation_solving-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
    target_compile_options(pyariadne-equation_solving-submodule-obj PRIVATE
        /bigobj
        /wd4702
    )
endif()

add_library(pyariadne-equation_solving INTERFACE)
target_sources(pyariadne-equation_solving INTERFACE
    $<TARGET_OBJECTS:pyariadne-equation_solving-submodule-obj>
)
target_include_directories(pyariadne-equation_solving INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-function,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-equation_solving INTERFACE
    pyariadne-function
)

ariadne_add_python_module(
    NAME equation_solving
    MODULE_SOURCE equation_solving_module.cpp
)
