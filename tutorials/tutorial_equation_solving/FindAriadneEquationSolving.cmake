find_library(ARIADNE_EQUATION_SOLVING_LIBRARY NAMES ariadne-equation-solving)

find_package(Threads REQUIRED)

find_path(ARIADNE_EQUATION_SOLVING_INCLUDE_DIR ariadne-equation-solving.hpp PATH_SUFFIXES ariadne-equation-solving)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneEquationSolving DEFAULT_MSG
  ARIADNE_EQUATION_SOLVING_LIBRARY
  ARIADNE_EQUATION_SOLVING_INCLUDE_DIR
)

if(AriadneEquationSolving_FOUND)
  get_filename_component(ARIADNE_EQUATION_SOLVING_INCLUDE_PARENT_DIR ${ARIADNE_EQUATION_SOLVING_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_EQUATION_SOLVING_INCLUDE_DIRS
    ${ARIADNE_EQUATION_SOLVING_INCLUDE_PARENT_DIR}
    ${ARIADNE_EQUATION_SOLVING_INCLUDE_DIR}
  )
  if(NOT TARGET AriadneEquationSolving::ariadne-equation-solving)
    add_library(AriadneEquationSolving::ariadne-equation-solving UNKNOWN IMPORTED)
    set_target_properties(AriadneEquationSolving::ariadne-equation-solving PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_EQUATION_SOLVING_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_EQUATION_SOLVING_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "Threads::Threads"
    )
  endif()
endif()
