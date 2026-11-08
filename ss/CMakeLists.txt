if(NOT TARGET function)
    ariadne_add_submodule(function)
endif()

# Testing is used only by the test executable.
if(NOT TARGET logging)
    ariadne_add_submodule(logging)
endif()

ariadne_setup_testing()
