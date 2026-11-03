if(NOT TARGET function)
    ariadne_add_submodule(function)
endif()

# Testing is used only by the test executable.
ariadne_setup_testing()
