#CMake file for Adding executables

function( ADD_RBM_TESTCASE arg1)

    add_executable(${arg1} ${arg1}.cpp)
    target_link_libraries(${arg1} RAZOR ${MANDATORY_LIBRARIES})

endfunction()

ADD_RBM_TESTCASE(razor)

