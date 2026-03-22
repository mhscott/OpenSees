#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_SysOfEqn
    PRIVATE
        BennySparseLinSOE.cpp
        BennySparseLinSolver.cpp
        benny_sparse.cpp
        benny_symbolic.cpp
        benny_numeric.cpp

    PUBLIC
        BennySparseLinSOE.h
        BennySparseLinSolver.h
        benny_sparse.hpp
        benny_symbolic.hpp
        benny_numeric.hpp

)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

