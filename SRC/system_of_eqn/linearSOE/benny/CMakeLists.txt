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

    PUBLIC
        BennySparseLinSOE.h
        BennySparseLinSolver.h

)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

