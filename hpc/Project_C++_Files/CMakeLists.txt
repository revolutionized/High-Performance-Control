cmake_minimum_required(VERSION 3.5)
project(BebProject)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        EulersMethod1D.cpp
        EulersMethod1D.h
        MarkovChainApproximation1D.cpp
        MarkovChainApproximation1D.h
#        DotMatCreator.cpp
#        DotMatCreator.h
        Functions.h
        Functions.cpp 
        OneDimension_Script.cpp
        MarkovChainApproximation.cpp
        MarkovChainApproximation.h
        Functions2.h
        Functions2.cpp MarkovChainParameters.h MarkovChainParameters.cpp OneDimensionPlus_Script.cpp EulersMethod.cpp EulersMethod.h EulerParameters.cpp EulerParameters.h)

# Include header files for MATLAB (to create .mat files)
#set(MATLAB_DIR /usr/local/MATLAB/R2017a)
#include_directories(${MATLAB_DIR}/extern/include)
# Incude MATLAB libraries directory
#include_directories(${MATLAB_DIR}/bin/glnxa64/)
#link_directories(${MATLAB_DIR}/bin/glnxa64)

# Include all header files of the C3 lib
set(C3_DIR ../../c3)
set(C3_DIR_SRC ${C3_DIR}/src)
include_directories(${C3_DIR_SRC}/lib_array)
include_directories(${C3_DIR_SRC}/lib_clinalg)
include_directories(${C3_DIR_SRC}/lib_fft)
include_directories(${C3_DIR_SRC}/lib_funcs)
include_directories(${C3_DIR_SRC}/lib_interface)
include_directories(${C3_DIR_SRC}/lib_linalg)
include_directories(${C3_DIR_SRC}/lib_optimization)
include_directories(${C3_DIR_SRC}/lib_probability)
include_directories(${C3_DIR_SRC}/lib_quadrature)
include_directories(${C3_DIR_SRC}/lib_stringmanip)
include_directories(${C3_DIR_SRC}/lib_superlearn)
include_directories(${C3_DIR_SRC}/lib_tensdecomp)
include_directories(${C3_DIR_SRC}/lib_tensor)
# Include library and main c3 locations
include_directories(${C3_DIR}/include)
# Include link of C3 library
link_directories(${C3_DIR}/lib)

# Include all header files of the C3SC lib
set(C3SC_DIR ../../c3sc)
include_directories(${C3SC_DIR}/src)
# Include library for c3sc
link_directories(${C3SC_DIR}/lib)

# Include all header files of the CDYN lib
set(CDYN_DIR ../../cdyn)
include_directories(${CDYN_DIR}/src)
# Include library for cdyn
link_directories(${CDYN_DIR}/lib)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${C_SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} c3)
target_link_libraries(${PROJECT_NAME} c3sc)
target_link_libraries(${PROJECT_NAME} cdyn)
#target_link_libraries(${PROJECT_NAME} m)
#target_link_libraries(${PROJECT_NAME} mat)
#target_link_libraries(${PROJECT_NAME} mx)
