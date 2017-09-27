cmake_minimum_required(VERSION 3.5)
project(BebProject)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

# For 1-D only example uncomment below
#set(SOURCE_FILES src/mca_1Donly/EulersMethod1D.cpp src/mca_1Donly/EulersMethod1D.h
#        src/mca_1Donly/MarkovChainApproximation1D.cpp src/mca_1Donly/MarkovChainApproximation1D.h
#        examples/mca_1Donly/Functions.h examples/mca_1Donly/Functions.cpp
#        examples/mca_1Donly/OneDimension_Script.cpp)

include_directories(src)
include_directories(examples)

set(SOURCE_FILES
        src/mca/MarkovChainApproximation.cpp
        src/mca/MarkovChainApproximation.h
        src/mca/MarkovChainParameters.h
        src/mca/MarkovChainParameters.cpp
        src/eulersmethod/EulersMethod.cpp
        src/eulersmethod/EulersMethod.h
        src/eulersmethod/EulerParameters.cpp
        src/eulersmethod/EulerParameters.h
        src/multidim/GridIndex.cpp
        src/multidim/GridIndex.h
        examples/mca/Functions2.h
        examples/mca/Functions2.cpp
        examples/mca/OneDimensionPlus_Script.cpp src/multidim/GridParameters.cpp src/multidim/GridParameters.h)

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
