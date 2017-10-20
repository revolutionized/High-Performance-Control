cmake_minimum_required(VERSION 3.5)
project(hpc)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

include_directories(examples)
set(EXAMPLE_FILES
        examples/onedimensional/OneDimension_Functions.h
        examples/onedimensional/OneDimension_Functions.cpp
        examples/onedimensional/OneDimensionPlus_Script.cpp
        examples/onedimensional/OneDimensionPlus_Script.h
        examples/dubinscar/dubinscar.cpp
        examples/dubinscar/dubinscar.h)

include_directories(3rdparty)

# For 1-D only example uncomment below
#set(SOURCE_FILES src/onedimensional_OLD/EulersMethod1D.cpp src/onedimensional_OLD/EulersMethod1D.h
#        src/onedimensional_OLD/MarkovChainApproximation1D.cpp src/onedimensional_OLD/MarkovChainApproximation1D.h
#        examples/onedimensional_OLD/Functions.h examples/onedimensional_OLD/Functions.cpp
#        examples/onedimensional_OLD/OneDimension_Script.cpp)

include_directories(src)
set(SOURCE_FILES
        src/main/main.cpp
        src/mca/MarkovChainApproximation.cpp
        src/mca/MarkovChainApproximation.h
        src/mca/MarkovChainParameters.h
        src/mca/MarkovChainParameters.cpp
        src/eulersmethod/EulersMethod.cpp
        src/eulersmethod/EulersMethod.h
        src/multidim/GridIndex.cpp
        src/multidim/GridIndex.h
        src/multidim/GridParameters.cpp
        src/multidim/GridParameters.h
        src/utilities/ProgressBar.h)



# To use the C3 library, uncomment below
#set(USEC3 TRUE)
if (${USEC3})
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
endif (${USEC3})

# To convert data to MATLAB files uncomment below and library links for m, mat and mx
#set(MATLAB_DIR /usr/local/MATLAB/R2017a)
#include_directories(${MATLAB_DIR}/extern/include)
#include_directories(${MATLAB_DIR}/bin/glnxa64/)
#link_directories(${MATLAB_DIR}/bin/glnxa64)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${EXAMPLE_FILES} ${C_SOURCE_FILES})

if (${USEC3})
    target_link_libraries(${PROJECT_NAME} c3)
    target_link_libraries(${PROJECT_NAME} c3sc)
    target_link_libraries(${PROJECT_NAME} cdyn)
endif (${USEC3})

#target_link_libraries(${PROJECT_NAME} m)
#target_link_libraries(${PROJECT_NAME} mat)
#target_link_libraries(${PROJECT_NAME} mx)