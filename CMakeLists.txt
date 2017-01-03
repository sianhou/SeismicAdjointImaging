cmake_minimum_required(VERSION 3.0)
project(SeismicAdjointImaging)
file(MAKE_DIRECTORY ${SeismicAdjointImaging_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${SeismicAdjointImaging_BINARY_DIR}/bin)

#! Output path
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SeismicAdjointImaging_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SeismicAdjointImaging_BINARY_DIR}/bin)

#--------------------------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------------------------
add_definitions("-O3")

#! Library
set(SJINC lib/sjinc.h)
set(SJLIB lib/sjabc.c lib/sjfile.c lib/sjmalloc.c lib/sjmath.c lib/sjsimulation.c)
add_library(sj ${SJLIB} ${SGLIB})

#! RUNTIME
find_library(sj ${SeismicAdjointImaging_BINARY_DIR}/lib)

#! Test
set(TEST_ASGFD test/test_asgfd.c)
add_executable(test_asgfd ${TEST_ASGFD} ${SJINC})
target_link_libraries(test_asgfd sj m)

#! sgbin2su
set(SOURCEFILE src/sjbin2su.c)
add_executable(sjbin2su ${SOURCEFILE} ${SJINC})
target_link_libraries(sjbin2su sj m)
#! sgdisplaysurvey2d
set(SOURCEFILE src/sjdisplaysurvey2d.c)
add_executable(sjdisplaysurvey2d ${SOURCEFILE} ${SJINC})
target_link_libraries(sjdisplaysurvey2d sj m)
#! sgsurvey2d
set(SOURCEFILE src/sjsurvey2d.c)
add_executable(sjsurvey2d ${SOURCEFILE} ${SJINC})
target_link_libraries(sjsurvey2d sj m)

#--------------------------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_C_COMPILER "mpicc")

find_package ( MPI )
if ( MPI_FOUND )
    #! Library
    set(SJINC lib/sjinc.h)
    set(SJLIB lib/sjabc.c lib/sjfile.c lib/sjmalloc.c lib/sjmath.c lib/sjsimulation.c)
    add_library(mpisj ${SJINC} ${SJLIB})

    #! RUNTIME
    find_library(mpisj ${SeismicAdjointImaging_BINARY_DIR}/lib)

    #! sgmpiawsgfd2d
    set(SOURCEFILE src/sjmpiawsgfd2d.c)
    add_executable(sjmpiawsgfd2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiawsgfd2d mpisj m)

    #! sgmpiawfd2d
    set(SOURCEFILE src/sjmpiawfd2d.c)
    add_executable(sjmpiawfd2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiawfd2d mpisj m)

    #! sgmpiartm2d
    set(SOURCEFILE src/sjmpiartm2d.c)
    add_executable(sjmpiartm2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpiartm2d mpisj m)

    #! sgmpilsartm2d
    set(SOURCEFILE src/sjmpilsartm2d.c)
    add_executable(sjmpilsartm2d ${SOURCEFILE} ${SJINC})
    target_link_libraries(sjmpilsartm2d mpisj m)
endif( MPI_FOUND )