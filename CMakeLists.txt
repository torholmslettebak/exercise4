project(commonlib)
cmake_minimum_required(VERSION 2.6)

enable_language(C)
SET(CMAKE_COMPILER_IS_GNUCC) 
# SET(CMAKE_C_FLAGS "-lm")


INCLUDE(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-std=c99" HAVE_C99)

IF (HAVE_C99)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 ")

ELSE(HAVE_C99)
  MESSAGE(WARNING "C99 support not available, program may fail to compile!")
ENDIF(HAVE_C99)

find_package(MPI REQUIRED)
OPTION(USE_OPENMP "Use OpenMP?" ON)
OPTION(USE_MPI "Use MPI?" ON)

IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP)
ENDIF(USE_OPENMP)
IF(USE_MPI)
  FIND_PACKAGE(MPI)
ENDIF(USE_MPI)


include_directories( ${MPI_C_INCLUDE_PATH})
if (MPI_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_FLAGS}")
endif()
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()



add_library(common common.c)
add_executable(gensumserial genSumVec.c)
target_link_libraries(gensumserial common)

add_executable(sumvec_openmp sumvec_openmp.c)
find_library(M_LIB m)
target_link_libraries(sumvec_openmp ${M_LIB})
target_link_libraries(sumvec_openmp common)

add_executable(sumvec_MPI sumvec_MPI.c)
target_link_libraries(sumvec_MPI ${MPI_C_LIBRARIES})
target_link_libraries(sumvec_MPI common)
target_link_libraries(sumvec_MPI ${M_LIB})

add_executable(sumvec_comb sumvec_comb.c)
target_link_libraries(sumvec_comb ${MPI_C_LIBRARIES})
target_link_libraries(sumvec_comb common)
target_link_libraries(sumvec_comb ${M_LIB})

enable_testing()


add_test(gensumserial ${PROJECT_BINARY_DIR}/gensumserial)
add_test(sumvec_openmp ${PROJECT_BINARY_DIR}/sumvec_openmp)
add_test(sumvec_MPI ${PROJECT_BINARY_DIR}/sumvec_MPI)
add_test(sumvec_comb ${PROJECT_BINARY_DIR}/sumvec_comb)
set_tests_properties(sumvec_openmp PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)