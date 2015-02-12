project(commonlib)
cmake_minimum_required(VERSION 2.6)

enable_language(C)
SET(CMAKE_COMPILER_IS_GNUCC) 
SET(CMAKE_CXX_FLAGS "-lm")

INCLUDE(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-std=c99" HAVE_C99)

IF (HAVE_C99)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
ELSE(HAVE_C99)
  MESSAGE(WARNING "C99 support not available, program may fail to compile!")
ENDIF(HAVE_C99)



add_library(common common.c)
add_executable(gensumserial genSumVec.c)
target_link_libraries(gensumserial common)

enable_testing()

add_test(gensumserial ${PROJECT_BINARY_DIR}/gensumserial)
#set_tests_properties(gensumserial PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)