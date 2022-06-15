project(ntHits)
cmake_minimum_required(VERSION 3.22)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories("vendor/btllib/include"
        "vendor/btllib/subprojects/cpptoml/include"
        "vendor/btllib/subprojects/sdsl-lite/include")

include_directories(Common)
include_directories(vendor/ntHash)
include_directories(.)

add_executable(ntHits Common/Fcontrol.cpp Common/SignalHandler.cpp nthits.cpp)
