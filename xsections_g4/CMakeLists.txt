cmake_minimum_required(VERSION 3.20)
project(xsections_thin_target LANGUAGES CXX)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

if(NOT Geant4_multithreaded_FOUND)
  message(FATAL_ERROR "This project requires a multithreaded Geant4 build.")
endif()

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
add_executable(thin_target thin_target.cc ${sources})
target_include_directories(thin_target PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(thin_target PRIVATE ${Geant4_LIBRARIES})
target_compile_features(thin_target PRIVATE cxx_std_17)

add_executable(denominator_sampler denominator_sampler.cc
  src/ThinTargetDetector.cc)
target_include_directories(denominator_sampler PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(denominator_sampler PRIVATE ${Geant4_LIBRARIES})
target_compile_features(denominator_sampler PRIVATE cxx_std_17)
