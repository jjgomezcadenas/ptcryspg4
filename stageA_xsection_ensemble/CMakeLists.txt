cmake_minimum_required(VERSION 3.20)
project(stageA_xsection_ensemble LANGUAGES CXX)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

if(NOT Geant4_multithreaded_FOUND)
  message(FATAL_ERROR "This project requires a multithreaded Geant4 build.")
endif()

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
add_executable(xsection_ensemble xsection_ensemble.cc ${sources})
target_include_directories(xsection_ensemble PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(xsection_ensemble PRIVATE ${Geant4_LIBRARIES})
target_compile_features(xsection_ensemble PRIVATE cxx_std_17)
