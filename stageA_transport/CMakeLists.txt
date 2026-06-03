cmake_minimum_required(VERSION 3.20)
project(stageA_transport LANGUAGES CXX)

# --- Geant4 (multithreaded; no UI/vis components needed for batch) ----------
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

if(NOT Geant4_multithreaded_FOUND)
  message(FATAL_ERROR "This project requires a multithreaded Geant4 build.")
endif()

# --- sources ----------------------------------------------------------------
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
# Shared isotope table + (later) HDF5 schema helpers live in ../common.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../common)

add_executable(proton_transport proton_transport.cc ${sources})
target_link_libraries(proton_transport PRIVATE ${Geant4_LIBRARIES})
target_compile_features(proton_transport PRIVATE cxx_std_17)

# Copy macros next to the binary for convenience.
file(GLOB macros ${CMAKE_CURRENT_SOURCE_DIR}/macros/*.mac)
foreach(_m ${macros})
  configure_file(${_m} ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
endforeach()
