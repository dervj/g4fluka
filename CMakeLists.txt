cmake_minimum_required(VERSION 3.16)
project(G4FlukaComparison)

set(CMAKE_CXX_STANDARD 17)

# Find Geant4 — adjust the path if needed
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# Headers
include_directories(${PROJECT_SOURCE_DIR}/include)

# Sources
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)

add_executable(g4fluka g4fluka.cc ${SOURCES})
target_link_libraries(g4fluka ${Geant4_LIBRARIES})

# Copy macros to build directory
file(GLOB MACROS macros/*.mac)
foreach(_mac ${MACROS})
    configure_file(${_mac} ${PROJECT_BINARY_DIR} COPYONLY)
endforeach()
