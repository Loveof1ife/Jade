cmake_minimum_required(VERSION 3.5)

project(geom_proc)

#设置不生成MinSizeRel;RelWithDebInfo
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

# Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#path
SET(EIGEN_PATH "./include/Eigen")

############################################################
# Create a library
############################################################

#Generate the static library from the library sources
add_library(ACAM_mesh STATIC
        src/IOManger.cpp
        src/PolyMesh.cpp
        src/PolyMesh_Base.cpp
        src/PolyMeshIterators.cpp
)
target_include_directories(ACAM_mesh
        PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

############################################################
#1.Dijk_Prim
############################################################

file(GLOB Dijk_Prim
        src/Dijk_Prim/*.cpp
)

add_executable(dj ${Dijk_Prim} )

target_include_directories(dj
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(dj
        PRIVATE
        ACAM_mesh
)

