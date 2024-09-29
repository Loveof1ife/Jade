cmake_minimum_required(VERSION 3.5)

project(geom_proc)

#设置不生成MinSizeRel;RelWithDebInfo
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

# Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#path
SET(EIGEN_PATH "./include/Eigen")

include_directories(${PROJECT_SOURCE_DIR}/include)

############################################################
# Create a library
############################################################

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

function(add_dgp_executable target_name src_path)
    file(GLOB TARGET_SRC
            ${src_path}/*.cpp
    )

    add_executable(${target_name} ${TARGET_SRC})

    target_include_directories(${target_name}
            PRIVATE ${PROJECT_SOURCE_DIR}/include
    )

    target_link_libraries(${target_name}
            PRIVATE ACAM_mesh
    )
endfunction()

############################################################
#add digital geometry process exe
############################################################
# Dijk_Prim executable
add_dgp_executable(dj src/Dijk_Prim)

# Curvature executable
add_dgp_executable(curvature src/Curvature)

# Filter executable
add_dgp_executable(filter src/Mesh-Filter)

#
add_dgp_executable(parametrization src/Parametrization)



