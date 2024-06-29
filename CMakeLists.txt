cmake_minimum_required(VERSION 3.5)

project(geom_proc)

#设置不生成MinSizeRel;RelWithDebInfo
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")


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

add_executable(dijk_prim ${Dijk_Prim}
        src/Dijk_Prim/main.cpp)

target_include_directories(dijk_prim
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include

)

target_link_libraries(dijk_prim
        PRIVATE
        ACAM_mesh
)
############################################################
#2.Curvature
############################################################

file(GLOB Curvature
        src/Curvature/*.cpp
)

add_executable(curvature ${Curvature}
)

target_include_directories(curvature
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include

)
target_link_libraries(curvature
        PRIVATE
        ACAM_mesh
)

############################################################
#3.B_Filter
############################################################

file(GLOB Bilateral_Normal_Filtering
        src/Bilateral_Normal_Filtering/*.cpp
)

add_executable(bilateral_normal_filtering ${Bilateral_Normal_Filtering}
)

target_include_directories(bilateral_normal_filtering
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include

)
target_link_libraries(bilateral_normal_filtering
        PRIVATE
        ACAM_mesh
)

############################################################
#4.parameter
############################################################

file(GLOB Parameterization
        src/Parameterizations/*.cpp
)

add_executable(parameterization ${Parameterization}
)

target_include_directories(parameterization
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include

)
target_link_libraries(parameterization
        PRIVATE
        ACAM_mesh
)

############################################################
#5.As_Rigid_As_Posibile deformation
############################################################
file(GLOB Asag_Deformation
        src/Asag_Deformation/*.cpp
)
add_executable(asag_deformation ${Asag_Deformation})
target_include_directories(asag_deformation
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${EIGEN_PATH}
)

target_link_libraries(asag_deformation
        PRIVATE
        ACAM_mesh
)
############################################################
#5.As-Rigid-As-Possible Shape Interpolation
############################################################\
file(GLOB rigid_shape_interpolation
        src/rigid_shape_interpolation/*.cpp
)
add_executable(rigid_shape_interpolation_ ${rigid_shape_interpolation})
target_include_directories(rigid_shape_interpolation_
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include/
        ${EIGEN_PATH}
)

target_link_libraries(rigid_shape_interpolation_
        PRIVATE
        ACAM_mesh
)

############################################################
#5.QEM, Incremental_Remeshing, VSA
############################################################\
file(GLOB MeshSimplification
        src/Simplification_Remeshing/QEM.cpp
)
add_executable(MeshSimplification_ ${MeshSimplification})
target_include_directories(MeshSimplification_
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include/
        ${EIGEN_PATH}
)

target_link_libraries(MeshSimplification_
        PRIVATE
        ACAM_mesh
)

file(GLOB Remeshing
        src/Simplification_Remeshing/Incremental_Remeshing.cpp
)
add_executable(Remeshing_ ${Remeshing})
target_include_directories(Remeshing_
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include/
        ${EIGEN_PATH}
)

target_link_libraries(Remeshing_
        PRIVATE
        ACAM_mesh
)

file(GLOB Variational_Shape_Approximation
        src/Simplification_Remeshing/Variational_Shape_Approximation.cpp
)
add_executable(Variational_Shape_Approximation_ ${Variational_Shape_Approximation})
target_include_directories(Variational_Shape_Approximation_
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include/
        ${EIGEN_PATH}
)

target_link_libraries(Variational_Shape_Approximation_
        PRIVATE
        ACAM_mesh
)

file(GLOB Optimal_Delaunay_Triangulation
        src/Simplification_Remeshing/Optimal_Delaunay_Triangulation.cpp
)
add_executable(Optimal_Delaunay_Triangulation_ ${Optimal_Delaunay_Triangulation})
target_include_directories(Optimal_Delaunay_Triangulation_
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include/
        ${EIGEN_PATH}
)

target_link_libraries(Optimal_Delaunay_Triangulation_
        PRIVATE
        ACAM_mesh
)

############################################################
#6.Field
############################################################\
file(GLOB Field
        src/Field/*.cpp
)
add_executable(Field_ ${Field})
target_include_directories(Field_
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include/
        ${EIGEN_PATH}
)

target_link_libraries(Field_
        PRIVATE
        ACAM_mesh
)

