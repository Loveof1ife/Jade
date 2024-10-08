#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "PolyMesh/IOManager.h"
#include <string>
#include <memory>

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

using namespace acamcad;
using namespace polymesh;

namespace Jade::Param{
    using SparseMatrix = Eigen::SparseMatrix<double>;
    std::shared_ptr<PolyMesh> Local_Global_Approach(const std::shared_ptr<const PolyMesh>& mesh );

    class Parameterizer{
    public:
        explicit Parameterizer(const std::string& meshPath);
        void tutteEmbedding(){
            area = calPolymeshArea();
            auto boundary_of_uv = setBoundaryVerticesToCircle();
            solveUV(boundary_of_uv);
        }

        void arapEmbedding() {
            // Step 1: Tutte initialization
            tutteEmbedding();
            // Step 2: Compute local coordinates (cotangent weights)
            computeLocalCoordinates();
            computeCotangentWeights();  // Compute cotangent weights
            // Step 3: Iteratively solve local-global ARAP
            localGlobalArapIterations(100);  // Fixed number of iterations
            // After solving, update the mesh with new UV positions
            auto uvMesh = applyUVtoMesh();
        }

        Eigen::VectorXd  setBoundaryVerticesToCircle();
        void solveUV(Eigen::VectorXd boundary_of_uv);

        void computeCotangentWeights();
        void computeLocalCoordinates();

        void localGlobalArapIterations(int num_iterations);
        void solveLocalStep();
        void solveGlobalStep();

        [[nodiscard]] double calPolymeshArea() const;
        [[nodiscard]] auto getBoundaryVertices() {return  mesh->boundaryVertices();}
        [[nodiscard]] auto getNumofVertices() {return  mesh->numVertices();}
        [[nodiscard]] std::shared_ptr<PolyMesh> applyUVtoMesh() const;
    private:
        std::shared_ptr<PolyMesh> mesh{nullptr};
        std::shared_ptr<PolyMesh> param_mesh{nullptr};
        Eigen::MatrixX2d tutteParameterization{};
        Eigen::MatrixX2d iterParameterization{};
        std::vector<MVert* > boundary_vertices{};
        Eigen::MatrixXd localCoord{};  // Local coordinate matrix
        std::vector<Eigen::Matrix2d> rotation_matrices{};  // Rotation matrices for ARAP

        std::vector<double> cotangent_weights{};  // Cotangent weights for mesh
        Eigen::SparseMatrix<double> laplaceMatrix{};
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver{} ;
        double area{};
        size_t n_vertices{};
        size_t n_faces{};
    };
}
