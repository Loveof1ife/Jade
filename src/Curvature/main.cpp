#include <memory>
#include "Cuvature.h"
#include "BaseFunc.h"

using namespace Jade::Curvature;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "========== Curvature  ==========\n";
        std::cout << std::endl;
        std::cout << "Input:	curvature.exe	mesh.obj\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }

    std::string mesh_path = argv[1];
    auto mesh = std::make_shared<acamcad::polymesh::PolyMesh>();
    if (!loadMesh(mesh_path, mesh.get())) {
        std::cerr << "Error loading mesh: " << mesh_path << std::endl;
        return 1;
    }

    std::cout << "The curvature has area weight" << std::endl;
    size_t N = mesh->numVertices();
    std::vector<double> vertexLocalAveRegion(N, 0.0);
    std::vector<CurvatureInfo> curvatures(N, CurvatureInfo{});

    CalLocalAveRegion(mesh, vertexLocalAveRegion );
    CalMeanCurvature(mesh, vertexLocalAveRegion, curvatures);
    CalGaussianCurvature(mesh, vertexLocalAveRegion, curvatures);
}
