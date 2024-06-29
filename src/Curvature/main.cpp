#include <memory>
#include "Cuvature.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "========== Hw2 Usage  ==========\n";
        std::cout << std::endl;
        std::cout << "Input:	ACAM_mesh_HW2.exe	mesh.obj\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }

    //¶ÁÈëÍø¸ñ
    std::string mesh_path = argv[1];
    auto mesh = new PolyMesh();
    loadMesh(mesh_path, mesh);
    //create a basic mesh
//	testBasicTriangle();

    //mesh load an write , now only support obj/off
    //PolyMesh* mesh = new PolyMesh();
    //loadMesh("PumpkinMesh.obj", mesh);

    std::cout << "The curvature has area weight" << std::endl;
    std::vector<double> vertexLAR(mesh->numVertices(), 0.0);
    cal_local_ave_region(mesh, vertexLAR);
    cal_mean_curvature(mesh, vertexLAR);
    cal_gaussian_curvature(mesh, vertexLAR);

}
