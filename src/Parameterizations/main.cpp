#include <memory>
#include "parameterizations.h"

using namespace acamcad;
using namespace polymesh;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "========== Hw4 Usage  ==========\n";
        std::cout << std::endl;
        std::cout << "Input:	ACAM_mesh_HW4.exe	mesh.obj\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }
    std::string mesh_path = argv[1];
    auto mesh = new PolyMesh();
    loadMesh(mesh_path, mesh);
    Tutte_Embedding(mesh);

    std::cout<<"Finish Parameterization "<<std::endl;
}