#include "parametrization.h"

using namespace acamcad;
using namespace polymesh;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "========== parametrization  ==========\n";
        std::cout << std::endl;
        std::cout << "Input:	parametrization.exe	 mesh.obj\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }
    std::string mesh_path = argv[1];
    auto mesh {std::make_shared<PolyMesh> () } ;
    loadMesh(mesh_path, mesh.get());

    auto para_mesh = Jade::Param::Tutte_Embedding(mesh);

    std::cout<<"Finish Parameterization "<<std::endl;
}