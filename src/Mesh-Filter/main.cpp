#include <memory>
#include "Filter.h"

using namespace acamcad;
using namespace polymesh;


int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "========== Mesh Filter ==========\n";
        std::cout << std::endl;
        std::cout << "Input: Filter.exe	mesh.obj\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }
        std::string mesh_path = argv[1];

        auto in_mesh{ std::make_shared<PolyMesh>()} ;
        auto out_mesh{ std::make_shared<PolyMesh>()};

        std::vector<MVector3> filterNormalField;
        std::vector<MPoint3> vertexPositions;
        std::vector<MPoint3> faceCenter;
        std::vector<double> faceArea;

        loadMesh(mesh_path, in_mesh.get());

        size_t fN = in_mesh->numPolygons();
        size_t fV = in_mesh->numVertices();
        filterNormalField.resize(fN), faceCenter.resize(fN), faceArea.resize(fN), vertexPositions.resize(fV);

        Jade::Filter::BilateralNormalFiltering(in_mesh, filterNormalField, faceCenter, faceArea);
        Jade::Filter::fitVertexToNormal(vertexPositions, in_mesh, filterNormalField, faceCenter, faceArea);

        std::cout<<"Filtering done"<<std::endl;
}