#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#include "Algorithm.h"

int main(int argc, char** argv){
        //alien.obj 1 v.txt path.txt
        if (argc != 4)
        {
            std::cout << "========== ShortestPath && SteinerTree ==========\n";
            std::cout << std::endl;
            std::cout << "ShortestPath:	dj.exe	alien.obj	vertices.txt	path.txt\n";
            std::cout << "SteinerTree:	dj.exe	alien.obj	vertices.txt	path.txt\n";
            std::cout << std::endl;
            std::cout << "=================================================\n";
            return 0;
        }
        std::string mesh_path = argv[1];
        auto mesh = std::make_unique<acamcad::polymesh::PolyMesh>();
        if (!loadMesh(mesh_path, mesh.get())) {
            std::cerr << "Error loading mesh: " << mesh_path << std::endl;
            return 1;
        }

        std::cout << "Number of vertices in mesh: " << mesh->numVertices() << std::endl;

        std::vector<int> landmarks;

        std::ifstream vertices_in(argv[2]);

        if (!vertices_in.is_open()) {
            std::cerr << "Error opening vertices file: " << argv[2] << std::endl;
            return 1;
        }
        std::string line;

        while (std::getline(vertices_in, line)){
            std::stringstream ss;
            ss << line;
            int vid;
            ss >> vid;
            landmarks.push_back(vid);
        }

        vertices_in.close();

        std::vector<int> vertices{}, edges{};
        std::vector<std::pair<int, int> > tree_edge{};
        if (landmarks.size() == 2){
            std::vector<int> vertexPath;

            Jade::ShortPath::Dijkstra(mesh, landmarks[0], landmarks[1], vertexPath);

            vertices.push_back(vertexPath[0]);

            for (int i = 0; i < vertexPath.size() - 1; i++){
                vertices.push_back(vertexPath[i + 1]);
                edges.push_back(mesh->edgeBetween(mesh->vert(vertexPath[i]), mesh->vert(vertexPath[i+1]))->index());
            }
        }
        else if (landmarks.size()>2){

//            Jade::ShortPath::SteinerTree(mesh, landmarks, tree_edge);

            std::vector<std::vector<int> > vertexPaths;

            Jade::ShortPath::DijkstraGroup(mesh, landmarks, vertexPaths);

            for (auto path : vertexPaths){
                if(path.empty()) continue;
                /* edges.push_back(mesh->edgeBetween(mesh->vert(path[0]), mesh->vert(path[1]))->index());*/
                for (int i = 0; i < path.size() - 1; i++)
                {
                    edges.push_back(mesh->edgeBetween(mesh->vert(path[i]), mesh->vert(path[i+1]))->index() );
                    vertices.push_back(path[i]);
                }
            }
            for (auto a : landmarks)
                vertices.push_back(a);
        }
        std::ofstream _out(argv[3]);
        _out << "VERTICES\n";
        for (auto a : vertices){
            _out << a << std::endl;
        }

        _out << "EDGES\n";
        for (auto a : edges){
            _out << a << std::endl;
        }
        _out.close();


    return 1;
}
