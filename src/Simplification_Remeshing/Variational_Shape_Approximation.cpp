#include "PolyMesh/IOManager.h"
#include <fstream>
#include <queue>
#include <sstream>
#include <ctime>
#include <memory>

using namespace acamcad;
using namespace polymesh;

struct Tri {
    int id;
    int tag;
    double weight;

    Tri(int id, int tag, double weight) :tag(tag), id(id), weight(weight) {};

    bool operator<(const Tri& a) const {
        return weight > a.weight;
    }
};

double calculateArea(std::unique_ptr<PolyMesh> const &mesh, MHalfedge* he) {
    MVector3 p0 = he->fromVertex()->position();
    MVector3 p1 = he->toVertex()->position();
    MVector3 p2 = he->next()->toVertex()->position();
    return abs(norm(cross(p0 - p1, p2 - p1)));
};

void vsa_algorithm(std::unique_ptr<PolyMesh> const &mesh, int k_num,  //聚类数量
                   std::vector<int>& partition, //  face_index to proxy_index
                   std::vector<MVector3>& planes) // Proxy normal
{
    int f_num = mesh->numPolygons();
    partition.clear();
    partition.resize(f_num, -1);
    planes.resize(k_num);

    srand(static_cast<unsigned>(time(0) ) );
    mesh->updateMeshNormal();
    std::vector<bool> conquered(f_num, false);
    std::vector<int> proxy_seed(k_num);

    //Initial seeding
    for (int i = 0; i < k_num; ) {
        int random_face = rand() % f_num;
        if(partition[random_face] == -1 ){
            partition[random_face] = i; // random_face -> proxy_index_i
            planes[i] = mesh->polyface(random_face)->normal();// proxy_index_i normal = partition[random_face] normal
            proxy_seed[i] = random_face; //Store the index of the seed face for the current proxy
            i++;
        }
    }
    double last_energy {0.0};
    for (int iteration = 0; true; iteration++) {
        std::vector<double> smallest_energy(k_num, 1e10);
        //identifies the face with the smallest distortion error (the most similar face) for each proxy in the current partition
        for(auto f_it = mesh->polyfaces_begin(); f_it!=mesh->polyfaces_end(); f_it++){
            MPolyFace* face = *f_it;
            int face_id = face->index();
            if (partition[face_id] != -1){
                double energy = norm(planes[partition[face_id]] - mesh->polyface(face_id)->normal()); // Calculate the energy (distortion error)
                if (energy < smallest_energy[partition[face_id]]) { // Update the smallest energy and seed face for the proxy
                    smallest_energy[partition[face_id]] = energy;
                    proxy_seed[partition[face_id]] = face_id;
                }
            }
        }

        std::fill(conquered.begin(), conquered.end(), false); // Reset the conquered vector for the next iteration
        std::priority_queue<Tri> pqueue; // Priority queue for region growing

        // Initialize queue with adjacent faces of seeds
        for (int i = 0; i < k_num; i++){
            partition[proxy_seed[i]] = i; // Ensure the seed face is assigned to its proxy
            conquered[proxy_seed[i]] = true; // Mark the seed face as conquered
            for (FaceFaceIter ff_it = mesh->ff_iter(mesh->polyface(proxy_seed[i])); ff_it.isValid(); ff_it++) {
                MPolyFace* face = *ff_it;
                int face_id = face->index();
                double energy = norm(planes[i] - face->normal());
                pqueue.emplace(face_id, i, energy);
            }
        }

        while(!pqueue.empty()){
            int face_id = pqueue.top().id; // Get the face ID with the smallest energy
            int proxy_id = pqueue.top().tag; // Get the proxy ID associated with this face
            pqueue.pop(); // Remove the top element from the queue
            if(!conquered[face_id]){
                partition[face_id] = proxy_id;
                conquered[face_id] = true;
                for (FaceFaceIter ff_it = mesh->ff_iter(mesh->polyface(face_id) ); ff_it.isValid(); ff_it++) {
                    MPolyFace* adjacent_face = *ff_it;
                    int adjacent_face_id = adjacent_face->index();
                    if(!conquered[face_id]){
                        double energy = norm(planes[proxy_id] - adjacent_face->normal());
                        pqueue.emplace(adjacent_face_id, proxy_id, energy);
                    }
                }
            }
        }

        std::vector<MVector3>    fitting_normal(k_num, MVector3(0, 0, 0));
        for(auto f_it = mesh->polyfaces_begin(); f_it!=mesh->polyfaces_end(); f_it++) {
            auto face = *f_it;
            double area = calculateArea(mesh, face->halfEdge());
            fitting_normal[partition[face->index()] ] +=  face->normal() * area;
        }

        for (int i = 0; i < k_num; i++) {
            fitting_normal[i].normalized();
            planes[i] = fitting_normal[i];
        }
    }
}
void writePartition(const std::string& filename, std::vector<int>& partition)
{
    std::fstream ofile(filename.c_str(), std::ios_base::out);
    ofile << partition.size()<<std::endl;
    for (int i = 0; i < partition.size(); i++)
    {
        ofile << partition[i] << std::endl;
    }
}


int main(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cout << "========== Hw13 Usage  ==========\n";
        std::cout << std::endl;
        std::cout << "ACAM_mesh_HW13.exe [model] [partion_num] [out_partion]\n";
        std::cout << "ACAM_mesh_HW13.exe	mesh.obj 20 Partition.txt\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0 ;
    }
    //读网格
    std::string mesh_path = argv[1];
    std::stringstream ss;
    std::string cmd2 = argv[2];
    ss << cmd2;
    double K;
    ss >> K;
    std::string output_partition = argv[3];

    //std::string input_mesh = "input_mesh.obj";
    //double K = 20;
    //std::string output_partition = "Partition.txt";

    auto mesh = std::make_unique<PolyMesh>();
    std::vector<int> partition;
    std::vector<MVector3> plane(K);

    loadMesh(mesh_path, mesh.get());
    vsa_algorithm(mesh, K, partition, plane);
    writePartition(output_partition, partition);
}