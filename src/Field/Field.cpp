#include <PolyMesh/IOManager.h>
#include<complex>
#include<Eigen/Eigen>
#include<fstream>
#include<cmath>
#include<string>
#include<memory>

using namespace acamcad;
using namespace polymesh;

#define  PI 3.14159265358979323846

void calculateMeshFaceBase(std::unique_ptr<PolyMesh>& mesh, std::vector<MVector3>& f_base)
{
    f_base.resize(mesh->numPolygons() * 2);
    for(auto face_it = mesh->polyfaces_begin(); face_it != mesh->polyfaces_end(); face_it++){
        MPolyFace* face = *face_it;
        assert(face->PolyNum() >= 3);

        int face_id = face->index();
        auto fv_it = mesh->fv_iter(face);

        MVert* v1, * v2, * v3;
        MVector3 p1, p2, p3;
        v1 = *fv_it; p1 = v1->position(); ++fv_it;
        v2 = *fv_it; p2 = v2->position(); ++fv_it;
        v3 = *fv_it; p3 = v3->position();

        f_base[face_id * 2] = (p2 - p1).normalized();
        MVector3 normal = cross((p3 - p1), f_base[face_id * 2]);
        normal.normalize();
        f_base[face_id * 2 + 1] = cross(normal, f_base[face_id * 2]).normalized();
    }
}

void crossfieldCreator(std::unique_ptr<PolyMesh>& mesh, std::vector<int>& cons_id, std::vector<MVector3>& cons_vec, std::vector<MVector3>& crossfield)
{
    using namespace std;
    using COMPLEX = complex<double> ;

    size_t f_num = mesh->numPolygons();
    vector<int> status(f_num, 0);
    vector<COMPLEX> f_direction(f_num);

    crossfield.clear();
    crossfield.resize(f_num);
    vector<MVector3> f_base(f_num * 2);

    calculateMeshFaceBase(mesh, f_base);

    for(int i = 0; i < cons_id.size(); i++)
    {
        int constrain_fid = cons_id[i];
        status[constrain_fid] = 1;
        MVector3 constrain_vec = cons_vec[i];
        f_direction[constrain_fid] = std::pow( COMPLEX(dot(f_base[2 * constrain_fid], constrain_vec), dot(f_base[2 * constrain_fid + 1], constrain_vec)), 4);
    }

    vector<int> id2sln(f_num, -1);
    vector<int> sln2id(0);
    int count = 0;

    for (int i = 0; i < f_num; i++) {
        if(status[i] == 0){
            sln2id.push_back(i); //row_id map back to face_id
            id2sln[i] = count; // face_id map to the row_id of X (AX = b)
            count++;
        }
    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> slu;
    Eigen::SparseMatrix<COMPLEX> A;
    Eigen::VectorXcd b_pre(mesh->numEdges());
    Eigen::VectorXcd b;
    b_pre.setZero();
    vector<Eigen::Triplet<COMPLEX>> tris;

    count = 0;
    for(auto f_it = mesh->polyfaces_begin(), f_end = mesh->polyfaces_end(); f_it != f_end; ++f_it){
        MPolyFace* face = *f_it;
        size_t f_id = face->index();

        for(auto fh_it = mesh->fhe_iter(face); fh_it.isValid(); fh_it++){
            auto f_halfedge = *fh_it;
            if(!mesh->isBoundary(f_halfedge)){
                MVector3 p1, p2;
                auto fo_halfedge = f_halfedge->pair();
                MPolyFace* face_oppo = fo_halfedge->polygon();
                size_t g_id = face_oppo->index();
                p1 = f_halfedge->toVertex()->position();
                p2 = f_halfedge->fromVertex()->position();

                MVector3 e = (p2 - p1).normalized();

                COMPLEX e_f = COMPLEX(dot(e, f_base[f_id * 2]), dot(e, f_base[f_id * 2 + 1]));
                COMPLEX e_g = COMPLEX(dot(e, f_base[g_id* 2]), dot(e, f_base[g_id * 2 + 1]));

                COMPLEX e_f_c_4 = pow(conj(e_f), 4);
                COMPLEX e_g_c_4 = pow(conj(e_g), 4);

                if(f_id < g_id){
                    if (status[f_id] == 1 && status[g_id] == 1) continue;
                    if (status[f_id] == 0){
                        tris.push_back((Eigen::Triplet<COMPLEX>(count, id2sln[f_id], e_f_c_4) ) );
                    }
                    else
                    {
                        b_pre[count] += -e_f_c_4 * f_direction[f_id];
                    }
                    if (status[g_id] == 0)
                    {
                        tris.push_back(Eigen::Triplet<COMPLEX>(count, id2sln[g_id], -e_g_c_4));
                    }
                    else
                    {
                        b_pre[count] += e_g_c_4 * f_direction[g_id];
                    }
                    count++;
                }
            }
        }
    }
    // no_boundary_edge * valid_face
    A.resize(count, sln2id.size() );
    b.resize(count);
    b = b_pre.head(count);
    A.setFromTriplets(tris.begin(), tris.end());
    Eigen::SparseMatrix<COMPLEX> AT = A.adjoint();
    slu.compute(AT * A);
    Eigen::VectorXcd x = slu.solve(AT * b);

    crossfield.resize(4 * f_num);
    for (int i = 0; i < f_num; i++)
    {
        if (status[i] == 0)
        {
            f_direction[i] = x(id2sln[i]);
        }
        double length = 1;
        double arg = std::arg(f_direction[i]) / 4;
        for (int j = 0; j < 4; j++)
        {
            crossfield[i * 4 + j] = f_base[i * 2] * length * cos(arg + j * PI / 2) + f_base[i * 2 + 1] * length * sin(arg + j * PI / 2);
        }
    }
}

void loadConstrains(const std::string& filename, std::vector<int>& cons_id, std::vector<MVector3>& cons_vec)
{
    std::fstream infile(filename.c_str(), std::ios_base::in);
    int num;
    infile >> num;
    cons_id.reserve(num);
    cons_vec.reserve(num);
    int tempi;
    for (int i = 0; i < num; i++)
    {
        infile >> tempi;
        cons_id.push_back(tempi);
        double p0, p1, p2;
        infile >> p0;	infile >> p1;	infile >> p2;
        cons_vec.emplace_back(p0, p1, p2);
    }
    infile.close();
}

void writeCrossField(const std::string& filename, std::vector<MVector3>& crossfield)
{
    std::fstream ofile(filename.c_str(), std::ios_base::out);
    int num = crossfield.size();
    ofile << num << std::endl;
    for (int i = 0; i < num; i++)
    {
        ofile << crossfield[i][0] << " " << crossfield[i][1] << " " << crossfield[i][2] << std::endl;
    }
    ofile.close();
}

int main(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cout << "========== Hw10 Usage  ==========\n";
        std::cout << std::endl;
        std::cout << "Input:	ACAM_mesh_HW10.exe	mesh.obj	constrains.txt	crossfield.txt\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }

    //¶ÁÈëÍø¸ñ
    std::string mesh_path = argv[1];
    auto mesh = std::make_unique<PolyMesh>();
    loadMesh(mesh_path, mesh.get());

    std::string input_constrains = argv[2];
    std::string output_crossfield = argv[3];

    //std::string input_mesh = "input_mesh.obj";
    //std::string input_constrains = "constrains.txt";
    //std::string output_crossfield = "crossfield.txt";

    //PolyMesh* mesh = new PolyMesh();
    std::vector<int> cons_id(0);
    std::vector<MVector3> cons_vec(0), crossfield(0);
    //loadMesh(input_mesh, mesh);
    loadConstrains(input_constrains, cons_id, cons_vec);
    crossfieldCreator(mesh,cons_id, cons_vec, crossfield);
    writeCrossField(output_crossfield, crossfield);
    return 1;
}