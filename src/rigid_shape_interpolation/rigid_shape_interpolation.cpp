#include "PolyMesh/IOManager.h"
#include "memory"
#include <Eigen/Sparse>
#include <Eigen/SVD>

using namespace acamcad;
using namespace polymesh;
using namespace Eigen;

void MeshInterpolation(const std::shared_ptr<PolyMesh>& source_mesh, const std::shared_ptr<PolyMesh>& target_mesh, double delta_t );

void ComputeTransferMatrix(const std::shared_ptr<PolyMesh>& source_mesh, const std::shared_ptr<PolyMesh>& target_mesh,
                           std::vector<Matrix2d>& S, std::vector<double>& angle, std::vector<double>& area, size_t nv, size_t nf);

std::unique_ptr<PolyMesh>  SolvePosition(const std::shared_ptr<PolyMesh>& source_mesh, const std::shared_ptr<PolyMesh>& target_mesh,
                       const std::vector<Matrix2d>& interpolation_A);

int main(int argc, const char **argv)
{
    if (argc != 4)
    {
        std::cout << "========== Hw8 Usage  ==========\n";
        std::cout << std::endl;
        std::cout << "blablabla.exe	source.obj target.obj delta_t\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }
    std::string source_path = argv[1];
    std::string target_path = argv[2];
    std::string t_string = argv[3];

    auto source_mesh {std::make_shared<PolyMesh>()};
    auto target_mesh {std::make_shared<PolyMesh>()};
    loadMesh(source_path, source_mesh.get());
    loadMesh(target_path, target_mesh.get());

    char* end;
    double delta_t = strtod(t_string.c_str(), &end);
    MeshInterpolation(source_mesh, target_mesh, delta_t);
    return 1;
}

void MeshInterpolation(const std::shared_ptr<PolyMesh>& source_mesh, const std::shared_ptr<PolyMesh>& target_mesh, double delta_t ){
    //for each triangle T compute J from src to tar
    //interpolate A(t)
    //reconstruct mesh(t)
    size_t nv = source_mesh->numVertices();
    size_t nf = source_mesh->numPolygons();

    std::vector<Matrix2d> S(nf);
    std::vector<double> area(nf);
    std::vector<double> angle(nf);

    ComputeTransferMatrix(source_mesh, target_mesh, S, angle, area, nv, nf);

    Matrix2d i = MatrixXd::Identity(2, 2);
    std::vector<Matrix2d>interpolation_A{};

    double interpolation_angle{};
    Matrix2d interRotation{};

    for(auto face_it = source_mesh->polyfaces_begin(); face_it != source_mesh->polyfaces_end(); face_it++){
        auto current_face = (*face_it);
        int face_index = current_face->index();

        interpolation_angle = (1 - delta_t) * angle[face_index];

        interRotation<< cos(interpolation_angle ), -sin(interpolation_angle ),
                        sin(interpolation_angle), cos(interpolation_angle );

        interpolation_A[face_index] =  interRotation * (delta_t * i + (1 - delta_t) * S[face_index]);
    }
    std::unique_ptr<PolyMesh> interPolyMesh{SolvePosition(source_mesh, target_mesh, interpolation_A)};
}


void ComputeTransferMatrix(const std::shared_ptr<PolyMesh>& source_mesh, const std::shared_ptr<PolyMesh>& target_mesh,
                           std::vector<Matrix2d>& S, std::vector<double>& angle, std::vector<double>& area, size_t nv, size_t nf)
{
    std::vector<Vector3d> srcX(nv), srcY(nv);
    std::vector<std::vector<int> > v_id(nf);

    for(auto face_it = source_mesh->polyfaces_begin(); face_it != source_mesh->polyfaces_end(); face_it++) {
        auto src_face = (*face_it);
        int face_index = src_face->index();
        MHalfedge *he = src_face->halfEdge();

        MVert *v0 = he->fromVertex();
        MVert *v1 = he->toVertex();
        MHalfedge *next_he = he->next();
        MVert *v2 = next_he->toVertex();

        v_id[face_index].push_back(v0->index());
        v_id[face_index].push_back(v1->index());
        v_id[face_index].push_back(v2->index());
        area[face_index] = cross(v1->position() - v0->position(), v2->position() - v0->position()).norm() / 2.0;

        srcX[face_index] = Vector3d(v2->position()[0] - v1->position()[0], v0->position()[0] - v2->position()[0],
                                    v1->position()[0] - v0->position()[0]);
        srcY[face_index] = Vector3d(v2->position()[1] - v1->position()[1], v0->position()[1] - v2->position()[1],
                                    v1->position()[1] - v0->position()[1]);

        Vector3d tarX(target_mesh->vert(v0->index())->position()[0], target_mesh->vert(v1->index())->position()[0],
                      target_mesh->vert(v2->index())->position()[0]);
        Vector3d tarY(target_mesh->vert(v0->index())->position()[1], target_mesh->vert(v1->index())->position()[1],
                      target_mesh->vert(v2->index())->position()[1]);

        Matrix2d Ai;
        Ai << srcY[face_index].dot(tarX) / (2 * area[face_index]), srcX[face_index].dot(tarX) / (2 * area[face_index]),
                srcY[face_index].dot(tarY) / (2 * area[face_index]), srcX[face_index].dot(tarY) /
                                                                     (2 * area[face_index]);
        Eigen::JacobiSVD<Matrix2d> svdAi(Ai, ComputeThinU | ComputeFullV);

        const Matrix2d& u = svdAi.matrixU();
        const Matrix2d& v = svdAi.matrixV();
        Matrix2d rotation = u * v.transpose();

        Matrix2d sigma;
        sigma << svdAi.singularValues()[0], 0,
                0, svdAi.singularValues()[1];
        S[face_index] = sigma;
        angle[face_index] = atan2(rotation(1, 0), rotation(1, 1));
    }
}

std::unique_ptr<PolyMesh> SolvePosition(const std::shared_ptr<PolyMesh>& source_mesh, const std::shared_ptr<PolyMesh>& target_mesh,
                       const std::vector<Matrix2d>& interpolation_A)
{
    std::unique_ptr<PolyMesh> interPolyMesh{std::make_unique<PolyMesh>()};



    return interPolyMesh;
}


