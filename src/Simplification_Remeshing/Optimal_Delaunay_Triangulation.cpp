#include "PolyMesh/IOManager.h"
#include <string>
#include <memory>

constexpr auto M_PI = 3.14159265358979323846;

using namespace acamcad;
using namespace polymesh;

double  get_triface_area(MPolyFace* f)
{
    MHalfedge* he = f->halfEdge();
    MPoint3 v0 = he->fromVertex()->position();
    MPoint3 v1 = he->toVertex()->position();
    MPoint3 v2 = he->next()->toVertex()->position();
    return 0.5 * ((v1 - v0) % (v2 - v0)).norm();
}

MPoint3 get_triface_circumcenter(MPolyFace* f) {
    MHalfedge *he = f->halfEdge();
    MPoint3 v0 = he->fromVertex()->position();
    MPoint3 v1 = he->toVertex()->position();
    MPoint3 v2 = he->next()->toVertex()->position();

    double x1, y1, x2, y2, x3, y3;
    x1 = v0[0];
    y1 = v0[1];
    x2 = v1[0];
    y2 = v1[1];
    x3 = v2[0];
    y3 = v2[1];

    double a1, b1, c1, a2, b2, c2;
    a1 = 2 * (x2 - x1);
    a2 = 2 * (x3 - x2);
    c1 = x2 * x2 + y2 * y2 - x1 * x1 - y1 * y1;
    b1 = 2 * (y2 - y1);
    b2 = 2 * (y3 - y2);
    c2 = x3 * x3 + y3 * y3 - x2 * x2 - y2 * y2;

    MPoint3 circumcenter(0.0, 0.0, 0.0);
    circumcenter[0] = (b2 * c1 - b1 * c2) / (a1 * b2 - a2 * b1);
    circumcenter[1] = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);
    circumcenter[2] = 0;

    return circumcenter;
}
void optimal_delaunay_trianglation(int iter_num, std::unique_ptr<PolyMesh>const & mesh){
    for(int iter = 0; iter < iter_num; iter++){
        //update triangulations to build DT: empty circular property
        for(auto e_it = mesh->halfedge_begin(); e_it != mesh->halfedge_end(); e_it++){
            MHalfedge* edge = *e_it;
            if(edge->isBoundary()) continue;

            MHalfedge* next_edge = edge->next(), * pair_edge = edge->pair();
            MHalfedge* pair_next_edge = pair_edge->next();

            MVert* v0 = edge->fromVertex(), * v1 = edge->toVertex();
            MVert* v2 = next_edge->toVertex();
            MVert* v3 = pair_next_edge->toVertex();

            double alpha(0.0), alpha1(0.0), alpha2(0.0);
            alpha1 = acos((pow((v0->position() - v2->position()).norm(), 2) + pow((v1->position() - v2->position()).norm(), 2) - pow((v0->position() - v1->position()).norm(), 2))
                             /(2 * (v0->position() - v2->position()).norm() * ( v1->position() - v2->position()).norm()));
            alpha2 = acos((pow((v0->position() - v3->position()).norm(), 2) + pow((v1->position() - v3->position()).norm(), 2) - pow((v0->position() - v1->position()).norm(), 2))
                             /(2 * (v0->position() - v3->position()).norm() * (v1->position() - v3->position()).norm()));

            alpha = alpha1 + alpha2;
            if (alpha > M_PI) {
                mesh->flipEdgeTriangle(edge->edge());
            }
        }
            //update vertices
        for(auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++){
            MVert* vert = *v_it;
            MPoint3 center{0.0, 0.0, 0.0};
            double area{0.0};
            if(mesh->isBoundary(vert)) continue;
            for(auto vf_it = mesh->vf_iter(vert); vf_it.isValid(); vf_it++){
                MPolyFace* one_ring_face = *vf_it;
                center += get_triface_circumcenter(one_ring_face);
                area += get_triface_area(one_ring_face);
            }
        }
    }
}
int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cout << "========== Hw12 Usage  ==========\n";
        std::cout << std::endl;
        std::cout << "Input:	ACAM_mesh_HW12.exe	input_mesh.obj	output_mesh.obj\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }

    auto mesh = std::make_unique<PolyMesh>();
    //∂¡»ÎÕ¯∏Ò
    std::string mesh_path = argv[1];
    loadMesh(mesh_path, mesh.get());
    std::string out_path = argv[2];

    //loadMesh("rectangle.obj", mesh);
    int iter_num = 1000;	// iterative number
    optimal_delaunay_trianglation(iter_num, mesh);
    writeMesh(out_path, mesh.get());
}