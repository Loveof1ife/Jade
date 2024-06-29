#include <iostream>
#include <cstring>
#include <memory>
#include "AABB_Tree/AABB_Tree.h"
#include "PolyMesh/IOManager.h"

using namespace acamcad;
using namespace polymesh;

double calculateTargetEdgeLength(std::unique_ptr<PolyMesh>& mesh);
std::unique_ptr<AABB_Tree> get_AABB_tree();

void split_long_edges(double high, std::unique_ptr<PolyMesh>& mesh);
void collapse_short_edges(double high, double low, std::unique_ptr<PolyMesh>& mesh);
void tangential_relaxation(std::unique_ptr<PolyMesh>& mesh);
void equalize_valences(std::unique_ptr<PolyMesh>& mesh);

void project_to_surface( std::unique_ptr<AABB_Tree>& abtree, std::unique_ptr<PolyMesh>& mesh );
int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cout << "========== Hw11 Usage  ==========\n";
        std::cout << std::endl;
        std::cout << "Input:	ACAM_mesh_HW11.exe	input_mesh.obj	output_mesh.obj\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }

    auto mesh = std::make_unique<PolyMesh>();
    //∂¡»ÎÕ¯∏Ò
    std::string mesh_path = argv[1];
    loadMesh(mesh_path, mesh.get());

    std::string out_path = argv[2];

    //mesh load an write , now only support obj/off
    //loadMesh("small_bunny.obj", mesh);
    double target_edge_length, high, low;
    target_edge_length = calculateTargetEdgeLength(mesh)/2.0;
    high = 4.0 / 3.0 * target_edge_length;
    low = 4.0 / 5.0 * target_edge_length;
    std::unique_ptr<AABB_Tree> abtree {get_AABB_tree()};
    for (int i = 0; i < 10; i++)
    {
        split_long_edges(high, mesh);
        collapse_short_edges(high, low, mesh);
        equalize_valences(mesh);
        tangential_relaxation(mesh);
        project_to_surface(abtree, mesh);
    }
    writeMesh(out_path, mesh.get());
}


double calculateTargetEdgeLength(std::unique_ptr<PolyMesh>& mesh){
    double target_edge_length = {0.0};
    for(auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++){
        target_edge_length += (*e_it)->length();
    }
    target_edge_length /= mesh->numEdges();
    return target_edge_length;
}

std::unique_ptr<AABB_Tree> get_AABB_tree(std::unique_ptr<PolyMesh>& mesh){
    std::vector<Vector3f> point_set;
    point_set.clear();
    for(auto f_it = mesh->polyfaces_begin(); f_it != mesh->polyfaces_end(); ++f_it){
        MPolyFace* face = *f_it;
        for(auto fv_it = mesh->fv_iter(face); fv_it.isValid(); ++fv_it){
            MVert* fv = *fv_it;
            Vector3f p{static_cast<float>(fv->x() ),
                       static_cast<float>(fv->y() ),
                       static_cast<float>(fv->z() )};
            point_set.push_back(p);
        }
    }
    return std::make_unique<AABB_Tree>(point_set);
}

void split_long_edges(double high, std::unique_ptr<PolyMesh>& mesh){
    for(auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++){
        MEdge* e = *e_it;
        double len = e->length();
        if (len > high)	mesh->splitEdgeTriangle(e);
    }
}

void collapse_short_edges(double high, double low, std::unique_ptr<PolyMesh>& mesh){
    for(auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++){
        MEdge* e = *e_it;
        MHalfedge* he = e->halfEdge();
        MVert* p0 = he->fromVertex();
        MVert* p1 = he->toVertex();
        if (!mesh->is_collapse_ok(he)) continue;
        if (mesh->isBoundary(p0) || mesh->isBoundary(p1)) continue;
        double len = e->length();
        if (len < low){
            bool is_collapse = true;
            for (auto vv_it = mesh->vv_iter(p0); vv_it.isValid(); ++vv_it)
            {
                MVert* vert = *vv_it;
                double length= (p1->position() - vert->position()).norm();
                if (length> high)
                {
                    is_collapse = false;
                    break;
                }
            }
            if(is_collapse) mesh->collapseTriangle(he);
        }
    }

}
void equalize_valences(std::unique_ptr<PolyMesh>& mesh){
    std::vector<int> target_valence;
    int deviation_pre, deviation_post;
    for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
    {
        if (mesh->isBoundary(*v_it))	target_valence.push_back(4);
        else target_valence.push_back(6);
    }
    for(auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it){
        if (mesh->isBoundary(*e_it) || !mesh->is_flip_ok_Triangle(*e_it)) continue;
        MHalfedge* he1 = (*e_it)->halfEdge();
        MVert* v1 = he1->fromVertex();
        MVert* v2 = he1->toVertex();

        MHalfedge* he2 = (*e_it)->halfEdge()->next();
        MVert* v3 = he2->toVertex();
        MHalfedge* he3 = (*e_it)->halfEdge()->pair()->next();
        MVert* v4 = he3->toVertex();

        deviation_pre =   std::abs( static_cast<int>(mesh->valence(v1) ) - target_valence[v1->index()])
                        + std::abs( static_cast<int>(mesh->valence(v2) ) - target_valence[v2->index()])
                        + std::abs( static_cast<int>(mesh->valence(v3) ) - target_valence[v3->index()])
                        + std::abs( static_cast<int>( mesh->valence(v4) ) - target_valence[v4->index()]);
        deviation_post =   std::abs( static_cast<int>(mesh->valence(v1) ) - 1 - target_valence[v1->index()])
                         + std::abs( static_cast<int>(mesh->valence(v2) ) - 1  - target_valence[v2->index()])
                         + std::abs( static_cast<int>(mesh->valence(v3) ) + 1 - target_valence[v3->index()])
                         + std::abs( static_cast<int>( mesh->valence(v4) ) + 1 - target_valence[v4->index()]);
        if (deviation_pre> deviation_post) mesh->flipEdgeTriangle(*e_it);
}

void tangential_relaxation(std::unique_ptr<PolyMesh>& mesh){
    for(auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++){
        MVert* vert = *v_it;
        if (mesh->isBoundary(vert)) continue;
        MVector3 one_ring_average{0, 0, 0};
        MVector3 normal = vert->normal();
        int count{0};
        for(auto vv_it = mesh->vv_iter(vert); vv_it.isValid(); ++vv_it){
            MVert* one_ring = *vv_it;
            one_ring_average += one_ring->position();
            count++;
        }
        one_ring_average /= count;
        vert->setPosition( (one_ring_average - normal.dot(one_ring_average.point()- vert->position()) * normal).point());
    }
}

void project_to_surface( std::unique_ptr<AABB_Tree>& abtree, std::unique_ptr<PolyMesh>& mesh ){
    for(auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it){
        Vector3f p;
        p[0] = float((*v_it)->x());
        p[1] = float((*v_it)->y());
        p[2] = float((*v_it)->z());
        Vector3f ab_nearst_point;
        abtree->findNearstPoint(p, ab_nearst_point);
        MPoint3 new_point;
        new_point[0] = double(ab_nearst_point[0]);
        new_point[1] = double(ab_nearst_point[1]);
        new_point[2] = double(ab_nearst_point[2]);
        (*v_it)->setPosition(new_point);
    }
}

