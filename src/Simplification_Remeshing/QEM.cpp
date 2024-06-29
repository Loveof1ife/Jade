#include "PolyMesh/IOManager.h"
#include "Eigen/Dense"
#include <queue>
#include <map>
#include <ctime>
#include "memory"

#define lambda 1e3

using namespace acamcad;
using namespace polymesh;

struct edgePriority{
    double cost{};
    MHalfedge* halfEdge{};
    int state{};
    MPoint3 newPoint;
};

struct cmp{
    bool operator() (const edgePriority& a, const edgePriority& b){
        return a.cost > b.cost;
    }
};

void initialQ(std::map<MVert*, Eigen::Matrix4d>& qv, const std::unique_ptr<PolyMesh>& mesh);
void computeOptimalPositionAndEdgeError(std::priority_queue<edgePriority, std::vector<edgePriority>, cmp>& queueCost,
                                         std::map<MVert*, Eigen::Matrix4d>& qv, std::map<MHalfedge*, int>& state,
                                         const std::unique_ptr<PolyMesh>& mesh);

void calEdgeCost(MHalfedge* const & he, std::map<MVert*, Eigen::Matrix4d>& qv, std::map<MHalfedge*, int>& state,
                 std::priority_queue<edgePriority, std::vector<edgePriority>, cmp>& queueCost);

bool qemCollapse(const edgePriority& collapseEdge, const std::unique_ptr<PolyMesh>& mesh,
                 std::map<MVert*, Eigen::Matrix4d>& qv, std::map<MHalfedge*, int>& state,
                 std::priority_queue<edgePriority, std::vector<edgePriority>, cmp>& queueCost);

void updateQandCost(MVert* v_collapse, const std::unique_ptr<PolyMesh>& mesh,
                    std::map<MVert*, Eigen::Matrix4d>& qv, std::map<MHalfedge*, int>& state,
                    std::priority_queue<edgePriority, std::vector<edgePriority>, cmp>& queueCost);

void updateVertQ( MVert* const &v,  std::map<MVert*, Eigen::Matrix4d>& qv, const std::unique_ptr<PolyMesh>& mesh);

void QEM(const std::unique_ptr<PolyMesh>& mesh);

int main(int argc, const char **argv)
{
    if (argc != 3)
    {
        std::cout << "========== Hw9 Usage  ==========\n";
        std::cout << std::endl;
        std::cout << "Input:	ACAM_mesh_HW9.exe	input_mesh.obj	output_mesh.obj\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }

    std::string mesh_path = argv[1];
    auto mesh = std::make_unique<PolyMesh>();
    loadMesh(mesh_path, mesh.get());

    std::string out_path = argv[2];

    // read input mesh
    //PolyMesh* mesh = new PolyMesh();
    //loadMesh("cat_open.obj", mesh);

    clock_t start, end;
    std::cout << "Begin QEM" << std::endl;
    start = clock();
    QEM(mesh);
    end = clock();
    std::cout << "time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
    writeMesh(out_path, mesh.get());
}

void QEM(const std::unique_ptr<PolyMesh>& mesh){
    //initial Q of vertex , edge
    std::priority_queue<edgePriority, std::vector<edgePriority>, cmp> queueCost;
    std::map<MVert*, Eigen::Matrix4d> qv;
    std::map<MHalfedge*, int> state;

    initialQ(qv, mesh);
    //initial cost of contraction
    computeOptimalPositionAndEdgeError(queueCost, qv, state, mesh);

    //simplification
    int nv = mesh->numVertices();
    int min_nv = std::min((int)(0.5*nv), 1000);

    while(nv > min_nv) {
        edgePriority collapseEdge = queueCost.top();
        queueCost.pop();

        if (collapseEdge.state == state[collapseEdge.halfEdge]) {
            if (qemCollapse(collapseEdge, mesh, qv, state, queueCost)) {
                nv--;
            }
        }
    }
}

void initialQ(std::map<MVert*, Eigen::Matrix4d>& qv, const std::unique_ptr<PolyMesh>& mesh) {
    if (mesh->numVertices() == 3)
        return;

    for (MVert *v: mesh->vertices()) {
        MVector3 vPos = v->position();
        Eigen::Matrix4d tmpQ;
        double a, b, c, d;

        if (!mesh->isBoundary(v)) {
            for (auto voh_it = mesh->voh_iter(v); voh_it.isValid(); ++voh_it) {
                MHalfedge *voh = *voh_it;
                if (!voh->isBoundary()) {
                    MPolyFace *vohFace = voh->polygon();
                    MVector3 normal = vohFace->normal();
                    a = normal[0];
                    b = normal[1];
                    c = normal[2];
                    d = normal.dot(vPos);
                    Eigen::Matrix<double, 4, 1> qi{a, b, c, d};
                    tmpQ += qi * qi.transpose();
                }
            }
        }
        else
        {
            for(auto voh_it = mesh->voh_iter(v); voh_it.isValid(); ++voh_it){
                MHalfedge* voh = *voh_it;
                if(! voh->isBoundary() ){
                    MPolyFace *vohFace = voh->polygon();
                    MVector3 normal = vohFace->normal();
                    a = normal[0];
                    b = normal[1];
                    c = normal[2];
                    d = normal.dot(vPos);
                    Eigen::Matrix<double, 4, 1> qi{a, b, c, d};
                    tmpQ += qi * qi.transpose();
                }
                else{
                    MHalfedge* prev_voh = voh->prev();
                    MPolyFace* face_hh = voh->pair()->polygon(), *face_prev_hh = prev_voh->pair()->polygon();
                    MVector3 face_hh_nor = face_hh->normal(), face_prev_hh_nor = face_prev_hh->normal();
                    MVector3 vir_face_hh_nor = cross(voh->tangent(), face_hh_nor).normalized(), vir_face_prev_hh_nor = cross(prev_voh->tangent(), face_prev_hh_nor).normalized();
                    a = vir_face_hh_nor[0], b = vir_face_hh_nor[1], c = vir_face_hh_nor[2], d = -dot(vir_face_hh_nor, vPos);
                    Eigen::Matrix<double, 4, 1> p = { a,b,c,d };
                    tmpQ += lambda * (p * p.transpose());
                    a = vir_face_prev_hh_nor[0], b = vir_face_prev_hh_nor[1], c = vir_face_prev_hh_nor[2], d = -dot(vir_face_prev_hh_nor, vPos);
                    p = { a,b,c,d };
                    tmpQ  += lambda * (p * p.transpose());
                }
            }
        }
        qv.insert(std::make_pair(v, tmpQ));
    }
}
void computeOptimalPositionAndEdgeError(std::priority_queue<edgePriority, std::vector<edgePriority>, cmp>& queueCost,
                                 std::map<MVert*, Eigen::Matrix4d>& qv,
                                 std::map<MHalfedge*, int>& state,
                                 const std::unique_ptr<PolyMesh>& mesh)
{
    for(auto he : mesh->halfEdges()) {
        state.insert(std::make_pair(he, 0));
        calEdgeCost(he, qv, state, queueCost);
    }

}


bool qemCollapse(const edgePriority& collapseEdge, const std::unique_ptr<PolyMesh>& mesh,
                 std::map<MVert*, Eigen::Matrix4d>& qv, std::map<MHalfedge*, int>& state,
                 std::priority_queue<edgePriority, std::vector<edgePriority>, cmp>& queueCost)
{
    MHalfedge* he = collapseEdge.halfEdge;
    MHalfedge* he_oppo = he->pair();
    MVert* v0 = he->fromVertex(), * v1 = he->toVertex();
    MVert* v_collapse;
    bool is_collapse = false;
    if(mesh->is_collapse_ok(he)){
        v1->setPosition(collapseEdge.newPoint);
        mesh->collapse(he);
        is_collapse = true;
        v_collapse = v1;
    }
    else if(mesh->is_collapse_ok(he_oppo) ){
        v0->setPosition(collapseEdge.newPoint);
        mesh->collapse(he);
        is_collapse = true;
        v_collapse = v1;
    }
    if(is_collapse){
        updateQandCost(v_collapse, mesh, qv, state, queueCost);
    }
    return is_collapse;
}

void updateQandCost(MVert* v_collapse, const std::unique_ptr<PolyMesh>& mesh,
                    std::map<MVert*, Eigen::Matrix4d>& qv, std::map<MHalfedge*, int>& state,
                    std::priority_queue<edgePriority, std::vector<edgePriority>, cmp>& queueCost)
{
    //update v_collapse s qv
    updateVertQ(v_collapse, qv, mesh);
    //update the one_ring of v_collapse
    for(auto vv_it = mesh->vv_iter(v_collapse); vv_it.isValid(); ++vv_it){
        updateVertQ(*vv_it, qv, mesh);
    }
    for (auto voh_it = mesh->voh_iter(v_collapse); voh_it.isValid(); ++voh_it){
        auto he = *voh_it;
        state[he]++;
        calEdgeCost(he, qv, state, queueCost);
    }
}

void updateVertQ( MVert* const &v,  std::map<MVert*, Eigen::Matrix4d>& qv, const std::unique_ptr<PolyMesh>& mesh){
    MVector3 vPos = v->position();
    Eigen::Matrix4d qv_new;
    double a, b, c, d;
    MVector3 face_normal;
    Eigen::Matrix<double, 4, 1> plane;

    if (!mesh->isBoundary(v) ){
        for(auto vf_it = mesh->vf_iter(v); vf_it.isValid(); ++vf_it) {
            face_normal = (*vf_it)->normal();
            a = face_normal[0], b = face_normal[1], c = face_normal[2], d = -dot(face_normal, vPos);
            plane = {a, b, c, d};
            qv_new += plane * plane.transpose();
        }
        qv[v] = qv_new;
    }
    else{
        for(auto voh_it = mesh->voh_iter(v); voh_it.isValid(); ++voh_it){
            MHalfedge* voh = *voh_it;
            if(! voh->isBoundary() ){
                MPolyFace *vohFace = voh->polygon();
                MVector3 normal = vohFace->normal();
                a = normal[0];
                b = normal[1];
                c = normal[2];
                d = normal.dot(vPos);
                Eigen::Matrix<double, 4, 1> plane{a, b, c, d};
                qv_new += plane * plane.transpose();
            }
            else{
                MHalfedge* prev_voh = voh->prev();
                MPolyFace* face_hh = voh->pair()->polygon(), *face_prev_hh = prev_voh->pair()->polygon();
                MVector3 face_hh_nor = face_hh->normal(), face_prev_hh_nor = face_prev_hh->normal();
                MVector3 vir_face_hh_nor = cross(voh->tangent(), face_hh_nor).normalized(), vir_face_prev_hh_nor = cross(prev_voh->tangent(), face_prev_hh_nor).normalized();
                a = vir_face_hh_nor[0], b = vir_face_hh_nor[1], c = vir_face_hh_nor[2], d = -dot(vir_face_hh_nor, vPos);
                Eigen::Matrix<double, 4, 1> p = { a,b,c,d };
                qv_new += lambda * (p * p.transpose());
                a = vir_face_prev_hh_nor[0], b = vir_face_prev_hh_nor[1], c = vir_face_prev_hh_nor[2], d = -dot(vir_face_prev_hh_nor, vPos);
                p = { a,b,c,d };
                qv_new  += lambda * (p * p.transpose());
            }
        }
    }
    qv[v] = qv_new;
}

void calEdgeCost(MHalfedge* const & he, std::map<MVert*, Eigen::Matrix4d>& qv, std::map<MHalfedge*, int>& state,
                 std::priority_queue<edgePriority, std::vector<edgePriority>, cmp>& queueCost)
{
    MVert* v0{he->fromVertex()}, * v1{he->toVertex()};
    Eigen::Matrix4d qe = qv[v0] + qv[v1], qSolver = qe;
    qSolver(3, 0) = 0.0, qSolver(3, 1) = 0.0, qSolver(3, 2) = 0.0, qSolver(3, 3) = 1.0;

    MPoint3 newPoint;
    Eigen::Vector4d newVec;

    if (qSolver.determinant() == 0){
        MVector3 temp = 0.5 * (v0->position() + v1->position());
        newPoint = { temp[0], temp[1], temp[2] };
        newVec = {newPoint[0], newPoint[1], newPoint[2], 1.0 };
    }
    else{
        Eigen::Vector4d b = { 0.0,0.0,0.0,1.0 };
        newVec = qSolver.inverse()*b;
        newPoint = { newVec[0], newVec[1], newVec[2] };
    }

    edgePriority temp;
    temp.halfEdge = he;
    temp.cost = newVec.transpose() * qe * newVec;
    temp.newPoint = newPoint;
    temp.state = state[he];
    queueCost.push(temp);
}