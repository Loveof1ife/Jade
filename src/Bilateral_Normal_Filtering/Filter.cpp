#include "Filter.h"

double Gaussian(const MVector3& ni, const MVector3& nj)
{
    double delta_n = (ni - nj).norm();
    return exp(sqrt(delta_n) / 2 / sqrt(NormalSigma));
}
double Gaussian(const MPoint3& ci, const MPoint3& cj){
    double delta_c = (ci - cj).norm();
    return exp(sqrt(delta_c) / 2 / sqrt(CenterSigma));
}

void Bilateral_Normal_Filtering(const PolyMesh* in_mesh)
{
    auto mesh = new PolyMesh(*in_mesh);
    //计算面法向
    //法向滤波：面i的新的法向等于i与其周围1-邻域面的法向的加权（GS）平均
    for(auto face : mesh->polyfaces()){
        MVector3 ni = face->normal();
        MPoint3 ci = face->getFaceCenter();
        auto startEdge = face->halfEdge();
        auto currentEdge = startEdge;
        MVector3 totalNormal{0.0, 0.0, 0.0};

        double K{0.0};
        do{
            auto face_j = currentEdge->pair()->polygon();
            double Kj{0.0};

            if(face_j == nullptr) // boundary
                continue;
            MVector3 nj = face_j->normal();
            MPoint3 cj = face_j->getFaceCenter();

            MHalfedge* hedge1 = face_j->halfEdge();
            MHalfedge* hedge2 = hedge1->next();
            auto e1 = hedge1->toVertex()->position() - hedge1->fromVertex()->position();
            auto e2 = hedge2->toVertex()->position() - hedge2->fromVertex()->position();
            double Aj = e1.cross(e2).norm();

            Kj = Aj * Gaussian(ci, cj) * Gaussian(ni, nj);
            totalNormal += Kj * nj;
            K += Kj;
            currentEdge = currentEdge->next();
        }while(currentEdge != startEdge);
        totalNormal /= K;
        face->setNormal(totalNormal);
    }
    //根据滤波后的面法向重建PolyMesh

    for(size_t i = 0; i < IteratorTimes; i++)
    {

        for(auto vertex : mesh->vertices())
        {

            MHalfedge* startEdge = vertex->halfEdge();  // Get the starting half-edge
            MHalfedge* currentEdge = startEdge;
            MVector3 v_rec{0.0, 0.0, 0.0};
            size_t n{0};
            do{

                auto vertex_j = currentEdge->toVertex();
                MVector3 v_ij(vertex_j->position() - vertex->position());
                MVector3 n_ij = currentEdge->polygon()->normal();
                v_rec += n_ij * n_ij.dot(v_ij);
                n++;
                currentEdge = currentEdge->pair();    // Move to the opposite half-edge
                if (!currentEdge) break;       // If there's no pair, it's a boundary edge
                currentEdge = currentEdge->next();
            }while(currentEdge != startEdge);
            v_rec /= static_cast<double>(n);
            vertex->setPosition(v_rec.point());
        }
    }
    writeMesh("result.obj", mesh);

}

//void Face_Area(const PolyMesh *mesh, std::vector<double>& FaceArea)
//{
//
//}
//void Face_Center(const PolyMesh *mesh, std::vector<MPoint3>& Center)
//{
//
//}


