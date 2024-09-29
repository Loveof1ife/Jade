#include "Filter.h"

namespace Jade::Filter{

void fitVertexToNormal(std::vector<MPoint3>& vertexPositions,
                       const std::shared_ptr<PolyMesh>& in_mesh,
                       const std::vector<MVector3>& filterNormalField,
                       const std::vector<MPoint3>& faceCenter,
                       const std::vector<double>& faceArea,
                       int maxIterations , double convergenceThreshold )

{
    if ( filterNormalField.empty()) throw std::invalid_argument("Invalid input mesh or normal field.");

    size_t vN = in_mesh->numVertices();

    // Initialize vertex positions (copy from input mesh)
    vertexPositions.resize(vN);
    for (auto& v : in_mesh->vertices()) {
        vertexPositions[v->index()] = v->position();
    }

    // Iterate with Gauss-Seidel updates
    for (int iter = 0; iter < maxIterations; ++iter){
        double maxDisplacement = 0.0;  // Track maximum change in vertex positions for convergence check

        for (auto& v : in_mesh->vertices()){
            int vid { v->index()}, neiCount{0} ;
            MPoint3 newVertex{v->position() }, deltaV{0, 0, 0};

            for(auto vf_iter {in_mesh->vf_iter(v)}; vf_iter.isValid(); vf_iter++, neiCount++){
                int nei_fid = (*vf_iter)->index();
                const MVector3& neiFaceNormal {filterNormalField[nei_fid] };
                const MPoint3& neiFaceCenter {faceCenter[nei_fid] };
                deltaV = deltaV + neiFaceNormal * dot(neiFaceNormal, (neiFaceCenter - v->position() ) );
            }

            newVertex += deltaV / neiCount;
            v->setPosition(newVertex);
        }
        // Check for convergence
        if (maxDisplacement < convergenceThreshold) {
            break;  // Convergence reached
        }
    }

}

void BilateralNormalFiltering(const std::shared_ptr<PolyMesh>& in_mesh,
                              std::vector<MVector3>& filterNormalField,
                              std::vector<MPoint3>& faceCenter,
                              std::vector<double>& faceArea
                              )
{
        if (!in_mesh ) throw std::invalid_argument("Empty Mesh pointers ");

        size_t fN{in_mesh->numPolygons()};
        faceCenter.resize(fN), faceArea.resize(fN);
        std::vector<MVector3> initNormalField;

        initNormalField.resize(fN), filterNormalField.resize(fN);

        for (auto face : in_mesh->polyfaces() ) {
            int fid = face->index();
            MVector3 faceNormal = face->normal();

            initNormalField[fid] = faceNormal;

            std::vector<MVert* > vList;
            for(FaceVertexIter fv_iter = in_mesh->fv_iter(face); fv_iter.isValid(); fv_iter++ ){
                vList.emplace_back(*fv_iter);
            }
            MVector3 e12 = vList[1]->position() - vList[0]->position();
            MVector3 e13 = vList[2]->position() - vList[0]->position();
            faceCenter[fid] = in_mesh->calculatFaceCenter(face);
            faceArea[fid] = 0.5 * cross(e12, e13).norm() ;
        }

        double sigmaCenter{0.0}, sigmaNormal{0.0};

        for (auto face : in_mesh->polyfaces() ) {
            int fid{face->index()};
            MVector3 normal{face->normal()};

            for (auto faceIter = in_mesh->ff_iter(face); faceIter.isValid(); faceIter++) {
                int nei_fid = (*faceIter)->index();
                sigmaCenter += (faceCenter[fid] - faceCenter[nei_fid]).norm();
                sigmaNormal += ((*faceIter)->normal() - normal).norm();
            }

        }
        double numFaces {static_cast<double > (in_mesh->numPolygons() ) };
        sigmaCenter /= (numFaces * 3.0);
        sigmaNormal /= (numFaces * 3.0);

        for (auto face : in_mesh->polyfaces()) {
            int fid {face->index() };
            MVector3 filterNormal(0, 0, 0);
            double Kp {0.0};

            filterNormal += initNormalField[fid] * faceArea[fid];
            Kp += faceArea[fid];

            for(auto neighbor_face = in_mesh->ff_iter(face); neighbor_face.isValid(); neighbor_face++) {
                int nei_fid {(*neighbor_face)->index() };
                auto deltaCenter {faceCenter[fid] - faceCenter[nei_fid]};
                auto deltaNormal {initNormalField[fid] - initNormalField[nei_fid]};
                double Aj { faceArea[nei_fid]};

                double Ws {Jade::Gaussian(deltaCenter, sigmaCenter)};
                double Wr {Jade::Gaussian(deltaNormal, sigmaNormal)};
                filterNormal += Aj * Ws * Wr * initNormalField[nei_fid];
                Kp += Aj * Ws * Wr;
            }

            filterNormalField[fid] = filterNormal/ Kp;
            filterNormalField[fid].normalize();
        }
}
}

//void BilateralNormalFiltering(const std::shared_ptr<PolyMesh const>& in_mesh,
//                              const std::shared_ptr<PolyMesh>& out_mesh,
//                              std::vector<MVector3>& filterNormalField)
//{
//    if (!in_mesh || !out_mesh) {
//        throw std::invalid_argument("Mesh pointers cannot be null.");
//    }
//
//    //step1: 计算面法向,法向滤波：面i的新的法向等于i与其周围1-邻域面的法向的加权（GS）平均
//    for(auto face : in_mesh->polyfaces()){
//        MVector3 ni = face->normal();
//        MPoint3 ci = face->getFaceCenter();
//        auto startEdge = face->halfEdge();
//        auto currentEdge = startEdge;
//        MVector3 totalNormal{0.0, 0.0, 0.0};
//
//        double K{0.0};
//        do{
//            auto face_j = currentEdge->pair()->polygon();
//            double Kj{0.0};
//
//            if(face_j == nullptr) // boundary
//                continue;
//            MVector3 nj = face_j->normal();
//            MPoint3 cj = face_j->getFaceCenter();
//
//            MHalfedge* hedge1 = face_j->halfEdge();
//            MHalfedge* hedge2 = hedge1->next();
//            auto e1 = hedge1->toVertex()->position() - hedge1->fromVertex()->position();
//            auto e2 = hedge2->toVertex()->position() - hedge2->fromVertex()->position();
//            double Aj = e1.cross(e2).norm();
//
//            Kj = Aj * GaussianFunc(ci, cj) * GaussianFunc(ni, nj);
//            totalNormal += Kj * nj;
//            K += Kj;
//            currentEdge = currentEdge->next();
//        }while(currentEdge != startEdge);
//        totalNormal /= K;
//        face->setNormal(totalNormal);
//    }
//    //根据滤波后的面法向重建PolyMesh
//
//    for(size_t i = 0; i < IteratorTimes; i++)
//    {
//
//        for(auto vertex : mesh->vertices())
//        {
//
//            MHalfedge* startEdge = vertex->halfEdge();  // Get the starting half-edge
//            MHalfedge* currentEdge = startEdge;
//            MVector3 v_rec{0.0, 0.0, 0.0};
//            size_t n{0};
//            do{
//
//                auto vertex_j = currentEdge->toVertex();
//                MVector3 v_ij(vertex_j->position() - vertex->position());
//                MVector3 n_ij = currentEdge->polygon()->normal();
//                v_rec += n_ij * n_ij.dot(v_ij);
//                n++;
//                currentEdge = currentEdge->pair();    // Move to the opposite half-edge
//                if (!currentEdge) break;       // If there's no pair, it's a boundary edge
//                currentEdge = currentEdge->next();
//            }while(currentEdge != startEdge);
//            v_rec /= static_cast<double>(n);
//            vertex->setPosition(v_rec.point());
//        }
//    }
//    writeMesh("result.obj", mesh);
//
//}

//void Face_Area(const PolyMesh *mesh, std::vector<double>& FaceArea)
//{
//
//}
//void Face_Center(const PolyMesh *mesh, std::vector<MPoint3>& Center)
//{
//
//}


