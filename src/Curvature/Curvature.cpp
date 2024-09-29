#include <optional>
#include "Cuvature.h"

namespace Jade::Curvature {

void CalLocalAveRegion(const std::shared_ptr<PolyMesh const>& mesh,
                       std::vector<double> &vertexLocalAveRegion)
{
    if (!mesh) throw std::invalid_argument("Mesh pointer is null.");

    auto N = mesh->numVertices();
    vertexLocalAveRegion.resize(N, 0.0);

    for (auto fh: mesh->polyfaces()) {
        MHalfedge *he{fh->halfEdge()};
        MHalfedge *he_next{he->next()};
        MHalfedge *he_prev{he->prev()};

        MVert *v_from_he{he->fromVertex()};
        MVert *v_from_he_next{he_next->fromVertex()};
        MVert *v_from_he_prev{he_prev->fromVertex()};

        MVector3 he_tangent{he->tangent()};
        MVector3 he_next_tangent{he_next->tangent()};
        MVector3 he_prev_tangent{he_prev->tangent()};

        std::optional<MVert *> obtuseVertex = std::nullopt;

        // 检查是否存在钝角
        if (vectorAngle(he_tangent, -he_prev_tangent) > PI / 2.0) {
            obtuseVertex = v_from_he;
        } else if (vectorAngle(he_next_tangent, -he_tangent) > PI / 2.0) {
            obtuseVertex = v_from_he_next;
        } else if (vectorAngle(he_prev_tangent, -he_next_tangent) > PI / 2.0) {
            obtuseVertex = v_from_he_prev;
        }

        // 计算面面积并分配给顶点
        MVector3 triangle_edge0{v_from_he_next->position() - v_from_he->position()};
        MVector3 triangle_edge1{v_from_he->position() - v_from_he_prev->position()};
        double triangle_area{0.5 * norm(cross(triangle_edge0, triangle_edge1))};

        if (obtuseVertex) {
            // 如果有钝角，分配面面积给对应的顶点
            for (auto *fv: mesh->polygonVertices(fh)) {
                if (fv == *obtuseVertex) {
                    vertexLocalAveRegion[fv->index()] += triangle_area / 2.0;
                } else {
                    vertexLocalAveRegion[fv->index()] += triangle_area / 4.0;
                }
            }

            auto circusCenter = Jade::CalCircumCenter(v_from_he->position(), v_from_he_prev->position(),
                                                v_from_he_next->position());

            for (const auto &face_he: mesh->polygonHalfedges(fh)) {
                MVector3 edgeMidPoint{0.5 * (face_he->toVertex()->position() - face_he->fromVertex()->position())};
                double edgeLength{face_he->edge()->length()};
                double partArea{0.5 * edgeLength * (edgeMidPoint - circusCenter).norm()};
                vertexLocalAveRegion[face_he->fromVertex()->index()] += 0.5 * partArea;
                vertexLocalAveRegion[face_he->toVertex()->index()] += 0.5 * partArea;
            }
        }
    }
}

void CalMeanCurvature(const std::shared_ptr<PolyMesh const>& mesh,
                   std::vector<double> &vertexLocalAveRegion,
                   std::vector<CurvatureInfo>& curvatures )
{
    std::ofstream f1("./MeanCurvature.txt");
    if (!f1)
        return;

    for (const auto &vertex: mesh->vertices()) {
        const MHalfedge* startEdge = vertex->halfEdge();

        if (startEdge->isBoundary()) continue;
        const MHalfedge* currentEdge = startEdge;
        MPoint3 vi = vertex->position();
        MVector3 lb_vi{0.0, 0.0, 0.0};

        // Traverse the 1-ring neighborhood of the vertex
        do {
            const MHalfedge* nextEdge = currentEdge->pair()->next();
            MPoint3 vj = currentEdge->toVertex()->position();
            MPoint3 vk = nextEdge->toVertex()->position();

            // Compute angles using the cotangent formula
            double angle_vj = Jade::Cotangent(vi - vk, vj - vk);
            double angle_vk = Jade::Cotangent(vk - vj, vi - vj);
            lb_vi += (angle_vk * (vj - vi) +
                      angle_vj * (vk - vi) );
            currentEdge = nextEdge;
        }while(currentEdge!= startEdge);

        lb_vi /= 2.0 * vertexLocalAveRegion[vertex->index()];

        curvatures[vertex->index()].meanCurvature = 0.5 * lb_vi.norm();

        f1 << "Vertex " << vertex->index() << ": Mean Curvature = "
           << curvatures[vertex->index()].meanCurvature << std::endl;
    }
    // Write the Gaussian curvature (gsCurvature) to the file

    std::cout << "Calculate Mean Curvature Done" << std::endl;
}

void CalGaussianCurvature(const std::shared_ptr<PolyMesh const>& mesh,
                          std::vector<double> &vertexLocalAveRegion,
                          std::vector<CurvatureInfo>& curvatures )
{
    std::ofstream f1("./GaussianCurvature.txt");
    if (!f1)
        return;

    for (const auto& vertex : mesh->vertices()){
        int vid = vertex->index();
        curvatures[vid].gsCurvature = 2 * PI / vertexLocalAveRegion[vid];
    }

    for (auto vertex: mesh->vertices()) {
        MHalfedge *startEdge = vertex->halfEdge();

        if (startEdge->isBoundary()) continue;

        MHalfedge* currentEdge = startEdge;
        MPoint3 vi = vertex->position();
        double totalAngle{0.0};

        // Traverse the 1-ring neighborhood of the vertex
        do {
            const MHalfedge* nextEdge = currentEdge->pair()->next();
            MPoint3 vj = startEdge->toVertex()->position();
            MPoint3 vk = nextEdge->toVertex()->position();

            totalAngle += AngleBetween(vj - vi, vk - vi);

            currentEdge = currentEdge->pair()->next();
        } while (currentEdge != startEdge);

        curvatures[vertex->index()].gsCurvature -= totalAngle / vertexLocalAveRegion[vertex->index()];

        f1 << "Vertex " << vertex->index() << ": Gaussian Curvature = "
           << curvatures[vertex->index()].gsCurvature << std::endl;
    }
    std::cout << "Calculate Gaussian Curvature Done" << std::endl;
}
}
