#include "Cuvature.h"

MVector3 cal_circum_enter(const MVector3& a, const MVector3& b, const MVector3& c)
{
    MVector3 ac = c - a, ab = b - a;
    MVector3 abXac = cross(ab, ac), abXacXab = cross(abXac, ab), acXabXac = cross(ac, abXac);
    return a + (abXacXab * ac.normSq() + acXabXac * ab.normSq()) / (2.0 * abXac.normSq());
}

double cotangent(const MVector3& v1, const MVector3& v2) {
    double cosTheta = v1.dot(v2) / (v1.norm() * v2.norm());
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    return cosTheta / sinTheta;
}

double angleBetween(const MVector3& v1, const MVector3& v2) {
    double dotProduct = v1.dot(v2);
    double normsProduct = std::sqrt(v1.norm() * v2.norm());
    double cosTheta = dotProduct / normsProduct;
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta)); // 防止超出 acos 的输入范围
    return std::acos(cosTheta); // 返回弧度值
}


void cal_local_ave_region(const PolyMesh *mesh, std::vector<double> &vertexLAR)
{
    auto N = mesh->numVertices();
    vertexLAR.resize(N, 0.0);
    for(auto fh : mesh->polyfaces()){
        MHalfedge *he = fh->halfEdge();
        MHalfedge *he_next = he->next();
        MHalfedge *he_prev = he->prev();

        MVert *v_from_he = he->fromVertex();
        MVert *v_from_he_next = he_next->fromVertex();
        MVert *v_from_he_prev = he_prev->fromVertex();

        MVector3 vec_he_nor = he->tangent();
        MVector3 vec_he_next_nor = he_next->tangent();
        MVector3 vec_he_prev_nor = he_prev->tangent();

        bool isObtuseAngle = false;
        MVert *obtuseVertexHandle = nullptr;

        if (vectorAngle(vec_he_nor, -vec_he_prev_nor) > PI / 2.0)
        {
            isObtuseAngle = true;
            obtuseVertexHandle = v_from_he;
        }
        else if (vectorAngle(vec_he_next_nor, -vec_he_nor) > PI / 2.0)
        {
            isObtuseAngle = true;
            obtuseVertexHandle = v_from_he_next;
        }
        else if (vectorAngle(vec_he_prev_nor, -vec_he_next_nor) > PI / 2.0)
        {
            isObtuseAngle = true;
            obtuseVertexHandle = v_from_he_prev;
        }
        if(obtuseVertexHandle) {
            for (auto fhh: mesh->polygonHalfedges(fh)) {
                MVector3 fromVertexPos = fhh->fromVertex()->position();

                MVector3 edge1Midpoint = 0.5 * (fromVertexPos + fhh->toVertex()->position());
                MVector3 edge2Midpoint = 0.5 * (fromVertexPos + fhh->prev()->toVertex()->position());
                MVector3 oppsMidpoint = 0.5 * (fhh->next()->toVertex()->position() + fhh->next()->fromVertex()->position());

                MVector3 eMv1 = (edge1Midpoint - fromVertexPos), eMv2 = (edge2Midpoint - fromVertexPos);
                MVector3 oMv = (oppsMidpoint - fromVertexPos);

                vertexLAR[fhh->fromVertex()->index()] = cross(eMv1, oMv).norm() + cross(eMv2, oMv).norm();
            }
        }
        else
        {
            auto cCenter = cal_circum_enter(v_from_he->position(), v_from_he_next->position(), v_from_he_prev->position());

            for(auto fhh : mesh->polygonHalfedges(fh)){
                double edgeLength = fhh->edge()->length();
                auto edgeMidpoint =0.5 * (fhh->fromVertex()->position() + fhh->toVertex()->position());
                double partArea = 0.5 * edgeLength * (edgeMidpoint - cCenter).norm();
                vertexLAR[fhh->fromVertex()->index()] += partArea;
                vertexLAR[fhh->toVertex()->index()] += partArea;
            }
        }
    }
}
void cal_mean_curvature(const PolyMesh *mesh, std::vector<double> &vertexLAR)
{
    std::ofstream	f1("./MeanCurvature.txt");
    if (!f1)
        return;
    std::ofstream	f2("./AbsoluteMeanCurvature.txt");
    if (!f2)
        return;
    // for every vertex
    for(auto vertex : mesh->vertices()) {
        MHalfedge* startEdge = vertex->halfEdge();

        if (startEdge->isBoundary())
            continue;

        MHalfedge* currentEdge = startEdge;
        MPoint3 vi = vertex->position();
        MVector3 lb_vi{0.0,0.0,0.0};
        do
        {
//            MPoint3 vj = currentEdge->toVertex()->position();
//
//            MHalfedge* nextEdge = currentEdge->pair()->next();
//            MHalfedge* preEdge = currentEdge->pair()->prev();
//
//            MPoint3 vk1 = nextEdge->toVertex()->position();
//            MPoint3 vk2 = preEdge->fromVertex()->position();
//
//            double angle1_vj = cotangent(vj - vk1, vi - vk1);
//            double angle2_vj = cotangent(vj - vk2, vi - vk2);
//            lb_vi += (angle1_vj + angle2_vj) * (vj - vi);
//
//            currentEdge = nextEdge;
            MPoint3 vj = currentEdge->toVertex()->position();
            MHalfedge* nextEdge = currentEdge->pair()->next();
            MPoint3 vk = nextEdge->toVertex()->position();

            double angle_vj = cotangent(vi - vk, vj - vk);
            double angle_vk = cotangent(vk - vj, vi - vj);

            lb_vi += (angle_vj) * (vj - vi);
            lb_vi += (angle_vk) * (vk - vi);
            currentEdge = currentEdge->pair()->next();
        }while(currentEdge != startEdge);

        lb_vi /= 2 * vertexLAR[vertex->index()];
        if(lb_vi.dot(vertex->normal()))
            f1<< (-0.5) * lb_vi.norm()<<std::endl;
        else
            f1<< (0.5) * lb_vi.norm()<<std::endl;

        f2 << (0.5) * lb_vi.norm() << std::endl;
    }
    std::cout << "Calculate Mean Curvature Done" << std::endl;
    std::cout << "Calculate Absolute Mean Curvature Done" << std::endl;
}
void cal_gaussian_curvature(const PolyMesh *mesh, std::vector<double> &vertexLAR)
{
    std::ofstream f1("./GaussianCurvature.txt");
    if (!f1)
        return;

    std::vector<double> gaussian_curvature(mesh->numVertices(), 0.0);

    for(auto vertex : mesh->vertices()){
        MHalfedge* startEdge = vertex->halfEdge();

        if (startEdge->isBoundary())
            continue;

        MHalfedge* currentEdge = startEdge;
        MPoint3 vi = vertex->position();
        double totalAngle{0.0};

        do{
            MPoint3 vj = startEdge->toVertex()->position();
            MPoint3 vk = startEdge->next()->toVertex()->position();

            totalAngle += angleBetween(vj - vi, vk - vi);

            currentEdge = currentEdge->pair()->next();
        }while(currentEdge != startEdge);

        gaussian_curvature[vertex->index()] =  (2 * M_PI - totalAngle) / vertexLAR[vertex->index()];

        f1 << gaussian_curvature[vertex->index()]  << std::endl;  // Write Gaussian curvature to file
    }
    std::cout << "Calculate Gaussian Curvature Done" << std::endl;

}
