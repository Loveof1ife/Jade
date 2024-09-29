#include "parametrization.h"
namespace Jade::Param {
    std::shared_ptr<PolyMesh> Tutte_Embedding(const std::shared_ptr<const PolyMesh> &mesh) {
        if (!mesh) {
            throw std::invalid_argument("Mesh is null.");
        }
        auto n_vertices = static_cast<int>( mesh->numVertices());
        auto boundary_vertices = mesh->boundaryVertices();

        // const_iterator<MHalfedge* const *, vector<MHalfedge*> >
        auto he_it = mesh->const_halfedge_begin();
        while (!mesh->isBoundary(*he_it))
            ++he_it;

        auto bhe_start = *he_it;
        auto bhe_current = bhe_start;
        int boundary_nums = 0;
        // Traverse the boundary, the boundary of a mesh forms a loop
        do {
            boundary_nums++;
            bhe_current = bhe_current->next();
        } while (bhe_start != bhe_current);

        double delta_angle = 2 * M_PI / boundary_nums;
        double area_sum = Cal_Area(mesh);
        double area_factor = std::sqrt(area_sum / M_PI); // Factor for scaling
        Eigen::VectorXd boundary_of_uv(boundary_nums * 2);  // Allocate space for boundary UV coordinates

        // Set up the boundary UV positions in a circle
        bhe_current = bhe_start;
        for (int i = 0; i < boundary_nums; ++i) {
            auto boundary_index = bhe_current->toVertex()->index();
            boundary_of_uv[boundary_index] = area_factor * std::cos(i * delta_angle);  // x-coordinate (u)
            boundary_of_uv[boundary_index + boundary_nums] =
                    area_factor * std::sin(i * delta_angle);  // y-coordinate (v)
            bhe_current = bhe_current->next();
        }

        SMatrix coff(n_vertices, n_vertices);
        std::vector<Eigen::Triplet<int> > tripletlist;
        Eigen::VectorXd bu = Eigen::VectorXd::Zero(n_vertices);
        Eigen::VectorXd bv = Eigen::VectorXd::Zero(n_vertices);

        for (const auto &v: mesh->vertices()) {
            int v_index = v->index();
            // If the vertex is a boundary vertex, set fixed UV coordinates
            if (mesh->isBoundary(v)) {
                tripletlist.emplace_back(v_index, v_index, 1.0);
                bu[v_index] = boundary_of_uv[v_index];
                bv[v_index] = boundary_of_uv[v_index + boundary_nums];
            } else {
                for (VertexVertexIter v_ring_it = mesh->vv_iter(v); v_ring_it.isValid(); v_ring_it++) {
                    auto v_ring = *v_ring_it;
                    int v_ring_index = v_ring->index();
                    tripletlist.emplace_back(v_index, v_ring_index, 1.0);
                }
                tripletlist.emplace_back(v_index, v_index, static_cast<int>(mesh->valence(v)));
            }
        }
        coff.setFromTriplets(tripletlist.begin(), tripletlist.end());
        Eigen::SparseLU<SMatrix> solver;
        solver.compute(coff);
        Eigen::VectorXd u = solver.solve(bu);
        Eigen::VectorXd v = solver.solve(bv);
        Eigen::MatrixXd uv(n_vertices, 2);

        uv.col(0) = u;  // U coordinates
        uv.col(1) = v;  // V coordinates

        auto para_mesh = std::make_shared<PolyMesh>();
        for (int i = 0; i < n_vertices; ++i) {
            para_mesh->addVertex(uv(i, 0), uv(i, 1), 0.0);  // Add vertex in 2D space (u, v, 0)
        }
        return para_mesh;
    }

//    void Cal_Local_Cord(const PolyMesh *pmesh, Eigen::MatrixXd &local_cord) {
//        size_t face_num = pmesh->numPolygons();
//        local_cord.resize(face_num, 3 * 2);
//
//        for (size_t i = 0; i < face_num; i++) {
//            //MPolyFace* face_current = *face_it;
//            MPolyFace *face_current = pmesh->polyface(i);
//            MVector3 normal = face_current->normal();
//            int face_index = face_current->index();
//            auto fv_it = pmesh->fv_iter(face_current);
//
//            auto p0 = (*fv_it)->position();
//            fv_it++;
//            auto p1 = (*fv_it)->position();
//            fv_it++;
//            auto p2 = (*fv_it)->position();
//
//            MVector3 e0 = (p1 - p0);
//            MVector3 e1 = (p2 - p0);
//
//            auto x_ = e0.normalized();
//            auto y_ = cross(x_, normal);
//
//            local_cord.row(face_index) << 0, 0,
//                    e0.norm(), 0,
//                    x_.dot(e1), y_.dot(e1);
//        }
//
//    }
//
//
//    void Cal_Jacobi_Uv(PolyMesh *pmesh, Eigen::MatrixX2d &uv, Eigen::MatrixXd &local_cord, Eigen::Matrix2d &P,
//                  Eigen::Matrix2d &S, Eigen::Matrix2d &J, std::vector<Eigen::Matrix2d> &l_t) {
//        size_t face_num = pmesh->numPolygons();
//
//        for (size_t i = 0; i < face_num; i++) {
//            MPolyFace *face_current = pmesh->polyface(i);
////      MPolyFace *face_current = *face_it;
//            auto fv_it = pmesh->fv_iter(face_current);
//            int face_index = face_current->index();
//
//            auto i0 = (*fv_it)->index();
//            ++fv_it;
//            auto i1 = (*fv_it)->index();
//            ++fv_it;
//            auto i2 = (*fv_it)->index();
//
//
//            Eigen::Matrix2d P, S, J;
//            S << local_cord(face_index, 2) - local_cord(face_index, 0), local_cord(face_index, 4) -
//                                                                        local_cord(face_index, 0),
//                    local_cord(face_index, 3) - local_cord(face_index, 1), local_cord(face_index, 5) -
//                                                                           local_cord(face_index, 1);
//
//
//            P << uv(i1, 0) - uv(i0, 0), uv(i2, 0) - uv(i0, 0),
//                    uv(i1, 1) - uv(i0, 1), uv(i2, 1) - uv(i0, 1);
//
//            J = P * S.inverse();
//
//            Eigen::JacobiSVD<Eigen::MatrixXd> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
//            Eigen::MatrixXd U = svd.matrixU();
//            Eigen::MatrixXd V = svd.matrixV();
//            Eigen::Matrix2d R = U * V.transpose();
//            if (R.determinant() < 0) {
//                U(0, 1) = -U(0, 1);
//                U(1, 1) = -U(1, 1);
//                R = U * V.transpose();
//            }
//            //l_t best rotation: u v_transpose
//            l_t[face_index] = R;
//            // face_it++;
//        }
//    }
//
//    void Tutte_Embedding(PolyMesh *mesh) {
//        int n_vertices = mesh->numVertices();
//        auto boundary_vertices = mesh->boundaryVertices();
//
//        int boundary_nums{0};
//        auto he_it = mesh->halfedge_begin();
//        //Finding the First Boundary Half-Edge:
//        while (!mesh->isBoundary(*he_it))
//            he_it++;
//        //Setting the Start of the Boundary:
//        auto bhe_start = *he_it;
//        auto bhe_current = bhe_start;
//
//        //traverse the boundary. the boundary of a mesh forms a loop
//        do {
//            boundary_nums++;
//            bhe_current = bhe_start->next();
//
//        } while (bhe_start != bhe_current);
//
//        double delta_angle = 2 * M_PI / boundary_nums;
//        double area_sum = Cal_Area(mesh);
//        double area_1_factor = sqrt(area_sum / M_PI);
//        Eigen::VectorXd boundary_of_uv;
//        boundary_of_uv.resize(boundary_nums * 2);
//
//        bhe_current = bhe_start;
//        for (int i = 0; i < boundary_nums; ++i) {
//            auto boundary_index = bhe_current->toVertex()->index();
//            boundary_of_uv[boundary_index] = cos(i * delta_angle);
//            boundary_of_uv[boundary_index + boundary_nums] = sin(i * delta_angle);
//            bhe_current = bhe_current->next();
//        }
//
//        SMatrix coff(n_vertices, n_vertices);
//        std::vector<Eigen::Triplet<double> > tripletlist;
//        Eigen::VectorXd bu = Eigen::VectorXd::Zero(n_vertices);
//        Eigen::VectorXd bv = Eigen::VectorXd::Zero(n_vertices);
//
//        for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
//            int v_it_index = (*v_it)->index();
//            //boundary
//            if (mesh->isBoundary(*v_it)) {
//                Eigen::Triplet<double> triplet(v_it_index, v_it_index, 1);
//                tripletlist.push_back(triplet);
//                bu[v_it_index] = boundary_of_uv[v_it_index];
//                bu[v_it_index] = boundary_of_uv[v_it_index + n_vertices];
//            }
//            //1-ring
//            for (auto v_ring = mesh->vv_iter(*v_it); v_ring.isValid(); v_ring++) {
//                Eigen::Triplet<double> triplet(v_it_index, (*v_ring)->index(), 1);
//                tripletlist.push_back(triplet);
//            }
//            //vi itself
//            Eigen::Triplet<double> triplet(v_it_index, v_it_index, static_cast<double>(mesh->valence((*v_it))));
//            tripletlist.push_back(triplet);
//        }
//        //set coff matrix
//        coff.setFromTriplets(tripletlist.begin(), tripletlist.end());
//        //solve
//        Eigen::SparseLU<SMatrix> solver;
//        solver.compute(coff);
//        Eigen::VectorXd u = solver.solve(bu);
//        Eigen::VectorXd v = solver.solve(bv);
//
//        Eigen::MatrixX2d uv;
//        uv.resize(n_vertices, 2);
//
//        uv.col(0) = u;
//        uv.col(1) = v;
//
//        auto para_mesh = new PolyMesh();
//        for (int i = 0; i < n_vertices; i++) {
//            para_mesh->addVertex(uv(i, 0), uv(i, 1), 0);
//        }
//
//        writeMesh("bunny_tutte_para.obj", para_mesh);
//    }
//
//    void Tutte_Embedding(Eigen::MatrixX2d &uv, PolyMesh &mesh) {
//        int F_N = mesh.numPolygons();
//        int V_N = mesh.numVertices();
//        //calc surface area
//        double area_sum = 0;
//        for (int i = 0; i < F_N; ++i) {
//            auto f_h = mesh.polyface(i);
//            auto itfv = mesh.fv_iter(f_h);
//            auto v0 = (*itfv)->position();
//            ++itfv;
//            auto v1 = (*itfv)->position();
//            ++itfv;
//            auto v2 = (*itfv)->position();
//
//            auto e0 = v1 - v0;
//            auto e1 = v2 - v0;
//            auto avec = cross(e0, e1);
//            area_sum += avec.norm() / 2.0;
//        }
//
//        //set the boundary vertices to circle
//        int boundary_num = 0;
//        auto it1 = mesh.halfedge_begin();
//        while (!mesh.isBoundary(*it1))
//            it1++;
//        auto he_start = *it1;
//        auto he_it = he_start;
//        do {
//            he_it = (he_it)->next();
//            boundary_num++;
//        } while (he_it != he_start);
//
//        double delta_angle = 2 * M_PI / boundary_num;
//        double area_1_factor = sqrt(area_sum / M_PI);
//        Eigen::VectorXd position_of_mesh;
//        position_of_mesh.resize(2 * V_N);
//        for (int i = 0; i < boundary_num; ++i) {
//            auto v_h = he_start->toVertex();
//            position_of_mesh(v_h->index()) = area_1_factor * cos(i * delta_angle);
//            position_of_mesh(v_h->index() + V_N) = area_1_factor * sin(-i * delta_angle);
//            he_start = he_start->next();
//        }
//
//
//        //calc the matrix
//        typedef Eigen::Triplet<double> T;
//        typedef Eigen::SparseMatrix<double> SMatrix;
//
//        std::vector<T> tripletlist;
//        Eigen::VectorXd bu = Eigen::VectorXd::Zero(V_N);
//        Eigen::VectorXd bv = Eigen::VectorXd::Zero(V_N);
//        for (auto it1 = mesh.vertices_begin(); it1 != mesh.vertices_end(); it1++) {
//            int it1idx = (*it1)->index();
//            if (mesh.isBoundary(*it1)) {
//                tripletlist.push_back(T(it1idx, it1idx, 1));
//                auto point = (*it1)->position();
//                bu(it1idx) = position_of_mesh[it1idx];
//                bv(it1idx) = position_of_mesh[it1idx + V_N];
//            } else {
//                for (auto it2 = mesh.vv_iter(*it1); it2.isValid(); ++it2) {
//                    tripletlist.push_back(T(it1idx, (*it2)->index(), -1));
//                }
//                tripletlist.push_back(T(it1idx, it1idx, mesh.valence(*it1)));
//            }
//        }
//
//        SMatrix coff(V_N, V_N);
//        coff.setFromTriplets(tripletlist.begin(), tripletlist.end());
//        Eigen::SparseLU<SMatrix> solver;
//        solver.compute(coff);
//        Eigen::VectorXd xu = solver.solve(bu);
//        Eigen::VectorXd xv = solver.solve(bv);
//
//        //solve the inner positions
//        uv.col(0) = xu;
//        uv.col(1) = xv;
//    }
//
//    void Local_Global_Approach(PolyMesh *pmesh) {
//        //init: Tutte Embedding uv
//        //tutte initialization;
//        Eigen::MatrixX2d uv;
//        Tutte_Embedding(uv, *pmesh);
//
//
//        //Optimize l_t
//        //once get uv -> jacobi -> l_t = U V_T
//        //for each face(x,y,z)
//        //cal local cord(x_local, y_local),
//        auto face_it = pmesh->polyfaces_begin();
//        size_t face_num = pmesh->numPolygons();
//
//        Eigen::MatrixXd local_cord;
//        Cal_Local_Cord(pmesh, local_cord);
//
//        //(u,v) = jac * (x_local, y_local)
//        std::vector<Eigen::Matrix2d> l_t(face_num);
//        Eigen::Matrix2d P, S, J;
//        l_t = Cal_Jacobi_Uv(pmesh, uv, local_cord, P, S, J, l_t);
//
//        //optim j_t
//        //j_t(uv):(x_local, y_local) -> (uv)
//        // | j_t - l_t |(uv) linear equation
//        std::vector<Eigen::Triplet<double> > triplets;
//        Eigen::SparseMatrix<double> sparse_mat;
//        Eigen::VectorXd bu, bv;
//
//        size_t num_vertex = pmesh->numVertices();
//        sparse_mat.resize(num_vertex, num_vertex);
//        bu.setZero(num_vertex);
//        bv.setZero(num_vertex);
//
//        //calculate cots
//        std::vector<double> cots(pmesh->numHalfEdges());
//        for (auto he_it = pmesh->halfedge_begin(); he_it != pmesh->halfedge_end(); he_it++) {
//            if (pmesh->isBoundary(*he_it))
//                continue;
//
//            auto v1 = (*he_it)->fromVertex()->position();
//            auto v2 = (*he_it)->toVertex()->position();
//            auto v0 = (*he_it)->next()->toVertex()->position();
//
//            auto e0 = v1 - v0;
//            auto e1 = v2 - v0;
//            double cotangle = dot(e0, e1) / cross(e0, e1).norm();
//
//            cots[(*he_it)->index()] = cotangle;
////        int vid0 = (*he_it)->fromVertex()->index();
////        int vid1 = (*he_it)->toVertex()->index();
////        trivec.emplace_back(vid0, vid0, cotangle);
////        trivec.emplace_back(vid1, vid1, cotangle);
////        trivec.emplace_back(vid0, vid1, -cotangle);
////        trivec.emplace_back(vid1, vid0, -cotangle);
//        }
//
//        // face_it = pmesh->polyfaces_begin();
//        for (size_t i = 0; i < face_num; i++) {
//            MPolyFace *face_current = pmesh->polyface(i);
//            auto face_vertex_it = pmesh->fv_iter(face_current);
//            auto face_he = face_current->halfEdge();
//            int face_index = face_current->index();
//
//            for (size_t j = 0; j < 3; j++) {
//                auto v0 = (*face_vertex_it);
//                face_vertex_it++;
//                auto v1 = (*face_vertex_it);
//                //cot_ij+cot_ji
//                double cot_ij = cots[face_he->index()];
//                face_he = face_he->next();
//                int v0_id = v0->index();
//                int v1_id = v1->index();
//                triplets.emplace_back(v0_id, v0_id, cot_ij);
//                triplets.emplace_back(v0_id, v1_id, -cot_ij);
//                triplets.emplace_back(v1_id, v0_id, -cot_ij);
//                triplets.emplace_back(v1_id, v1_id, cot_ij);
//
//                Eigen::Vector2d local_edge;
//                switch (j) {
//                    case 0:
//                        local_edge << local_cord(face_index, 2), local_cord(face_index, 3);
//                        break;
//                    case 1:
//                        local_edge << local_cord(face_index, 4) - local_cord(face_index, 2), local_cord(face_index, 5) -
//                                                                                             local_cord(face_index, 3);
//                        break;
//                    case 2:
//                        local_edge << -local_cord(face_index, 4), -local_cord(face_index, 5);
//
//                }
//                Eigen::Vector2d b = cot_ij * l_t[face_index] * (local_edge);
//                bu[v0_id] -= b[0];
//                bv[v0_id] -= b[1];
//                bu[v1_id] += b[0];
//                bv[v1_id] += b[1];
//            }
//
//        }
//        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//        sparse_mat.setFromTriplets(triplets.begin(), triplets.end());
//        solver.compute(sparse_mat);
//        uv.col(0) = solver.solve(bu);
//        uv.col(1) = solver.solve(bv);

//    for (int i = 0; i < nv; i++)
//    {
//        auto v_h = mesh.vert(i);
//        v_h->setPosition(uv(i, 0), uv(i, 1), 0);
//    }
//
//    writeMesh(out_path, &mesh);
//}
    double Cal_Area(const std::shared_ptr<const PolyMesh>& mesh){
        auto f_it = mesh->const_polyfaces_begin();
        auto f_start = (*f_it);
        auto f_current = f_start;
        double area{0.0};

        do {
            auto e1 = f_current->halfEdge();
            auto e2 = e1->next();
            auto p1 = e1->fromVertex()->position();
            auto p2 = e1->toVertex()->position();
            auto p3 = e2->toVertex()->position();
            MVector3 v21(p2 - p1), v31(p3 - p1);
            area += v21.cross(v31).norm();
            f_it++;
        } while ((f_current = *f_it) != f_start);
        return 0.5 * area;
    }
}

