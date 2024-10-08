#include "parametrization.h"
namespace Jade::Param {

    Parameterizer::Parameterizer(const std::string& meshPath){
        mesh = std::make_shared<PolyMesh>();
        loadMesh(meshPath, mesh.get());
        if (!mesh) throw std::invalid_argument("Mesh is null.");
        n_vertices = getNumofVertices();
    }

    double Parameterizer::calPolymeshArea() const  {
        double area_sum{0.0};
        int F_N = mesh->numPolygons();
        for (int i = 0; i < F_N; ++i)
        {
            auto f_h = mesh->polyface(i);
            auto itfv = mesh->fv_iter(f_h);
            auto v0 = (*itfv)->position();
            ++itfv;
            auto v1 = (*itfv)->position();
            ++itfv;
            auto v2 = (*itfv)->position();

            auto e0 = v1 - v0;
            auto e1 = v2 - v0;
            auto avec = cross(e0, e1);
            area_sum += avec.norm() / 2.0;
        }
        return area_sum;
    }

    Eigen::VectorXd  Parameterizer::setBoundaryVerticesToCircle(){
        auto he_it = mesh->halfedge_begin();
        while (!mesh->isBoundary(*he_it) && he_it != mesh->halfedge_end())
            he_it++;

        auto bhe_start = *he_it;
        auto bhe_current = bhe_start;
        int boundary_nums = 0;
        // Traverse the boundary, the boundary of a mesh forms a loop
        do {
            boundary_nums++;
            bhe_current = bhe_current->next();
        } while (bhe_start != bhe_current);

        //if(boundary_nums != boundary_vertices.size()) std::cerr << "Error: Boundary Cal Error"  << std::endl;

        double delta_angle = 2 * M_PI / boundary_nums;
        double area_factor = std::sqrt(area / M_PI); // Factor for scaling
        Eigen::VectorXd boundary_of_uv(boundary_nums * 2);  // Allocate space for boundary UV coordinates

        // Set up the boundary UV positions in a circle
        bhe_current = bhe_start;
        for (int i = 0; i < boundary_nums; ++i) {
            auto boundary_index = bhe_current->toVertex()->index();
            boundary_of_uv[boundary_index] = area_factor * std::cos(i * delta_angle);  // x-coordinate (u)
            boundary_of_uv[boundary_index + n_vertices] =
                    area_factor * std::sin(i * delta_angle);  // y-coordinate (v)
            bhe_current = bhe_current->next();
        }
        return boundary_of_uv;
    }

    void Parameterizer::solveUV(Eigen::VectorXd boundary_of_uv) {
        std::vector<Eigen::Triplet<int>> tripletlist;
        Eigen::VectorXd bu = Eigen::VectorXd::Zero(n_vertices);
        Eigen::VectorXd bv = Eigen::VectorXd::Zero(n_vertices);
        for (const auto &v: mesh->vertices()) {
            int v_index = v->index();
            // If the vertex is a boundary vertex, set fixed UV coordinates
            if (mesh->isBoundary(v)) {
                tripletlist.emplace_back(v_index, v_index, 1.0);
                bu[v_index] = boundary_of_uv[v_index];
                bv[v_index] = boundary_of_uv[v_index + boundary_vertices.size()];
            } else {
                for (VertexVertexIter v_ring_it = mesh->vv_iter(v); v_ring_it.isValid(); v_ring_it++) {
                    auto v_ring = *v_ring_it;
                    int v_ring_index = v_ring->index();
                    tripletlist.emplace_back(v_index, v_ring_index, 1.0);
                }
                tripletlist.emplace_back(v_index, v_index, static_cast<int>(mesh->valence(v)));
            }
        }
        SparseMatrix coff(n_vertices, n_vertices);
        coff.setFromTriplets(tripletlist.begin(), tripletlist.end());
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute(coff);
        Eigen::VectorXd xu = solver.solve(bu);
        Eigen::VectorXd xv = solver.solve(bv);
        tutteParameterization.resize(n_vertices, 2);
        tutteParameterization.col(0) = xu;
        tutteParameterization.col(1) = xv;
    }
    std::shared_ptr<PolyMesh> Parameterizer::applyUVtoMesh() const{
        auto uvMesh{std::make_shared<PolyMesh>() };
        for (int i = 0; i < n_vertices; i++)
        {
            auto v_h = uvMesh->vert(i);
            v_h->setPosition(iterParameterization(i, 0), iterParameterization(i, 1), 0);
        }
        return uvMesh;
    }

    void Parameterizer::computeLocalCoordinates(){
        n_faces = mesh->numPolygons();
        localCoord.resize(n_faces, 3*2);
        cotangent_weights.resize(mesh->numHalfEdges());

        #pragma omp parallel for
        for (int i = 0; i < n_faces; ++i){
            auto face = mesh->polyface(i);
            auto normal = face->normal();
            auto fv_iter = mesh->fv_iter(face);
            auto v0 = (*fv_iter)->position();
            ++fv_iter;
            auto v1 = (*fv_iter)->position();
            ++fv_iter;
            auto v2 = (*fv_iter)->position();

            MVector3 e1 = v1 - v0;
            MVector3 e2 = v2 - v0;
            MVector3 x_ = e1.normalized();
            MVector3 y_ = cross(x_, normal);
            localCoord.row(i) <<    0, 0,
                                    e1.norm(), 0,
                                    dot(e2, x_), dot(e2, y_);
        }
    }
    void Parameterizer::computeCotangentWeights() {
        cotangent_weights.resize(mesh->numHalfEdges());
        std::vector<Eigen::Triplet<double> > trivec;

        for (auto he_iter = mesh->halfedge_begin(); he_iter != mesh->halfedge_end(); ++he_iter) {
            if (mesh->isBoundary(*he_iter)) continue;

            auto v1 = (*he_iter)->fromVertex()->position();
            auto v2 = (*he_iter)->toVertex()->position();
            auto v0 = (*he_iter)->next()->toVertex()->position();

            auto e0 = v1 - v0;
            auto e1 = v2 - v0;
            double cotangle = dot(e0, e1) / cross(e0, e1).norm();

            cotangent_weights[(*he_iter)->index()] = cotangle;
            int vid0 = (*he_iter)->fromVertex()->index();
            int vid1 = (*he_iter)->toVertex()->index();

            trivec.emplace_back(vid0, vid0, cotangle);
            trivec.emplace_back(vid1, vid1, cotangle);
            trivec.emplace_back(vid0, vid1, -cotangle);
            trivec.emplace_back(vid1, vid0, -cotangle);
        }

        laplaceMatrix.resize(n_vertices, n_vertices);
        laplaceMatrix.setFromTriplets(trivec.begin(), trivec.end());

        solver.compute(laplaceMatrix);
    }

    void Parameterizer::solveLocalStep() {
        #pragma omp parallel for
        for (int i = 0; i < n_faces; ++i) {
            auto f_h = mesh->polyface(i);
            auto fv_iter = mesh->fv_iter(f_h);
            auto i0 = (*fv_iter)->index();
            ++fv_iter;
            auto i1 = (*fv_iter)->index();
            ++fv_iter;
            auto i2 = (*fv_iter)->index();

            Eigen::Matrix2d P, S, J;

            P <<    iterParameterization(i1, 0) - iterParameterization(i0, 0), iterParameterization(i2, 0) - iterParameterization(i0, 0),
                    iterParameterization(i1, 1) - iterParameterization(i0, 1), iterParameterization(i2, 1) - iterParameterization(i0, 1);
            // i face, i1(0,1) i2(2,3) i3(4,5)
            S <<    localCoord(i, 2) - localCoord(i, 0), localCoord(i, 4) - localCoord(i, 0),
                    localCoord(i, 3) - localCoord(i, 1), localCoord(i, 5) - localCoord(i, 1);

            J = P * S.inverse();
            Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);

            Eigen::Matrix2d U = svd.matrixU();
            Eigen::Matrix2d V = svd.matrixV();
            Eigen::Matrix2d R = U * V.transpose();

            if(R.determinant() < 0 ){
                V(0, 1) = -V(0, 1);
                V(1, 1) = -V(1, 1);
                R = U * V.transpose();
            }
            rotation_matrices[i] = R;
        }
    }

        void Parameterizer::solveGlobalStep() {
            Eigen::VectorXd bu{};
            Eigen::VectorXd bv{};
            for (int i = 0; i < n_faces; ++i) {
                auto face = mesh->polyface(i);
                auto fv_iter = mesh->fv_iter(face);
                auto i0 = (*fv_iter)->index();
                ++fv_iter;
                auto i1 = (*fv_iter)->index();
                ++fv_iter;
                auto i2 = (*fv_iter)->index();

                auto he2 = face->halfEdge();
                auto he0 = he2->next();
                auto he1 = he0->next();

                // note v0 is original in localspace
                Eigen::Vector2d e0, e1, e2;
                e0 << localCoord(i, 2) , localCoord(i, 3);
                e1 << localCoord(i, 4) - localCoord(i, 2), localCoord(i, 5) - localCoord(i, 3);
                e2 << -localCoord(i, 4), -localCoord(i, 5);

                Eigen::Vector2d b0 = cotangent_weights[he0->index()] * rotation_matrices[i] * e0;
                Eigen::Vector2d b1 = cotangent_weights[he1->index()] * rotation_matrices[i] * e1;
                Eigen::Vector2d b2 = cotangent_weights[he2->index()] * rotation_matrices[i] * e2;
                bu[i0] -= b0[0];
                bv[i0] -= b0[1];
                bu[i1] += b0[0];
                bv[i1] += b0[1];
                bu[i1] -= b1[0];
                bv[i1] -= b1[1];
                bu[i2] += b1[0];
                bv[i2] += b1[1];
                bu[i2] -= b2[0];
                bv[i2] -= b2[1];
                bu[i0] += b2[0];
                bv[i0] += b2[1];
        }
            iterParameterization.col(0) = solver.solve(bu);
            iterParameterization.col(1) = solver.solve(bv);
    }

    void Parameterizer::localGlobalArapIterations(int num_iterations) {
        rotation_matrices.resize(n_faces);

        for (int iter = 0; iter < num_iterations; ++iter) {
            // Step 1: Local step: Solve for rotation matrices
            solveLocalStep();
            // Step 2: Global step: Solve for new UV positions
            solveGlobalStep();
        }
    }


