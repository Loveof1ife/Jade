#include "Asag_deformation.h"

void Calc_Cot_Weight(std::vector<double>& cots)
{
    cots.clear();
    cots.resize(mesh.numHalfEdges(), 0.);
    for (auto ithe = mesh.halfedge_begin(); ithe != mesh.halfedge_end(); ithe++)
    {
        auto he = *ithe;
        if (mesh.isBoundary(*ithe))
            continue;
        auto v0=he->fromVertex()->position();
        auto v1 = he->toVertex()->position();
        auto v2 = he->next()->toVertex()->position();
        auto e0 = v0 - v2;
        auto e1 = v1 - v2;
        double cotangle = dot(e0,e1) / cross(e0,e1).norm();
//		cots[ithe->idx()] = cotangle;
        cots[he->index()] = cotangle;
    }
}
void Build_Origin_Positon(PolyMesh* pmesh,
                          std::vector<MVector3>& pos_mesh_ref,
                          int nv)
{
    pos_mesh_ref.resize(nv);
    for (auto itv = mesh.vertices_begin(); itv !=mesh.vertices_end() ; itv++)
    {
        pos_mesh_ref[(*itv)->index()] = (*itv)->position();
    }
}

void Cal_Lp_Mat(PolyMesh* pmesh,
                Eigen::SparseMatrix<double>& lp_mat,
                const std::vector<double>& cots,
                const std::set<int>& handles,
                const int nv)

{
    std::vector<Eigen::Triplet<double> > triplets;
    triplets.resize(nv);
    for(size_t i =0; i < nv; i++){
        MVert* v_i = pmesh->vert(i);
        int vi_index = v_i->index();

        if (handles.count(vi_index) > 0){
            Eigen::Triplet(vi_index, vi_index, 1);
            continue;
        }
        double weight_i{0.0};
        for(auto voh_it = pmesh->voh_iter(v_i); voh_it.isValid(); voh_it++){
            MHalfedge* he = *voh_it;
            double weight_ij = (cots[he->index()] + cots[he->pair()->index()]) ;
            auto v_j = he->toVertex();
            int vj_index = v_j->index();
            triplets.emplace_back(vi_index, vj_index, -weight_ij);
            weight_i += weight_ij;
        }
        triplets.emplace_back(vi_index, vi_index, weight_i);
    }

    lp_mat.resize(nv, nv);
    lp_mat.setFromTriplets(triplets.begin(),triplets.end());

}

void Global_Step(PolyMesh* pmesh,
                 std::vector<Eigen::Matrix3d>& R_list,
                 const std::vector<MVector3>& pos_mesh_ref,
                 const std::vector<double> cots,
                 const int nv)

{

    for (int i = 0; i < nv; i++){
        MVert* v_i= mesh.vert(i);
        auto vi_index = v_i->index();
        Eigen::Matrix3d R_i;
        for(auto voh_it = pmesh->voh_iter(v_i); voh_it.isValid(); voh_it++){
            auto he = *voh_it;
            auto v_j = he->toVertex();
            int vj_index = v_j->index();
            auto e_ij = pos_mesh_ref[vi_index] -  pos_mesh_ref[vj_index];
            auto e_ij_d = v_i->position() - v_j->position();
            Eigen::Vector3d e_ij_{e_ij[0], e_ij[1], e_ij[2]};
            Eigen::Vector3d e_ij_d_{e_ij_d[0], e_ij_d[1], e_ij_d[2]};
            R_i += cots[he->index()] * (e_ij_ * e_ij_d_.transpose());
        }
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(R_i, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();

        if (R_i.determinant() < 0)
        {
            U(0, 2) *= -1;
            U(1, 2) *= -1;
            U(2, 2) *= -1;
            R_i = V * U.transpose();
        }
        R_list[vi_index] = R_i;
    }
}
void Local_Step(PolyMesh* pmesh,
                const Eigen::SparseMatrix<double>& lp_mat,
                const std::vector<Eigen::Matrix3d>& R_list,
                const std::vector<MVector3>& pos_mesh_ref,
                const std::vector<double>& cots,
                const std::set<int>& handles_f,
                const std::vector<int>& handles_m,
                const std::vector<MVector3>& handles_m_pos,
                const int& nv)
{
    Eigen::MatrixX3d  xyz, b;
    xyz.resize(nv, 3);
    b.resize(nv, 3);

    for (int i = 0; i < nv; i++)
    {
        auto v_i = mesh.vert(i);
        Eigen::Vector3d b_tmp(0., 0., 0.);
        int vi_index = v_i->index();

        if(handles_f.count(vi_index) > 0)
            continue;
        else if(std::find(handles_m.begin(), handles_m.end(), vi_index) != handles_m.end())
            continue;

        for (auto voh_it = mesh.voh_iter(v_i); voh_it.isValid(); ++voh_it)
        {
            MVert* v_j = (*voh_it)->toVertex();
            int vj_index = v_j->index();

            Eigen::Matrix3d R_tmp = R_list[vi_index] + R_list[vj_index];
            auto e_ij_ = pos_mesh_ref[vi_index] - pos_mesh_ref[vj_index];
            Eigen::Vector3d e_ij{e_ij_[0], e_ij_[1], e_ij_[2]};
            double weight{(cots[vj_index] + cots[vi_index]) / 2.0 };
            b_tmp += weight * R_tmp * e_ij;
        }
        b(i, 0) = b_tmp(0);
        b(i, 1) = b_tmp(1);
        b(i, 2) = b_tmp(2);

        //set handles
        for (int index_hf:handles_f)
        {
            auto b_tmp = pos_mesh_ref[index_hf];
            b(index_hf, 0) = b_tmp[0];
            b(index_hf, 1) = b_tmp[1];
            b(index_hf, 2) = b_tmp[2];
        }

        for (int index_hm = 0; index_hm < handles_m.size();index_hm++)
        {
            auto b_tmp = handles_m_pos[index_hm];
            b(handles_m[index_hm], 0) = b_tmp[0];
            b(handles_m[index_hm], 1) = b_tmp[1];
            b(handles_m[index_hm], 2) = b_tmp[2];
        }

    }

    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    solver.compute(lp_mat);

    xyz.col(0) = solver.solve(b.col(0));
    xyz.col(1) = solver.solve(b.col(1));
    xyz.col(2) = solver.solve(b.col(2));


    for (int i = 0; i < nv; i++)
    {
        MVert* v_i = mesh.vert(i);
        int vi_index = v_i->index();
        v_i->setPosition(xyz(vi_index, 0), xyz(vi_index, 1), xyz(vi_index, 2));
    }
}