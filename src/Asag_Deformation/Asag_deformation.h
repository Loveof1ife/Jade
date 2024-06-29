#pragma once

#include <fstream>
#include"PolyMesh\IOManager.h"
#include"PolyMesh\PolyMesh.h"
#include<Eigen/SVD>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include <vector>
#include<set>
#include<iostream>
#include<math.h>
#include<algorithm>
#include<omp.h>
#include <fstream>
using namespace acamcad;
using namespace polymesh;

void Calc_Cot_Weight(std::vector<double>& cots);

void Build_Origin_Positon(PolyMesh* pmesh,
                          std::vector<MVector3>& pos_mesh_ref,
                          int nv);

void Cal_Lp_Mat(PolyMesh* pmesh,
                Eigen::SparseMatrix<double>& lp_mat,
                const std::vector<double>& cots,
                const std::set<int>& handles,
                const int nv);

void Global_Step(PolyMesh* pmesh,
                 std::vector<Eigen::Matrix3d>& R_list,
                 const std::vector<MVector3>& pos_mesh_ref,
                 const std::vector<double> cots,
                 const int nv);


void Local_Step(PolyMesh* pmesh,
                const Eigen::SparseMatrix<double>& lp_mat,
                const std::vector<Eigen::Matrix3d>& R_list,
                const std::vector<MVector3>& pos_mesh_ref,
                const std::vector<double>& cots,
                const std::set<int>& handles_f,
                const std::vector<int>& handles_m,
                const std::vector<MVector3>& handles_m_pos,
                const int& nv);