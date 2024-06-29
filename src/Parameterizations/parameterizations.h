#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "PolyMesh/IOManager.h"
#include <string>

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

using namespace acamcad;
using namespace polymesh;

typedef Eigen::SparseMatrix<double> SMatrix;

void Tutte_Embedding(PolyMesh *mesh);
void Tutte_Embedding(Eigen::MatrixX2d& uv, const PolyMesh& mesh);
void Tutte_Embedding_Mvc(PolyMesh* );
void Local_Global_Approach(PolyMesh*);

