#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "PolyMesh\IOManager.h"
#include <string>

#define PI 3.14159265359

using namespace acamcad;
using namespace polymesh;

void cal_local_ave_region(const PolyMesh *mesh, std::vector<double> &vertexLAR);
void cal_mean_curvature(const PolyMesh *mesh, std::vector<double> &vertexLAR);
void cal_gaussian_curvature(const PolyMesh *mesh, std::vector<double> &vertexLAR);

