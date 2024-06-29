#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "PolyMesh\IOManager.h"
#include <string>

#define CenterSigma  1
#define NormalSigma  1
#define IteratorTimes 1000


using namespace acamcad;
using namespace polymesh;

void Bilateral_Normal_Filtering(const PolyMesh*  mesh);
