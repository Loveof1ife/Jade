#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <cmath>

#include "PolyMesh\IOManager.h"
#include "BaseParameter.h"
#include "BaseFunc.h"

using namespace acamcad;
using namespace polymesh;

namespace Jade::Filter{
    void BilateralNormalFiltering(const std::shared_ptr<PolyMesh>& in_mesh,
                                  std::vector<MVector3>& filterNormalField,
                                  std::vector<MPoint3>& faceCenter,
                                  std::vector<double>& faceArea);

    void fitVertexToNormal(std::vector<MPoint3>& vertexPositions,
                           const std::shared_ptr<PolyMesh>& in_mesh,
                           const std::vector<MVector3>& filterNormalField,
                           const std::vector<MPoint3>& faceCenter,
                           const std::vector<double>& faceArea,
                           int maxIterations = 100,
                           double convergenceThreshold = 1e-6);
}


