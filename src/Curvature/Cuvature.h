#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <execution>  // For parallel execution policies
#include <mutex>      // For thread-safe file writing

#include "PolyMesh\IOManager.h"
#include "BaseFunc.h"
#define PI 3.14159265359

using namespace acamcad;
using namespace polymesh;

namespace Jade::Curvature {

    struct CurvatureInfo {
        double meanCurvature{0.0f};
        double gsCurvature{0.0f};

    };


    void CalLocalAveRegion(const std::shared_ptr<PolyMesh const> &mesh,
                           std::vector<double> &vertexLocalAveRegion);

    void CalMeanCurvature(const std::shared_ptr<PolyMesh const> &mesh,
                          std::vector<double> &vertexLocalAveRegion,
                          std::vector<CurvatureInfo> &curvatures);

    void CalGaussianCurvature(const std::shared_ptr<PolyMesh const> &mesh,
                              std::vector<double> &vertexLocalAveRegion,
                              std::vector<CurvatureInfo> &curvatures);
}