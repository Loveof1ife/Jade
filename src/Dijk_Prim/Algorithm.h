#pragma once
#include "PolyMesh\IOManager.h"
#include <PolyMesh/PolyMesh.h>
#include <PolyMesh/PolyMesh_Base.h>

using namespace acamcad;
using namespace polymesh;

namespace Algo
{
    void Dijkstra( PolyMesh &mesh, int vid1, int vid2, std::vector<int> &);

    void Dijkstra_group( PolyMesh &mesh, const std::vector<int> &landmarks,  std::vector<std::vector<int>> &);

    void SteinerTree(PolyMesh &mesh, const std::vector<int> &terminal_ids);
}