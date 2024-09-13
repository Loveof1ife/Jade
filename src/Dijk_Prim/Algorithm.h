#pragma once
#include "PolyMesh\IOManager.h"
#include <PolyMesh/PolyMesh.h>
#include <PolyMesh/PolyMesh_Base.h>

#include <unordered_map>
#include <algorithm>
#include <limits>
#include <iostream>
#include <memory>


namespace Jade::ShortPath{

    using PolyMesh = acamcad::polymesh::PolyMesh;

    void Dijkstra(std::unique_ptr<PolyMesh>& mesh, int vid1, int vid2, std::vector<int>& path);

    void DijkstraGroup(std::unique_ptr<PolyMesh>& mesh, const std::vector<int> &landmarks,  std::vector<std::vector<int>>& paths);

    void SteinerTree(std::unique_ptr<PolyMesh>& mesh, const std::vector<int> &terminal_ids, std::vector<std::pair<int, int> >& tree_edge);
}
