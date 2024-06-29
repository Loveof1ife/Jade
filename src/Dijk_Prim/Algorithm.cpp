#include <algorithm>
#include "Algorithm.h"

namespace Algo {
//to fine the shorest path from vid1 to vid2
void Dijkstra( PolyMesh& mesh, int vid1, int vid2, std::vector<int>& path)
{
    size_t N = mesh.numVertices();
    std::vector<bool> inS(N, false);
    std::vector<int> dis(N, std::numeric_limits<int>::infinity());
    std::vector<int> pre(N, -1);
    dis[vid1] = 0;
    pre[vid1] = vid1;

    while(true){
        int min_dist = std::numeric_limits<int>::infinity();
        int u = -1;
        // 更新相邻顶点的距离
        for(auto edge : mesh.vertAdjacentEdge(mesh.vert(u))){
            int v = edge->getVert(0)->index();
            if (v == u) v = edge->getVert(1)->index();
            auto dis_uv = edge->length();
            if(!inS[v] && dis[u] + dis_uv < dis[v])
                dis[v] = dis[u] + dis_uv;
            pre[v] = u;
        }
        // 找到距离最近的未处理顶点
        for(size_t v = 0; v < N; ++v ){
            if(!inS[v] && dis[v] < min_dist){
                min_dist = dis[v];
                u = v;
            }
        }
        if (u == -1 || min_dist == std::numeric_limits<int>::infinity())
            break;
        else
            inS[u] = true;
    }
    std::vector<MEdge* > Epath{0};

    int current = vid2;
    while(current != vid1){
        if (current == -1) {
            std::cerr << "No path found!" << std::endl;
            return;
        }
        path.push_back(vid2);
        Epath.push_back(mesh.edgeBetween(mesh.vert(current), mesh.vert(pre[current])));
        current = pre[current];
    }
    path.push_back(vid1);
    std::reverse(path.begin(), path.end());
    std::reverse(Epath.begin(), Epath.end());

}


void Dijkstra_group( PolyMesh& mesh, const std::vector<int>& landmarks, std::vector< std::vector<int> >& paths)
{
    auto N = landmarks.size();
    paths.clear();
    paths.resize(N - 1);

    for (size_t i = 0; i < N - 1; ++i) {
        int sp = landmarks[i];
        int ep = landmarks[i + 1];
        Algo::Dijkstra(mesh, sp, ep, paths[i]);
    }

}

void SteinerTree(PolyMesh &mesh, const std::vector<int> &terminal_ids, std::vector<std::pair<int, int>> &tree_edges){
    size_t N = terminal_ids.size();

    std::vector<std::vector<std::vector<int >>> paths(N, std::vector<std::vector<int>>(N));
    std::vector<std::vector<int>> dist(N, std::vector<int>(N, std::numeric_limits<int>::max()));

    //MinTotalGraph
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            Dijkstra(mesh, terminal_ids[i], terminal_ids[j], paths[i][j]);
            dist[i][j] = dist[j][i] = paths[i][j].size() - 1;  // assuming edge weights are 1 for simplicity
        }
    }

    //construct MinSpanTree from MinTotalGraph
    std::vector<bool> inMST(N, false);
    std::vector<int> minEdge(N, std::numeric_limits<int>::max());
    std::vector<int> parent(N, -1);
    minEdge[0] = 0;

    for (size_t i = 0; i < N; ++i) {
        int u = -1;
        int shortest = std::numeric_limits<int>::max();

        // 找到未加入最小生成树的节点中，具有最小边的节点
        for (size_t j = 0; j < N; ++j) {
            if (!inMST[j] && (u == -1 || minEdge[j] < shortest)) {
                shortest = minEdge[j];
                u = j;
            }
        }

        if (minEdge[u] == std::numeric_limits<int>::max()) {
            std::cerr << "No MST found!" << std::endl;
            return;
        }

        for (size_t k = 0; k < paths[u][parent[u]].size() - 1; ++k) {
            tree_edges.emplace_back(paths[u][parent[u]][k], paths[u][parent[u]][k + 1]);
        }

        // 将节点u加入最小生成树
        inMST[u] = true;

        // 更新其他节点到最小生成树的最短边信息
        //iter1：1 // n-1
        //iter2: 2 // n-2 if dis_new[n-2] < dis_old[n-2]
        //update
        for (size_t v = 0; v < N; ++v) {
            if (!inMST[v] && dist[u][v] < minEdge[v]) {
                minEdge[v] = dist[u][v];
                parent[v] = u;
            }
        }

    }



}

}


