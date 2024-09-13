#include "Algorithm.h"
#include <queue>

namespace Jade::ShortPath {

        const int INF = std::numeric_limits<int>::max();
        const int one_ring_dist = 1;

        struct Compare {
            bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b) {
                return a.second > b.second;  // Min-heap based on distance
            }
        };

        //to fine the shorest path from vid1 to vid2
        void Dijkstra(std::unique_ptr<PolyMesh> &mesh, int vid1, int vid2, std::vector<int>& vertexPath) {
            auto N = mesh->numVertices();

            //std::unordered_map<int, int> distance;
            //std::unordered_map<int, int> preNode;
            //std::unordered_map<int, bool> inS;

            std::vector<bool> inS(N, false);    // Set of processed vertices
            std::vector<int> distance(N, INF);  // Distance to each vertex
            std::vector<int> preNode(N, -1);

            std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, Compare> pq;

            distance[vid1] = 0;
            preNode[vid1] = vid1;
            pq.emplace(vid1, 0);

            while (!pq.empty()) {
                auto [vid, tmpDist] = pq.top();
                pq.pop();

                if (inS[vid]) continue;
                if (vid == vid2) break;

                inS[vid] = true;

                //update the 1-ring dist
                for (auto he: mesh->vertAdjacentHalfEdge(mesh->vert(vid))) {
                    auto adj_vid = he->toVertex()->index();
                    if (vid == adj_vid) {
                        std::cout << "skip a vertex ";
                        continue;
                    }

                    if (distance[vid] + one_ring_dist < distance[adj_vid]) {
                        distance[adj_vid] = distance[vid] + one_ring_dist;
                        preNode[adj_vid] = vid;
                        pq.emplace(adj_vid, distance[adj_vid]);
                    }
                }
            }

            //search path backward
            std::vector<int > edgePath;
            int current{vid2};

            if (distance[vid2] == INF) {
                std::cerr << "No path found!" << std::endl;
                return;
            } else {
                while (preNode[current] != vid1) {
                    vertexPath.push_back(preNode[current]);
                    edgePath.push_back(mesh->edgeBetween(mesh->vert(current), mesh->vert(preNode[current]))->index() );
                    current = preNode[current];
                }
                vertexPath.push_back(vid1);
                edgePath.push_back(mesh->edgeBetween(mesh->vert(current), mesh->vert(vid1))->index() );

                std::reverse(vertexPath.begin(), vertexPath.end());
                std::reverse(edgePath.begin(), edgePath.end());
            }
        }

        void DijkstraGroup(std::unique_ptr<PolyMesh> &mesh, const std::vector<int> &landmarks,
                           std::vector<std::vector<int> > &vertexPaths) {
            auto N = landmarks.size();
            vertexPaths.clear();
            vertexPaths.resize(N - 1);

            for (size_t i = 0; i < N - 1; ++i) {
                int sp = landmarks[i];
                int ep = landmarks[i + 1];
                if(sp == ep) continue;
                Dijkstra(mesh, sp, ep, vertexPaths[i]);
            }
        }

        void SteinerTree(std::unique_ptr<PolyMesh> &mesh, const std::vector<int> &terminal_ids,
                         std::vector<std::pair<int, int>> &tree_edges) {
            size_t N = terminal_ids.size();

            //using 2D matrix represent graph
            std::vector<std::vector<int> > dists(N, std::vector<int>(N, INF));

            //the innermost vector is the indice sequence of path nodes
            std::vector<std::vector<std::vector<int> > > paths(N, std::vector<std::vector<int> >(N));

            for (size_t i = 0; i < N; i++) {
                for (size_t j = i + 1; j < N; j++) {
                    Dijkstra(mesh, terminal_ids[i], terminal_ids[j], paths[i][j]);
                    dists[i][j] = dists[j][i] = static_cast<int>(paths[i][j].size());
                }
            }
            std::vector<bool> inMST(N, false);
            std::vector<int> minEdge(N, INF);
            std::vector<int> parent(N, -1);  // Store the MST
            minEdge[0] = 0;

            for (size_t i = 0; i < N; ++i) {
                int shortestId{-1};
                int shortestEdge{INF};

                //step1: Find the vertex with the smallest edge not in MST
                for (size_t j = 0; j < N; ++j) {
                    if (!inMST[j] && minEdge[j] < shortestEdge) {
                        shortestEdge = minEdge[j];
                        shortestId = j;
                    }

                    // If no valid vertex is found, the graph is disconnected
                    if (minEdge[shortestId] == std::numeric_limits<int>::max()) {
                        std::cerr << "No MST found!" << std::endl;
                        return;
                    }

                    // Step 3.2: Add the selected vertex to the MST
                    inMST[shortestId] = true;

                    if (parent[shortestId] != -1) {
                        size_t u = std::min(shortestId, parent[shortestId]);
                        size_t v = std::max(shortestId, parent[shortestId]);

                        for (size_t k = 0; k < paths[u][v].size() - 1; k++) {
                            tree_edges.emplace_back(paths[u][v][k], paths[u][v][k + 1]);
                        }
                    }

                    // Step 3.3: Update the minimum edges for the vertices not in MST
                    for (size_t i = 0; i < N; ++i) {
                        if (inMST[i]) continue;
                        else {
                            int newPathLength =
                                    i >= shortestId ? paths[i][shortestId].size() : paths[shortestId][i].size();
                            if (newPathLength < minEdge[i]) {
                                minEdge[i] = newPathLength;
                                parent[i] = shortestId;
                            }
                        }
                    }
                }
            }
        }
    }


    //void Dijkstra(std::unique_ptr<PolyMesh>& mesh, int vid1, int vid2, std::vector<int>& path)
//{
//    size_t N = mesh->numVertices();
//    std::vector<bool> inS(N, false);
//    std::vector<int> dis(N, std::numeric_limits<int>::infinity());
//    std::vector<int> pre(N, -1);
//    dis[vid1] = 0;
//    pre[vid1] = vid1;

//    while (true) {
//        int min_dist = std::numeric_limits<int>::infinity();
//        int u = -1;
//        // 更新相邻顶点的距离
//        for (auto edge : mesh->vertAdjacentEdge(mesh.vert(u))) {
//            int v = edge->getVert(0)->index();
//            if (v == u) v = edge->getVert(1)->index();
//            auto dis_uv = edge->length();
//            if (!inS[v] && dis[u] + dis_uv < dis[v])
//                dis[v] = dis[u] + dis_uv;
//            pre[v] = u;
//        }
//        // 找到距离最近的未处理顶点
//        for (size_t v = 0; v < N; ++v) {
//            if (!inS[v] && dis[v] < min_dist) {
//                min_dist = dis[v];
//                u = v;
//            }
//        }
//        if (u == -1 || min_dist == std::numeric_limits<int>::infinity())
//            break;
//        else
//            inS[u] = true;
//    }
//    std::vector<MEdge* > Epath{ 0 };

//    int current = vid2;
//    while (current != vid1) {
//        if (current == -1) {
//            std::cerr << "No path found!" << std::endl;
//            return;
//        }
//        path.push_back(vid2);
//        Epath.push_back(mesh.edgeBetween(mesh.vert(current), mesh.vert(pre[current])));
//        current = pre[current];
//    }
//    path.push_back(vid1);
//    std::reverse(path.begin(), path.end());
//    std::reverse(Epath.begin(), Epath.end());

//}
        
//
//    void SteinerTree(std::unique_ptr<PolyMesh>& mesh, const std::vector<int>& terminal_ids, std::vector<std::pair<int, int>>& tree_edges) {
//        size_t N = terminal_ids.size();
//
//        std::vector<std::vector<std::vector<int >>> paths(N, std::vector<std::vector<int>>(N));
//        std::vector<std::vector<int>> dist(N, std::vector<int>(N, std::numeric_limits<int>::max()));
//
//        //MinTotalGraph
//        for (size_t i = 0; i < N; ++i) {
//            for (size_t j = i + 1; j < N; ++j) {
//                Dijkstra(mesh, terminal_ids[i], terminal_ids[j], paths[i][j]);
//                dist[i][j] = dist[j][i] = paths[i][j].size() - 1;  // assuming edge weights are 1 for simplicity
//            }
//        }
//
//        //construct MinSpanTree from MinTotalGraph
//        std::vector<bool> inMST(N, false);
//        std::vector<int> minEdge(N, std::numeric_limits<int>::max());
//        std::vector<int> parent(N, -1);
//        minEdge[0] = 0;
//
//        for (size_t i = 0; i < N; ++i) {
//            int u = -1;
//            int shortest = std::numeric_limits<int>::max();
//
//            // 找到未加入最小生成树的节点中，具有最小边的节点
//            for (size_t j = 0; j < N; ++j) {
//                if (!inMST[j] && (u == -1 || minEdge[j] < shortest)) {
//                    shortest = minEdge[j];
//                    u = j;
//                }
//            }
//
//            if (minEdge[u] == std::numeric_limits<int>::max()) {
//                std::cerr << "No MST found!" << std::endl;
//                return;
//            }
//
//            for (size_t k = 0; k < paths[u][parent[u]].size() - 1; ++k) {
//                tree_edges.emplace_back(paths[u][parent[u]][k], paths[u][parent[u]][k + 1]);
//            }
//
//            // 将节点u加入最小生成树
//            inMST[u] = true;
//
//            // 更新其他节点到最小生成树的最短边信息
//            //iter1：1 // n-1
//            //iter2: 2 // n-2 if dis_new[n-2] < dis_old[n-2]
//            //update
//            for (size_t v = 0; v < N; ++v) {
//                if (!inMST[v] && dist[u][v] < minEdge[v]) {
//                    minEdge[v] = dist[u][v];
//                    parent[v] = u;
//                }
//            }
//        }
//    }
//}
//}

