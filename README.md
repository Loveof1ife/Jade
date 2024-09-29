
# Discrete Geometry Processing Project


## Introduction

This project focuses on discrete geometry processing algorithms，based on the USTC Digital Geometry Processing course and Games102。


## Algorithm 1: Dijkstra's Algorithm in Mesh


### Problem Description:

the mesh can be viewed as a graph, where the vertices of the mesh correspond to the nodes in the graph and the mesh topology corresponds to the connections of the graph.

- **input： mesh, id of 2 vertices**
- **output:  the shortest path on the mesh connecting the input vertices**

### Transformation:

1. **Mesh topology -> Graph connections**:
   The mesh's topology defines the neighboring relations between its vertices. This is represented as edges connecting the nodes in the graph.
2. **Mesh vertices -> Graph nodes**:
   Each vertex in the mesh corresponds to a node in the graph, with the edge weight typically representing the geometric distance or some other metric between vertices.
3. **Shortest path on the mesh -> Shortest path on a graph without negative weights**:
   In this case, the edges of the graph have no negative weights, making it suitable to use the classic **Dijkstra's algorithm** to find the shortest path.

### Data Structures:

- **d(v)**: The distance from the starting point `sp` to vertex `v`.
- **inS(v)**: Boolean indicating whether vertex `v` is included in the set `S` (the set of processed vertices).
- **N**: The total number of vertices in the grid `M`.
- **Length(u, v)**: The length of the edge connecting vertices `u` and `v`.
- **H**: A heap data structure used to efficiently find the next vertex to process.
- **pre(v)**: Stores the previous vertex in the shortest path leading to vertex `v`.


### Pseudocode:

```pseudo
//initialize
for i in [1, N] do
     d(i) ← +∞, inS(i) ← false
end
d(sp) ← 0, pre(sp) ← sp
H.push((sp, d(sp))

//main loop 
while H ≠ ø do
    tmp ← H.top()
    if inS(tmp) = true then
        continue
    end
    if tmp = ep then
        break
    else
        inS(tmp) ← true
        for each uv ∈ 1-ring of tmp do
            if d(tmp) + Length(tmp, v) < d(v) then
                d(vv) ← d(tmp) + Length(tmp, v)
                H.push((uv, d(vv)))
                pre(uv) ← tmp
            end
        end
    end
end

//search path backward
while ep ≠ pre(ep) do
    Epath.push((ep, pre(ep)))
    Vpath.push(ep)
    ep ← pre(ep)
end
Vpath.push(ep)

```
## Algorithm 2: Steiner's Tree

### Problem Description:

1. Computing the shortest paths between all pairs of vertices, resulting in a complete graph of the input vertices over the mesh.
2. Finding the MST of this complete graph, which serves as an approximate solution to Steiner's Tree.

**Input**: 
- `M`: Input mesh.
- `Lmk`: Input set of landmarks.

**Output**: 
- `MST`: Approximate Steiner's Tree.

**Data**:
- `complete_graph`: A heap data structure.
- `spanning_tree_current`: A disjoint-set (Union-Find) structure.

#### Steps

1. **Initialize and Compute the Complete Graph**:
     - For each landmark `i`:
          - Calculate the modified Dijkstra algorithm to find the shortest paths from `Lmk[i]` to all other landmarks.
          - Store the paths in the `complete_graph` heap.

     ```pseudo
     Lsize ← Lmk.size()
     for i ∈ [1, Lsize] do
          count ← Lsize - i
          # Modified Dijkstra Algorithm #
          while count ≠ 0 do
               if inS(tmp) = true then
                    continue
               end
               if tmp = Lmk[u] and u > i then
                    count ← count - 1
               end
               inS(tmp) ← true
          end
          
          for j ∈ [i + 1, Lsize] do
               calculate Vpath
               complete_graph.push(Vpath)
          end
     end
     ```

2. **Compute the Minimum Spanning Tree**:
   - Initialize the spanning tree.
   - While there are still vertices to process, extract the top element of the `complete_graph` and add the edge to the MST if it does not form a cycle.


    ```pseudo
    index ← 0
    for i ∈ [1, Lsize] do
        spanning_tree_current(i) ← i
    end
    
    while index < Lsize do
        E ← complete_graph.top()
        m, n ← the two endpoints of E
        m_id, n_id ← m, n in Lmk
        
        if m_id < n_id then
            swap(m_id, n_id)
        end
        
        if m_id ≠ n_id then
            MST.push(E)
            index ← index + 1
            for i ∈ [1, Lsize] do
                if spanning_tree_current(i) = n_id then
                    spanning_tree_current(i) ← m_id
                end
            end
        end
    end
    ```

---

## Algorithm 3: Curvature
### Problem Description:

- **vertexLocalAveRegion**: computes the average region for each vertex based on the areas of the adjacent faces. If an obtuse triangle exists, the area is divided accordingly to its vertices.
- **mean curvature**: calculates the mean curvature for each vertex using the Laplace-Beltrami operator. The cotangent of angles between edges in the vertex's 1-ring neighborhood is used in the calculation.
- **Gaussian curvature**: The angular deficit, which is the difference between 2π and the sum of angles at a vertex, is divided by the average area associated with the vertex to determine the Gaussian curvature.

## Algorithm 4: Filter

### Problem Description:


- **M**: The input 3D mesh with faces and vertices.
- **σ**: Parameter controlling the influence of spatial proximity (distance between neighboring faces).
- **ω**: Parameter controlling the influence of normal similarity (difference between normals of neighboring faces).
- **normal_iter**: Number of iterations to perform for normal smoothing.
- **vertex_iter**: Number of iterations to perform for updating vertex positions.

### Outputs:

- **normal field**
- **updated vertex**

### Process:

1. **Initialization of Normals**:
   - The normals of the mesh faces are initialized based on the original face normals.

2. **Bilateral Filtering on raw Normals**:
   - For each face in the mesh, its normal is updated using the normals of its neighboring faces, weighted by both the spatial distance between face centers and the normal differences.

3. **Vertex Position Update based on filterd normals**:
   - After filtering the normals, the positions of the vertices are adjusted based on the smoothed normals of the neighboring faces, iteratively updating the vertex positions(Gauss-Seidel iteration per vertex per iteration , other vertex fixed).

4. **Iterative Process**:
   - The algorithm iteratively refines the normals for a specified number of iterations (`normal_iter`).
   - Update vertex positions one by one, Fixed the positions of the remaining vertices. The reason for this is to make the loss function a quadratic function of one variable
     $$
     x_{\text{new}, i}= x_i + \sum_s{j \in \text{neighborhood}(x_i)} n_j \cdot( n_j^T (c_j - x_i))
     $$


## Algorithm 5: Tuttes's Embedding
Given a triangulated surface homeomorphic to a disk, if the(u,v)coordinates at the boundary vertices lie on a convex polygon inorder, and if the coordinates of the internal vertices are a convex combination of their neighbors, then the (u, v) coordinates form a valid parameterization (without self-intersections, bijective)

### Problem Description:

- **boundary**: Artificially put the boundary point UV On a convex polygon
- **interior**: Uniform laplacian, Linear system

### Process:
1. Read into a grid, identify the boundaries, place them in order on a polygon (bisected N bisected circles),
2. The inner points are convex combinations in it's one-ring, N points with N equations