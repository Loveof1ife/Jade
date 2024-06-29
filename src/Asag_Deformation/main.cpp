#include "Asag_deformation.h"

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cout << "========== Hw6 Usage  ==========\n";
        std::cout << std::endl;
        std::cout
                << "Input:	ACAM_mesh_HW6.exe	[model]	[fix vertex id] [move vertex id and position]\n";  //fix vertex id: number + vertex id ;move vertex id and position: number + id position
        std::cout << "Input:	ACAM_mesh_HW6.exe	mesh.obj	fix_handle.txt	move_handle.txt\n";
        std::cout << std::endl;
        std::cout << "=================================================\n";
        return 0;
    }
    PolyMesh* pmesh = new PolyMesh();
    //¶ÁÈëÍø¸ñ
    std::string mesh_path = argv[1];
    loadMesh(mesh_path, pmesh);

    //loadMesh("Bunny_head.obj",&mesh);
    pmesh->updateMeshNormal();
    int nf =pmesh->numPolygons();
    int nv = pmesh->numVertices();

    //position backup
    std::vector<MVector3> pos_mesh_ref;
    Build_Origin_Positon(pmesh, pos_mesh_ref, nv);

    std::vector<double> cots;
    Calc_Cot_Weight(cots);

    //set fix handle
    //set<int> handles_f = {12,505,381,712,296};
    std::set<int> handles_f;
    std::ifstream fix_handle(argv[2]);
    int fix_num;
    fix_handle >> fix_num;
    for (int i = 0; i < fix_num; i++) {
        int vertex_tmp;
        fix_handle >> vertex_tmp;
        handles_f.emplace(vertex_tmp);
    }

    //set move handle
    //vector<int> handles_m = {652};
    //vector<MVector3> handles_m_pos = { MVector3(0.05,0.2,0.05) };
    std::vector<int> handles_m;
    std::vector<MVector3> handles_m_pos;
    std::ifstream move_handle(argv[3]);
    int move_num;
    move_handle >> move_num;
    for (int i = 0; i < move_num; i++)
    {
        int vertex_tmp;
        move_handle >> vertex_tmp;
        double p1, p2, p3;
        move_handle >> p1;
        move_handle >> p2;
        move_handle >> p3;
        handles_m.push_back(vertex_tmp);
        handles_m_pos.push_back(MVector3(p1, p2, p3));
    }

    std::set<int> handles = handles_f;
    handles.insert(handles_m.begin(), handles_m.end());

    //local-global iteration
    std::vector<Eigen::Matrix3d> R_list;
    R_list.resize(nv);

    Eigen::SparseMatrix<double> lp_mat;
    lp_mat.resize(nv, nv);
    Cal_Lp_Mat(pmesh, lp_mat, cots, handles, nv);

    for (int iter = 0; iter < 10; iter++)
    {
        Global_Step(pmesh, R_list, pos_mesh_ref, cots, nv);
        Local_Step(pmesh, lp_mat, R_list, pos_mesh_ref, cots, handles_f, handles_m, handles_m_pos, nv);
    }

    writeMesh("deformation.obj", pmesh);

}



