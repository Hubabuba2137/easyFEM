#include "../include/geometry.h"
#include "../include/to_fem.h"
#include "mes/mes.h"
#include "../include/meshing.h"

#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

/*
void to_FEM(std::vector<go::Node> &nodes, 
    std::vector<go::Node> &bc_nodes, 
    std::vector<go::Vertex> &mesh, 
    Fem::GlobalData &conf)
{
    std::string file_name = "../Data/fem_data.txt";
    std::ofstream file(file_name, std::ios::out);

    if(!file.is_open()){
        std::cerr << "Error: Could not open file " << file_name << " for writing.\n";
        return;
    }

    remove_duplicate_nodes(nodes);
    remove_duplicate_nodes(bc_nodes);

    file << "SimulationTime "      << conf.total_time      << "\n";
    file << "SimulationStepTime "  << conf.time_step       << "\n";
    file << "Conductivity "        << conf.conductivity    << "\n";
    file << "InitialTemp "         << conf.init_temperature<< "\n";
    file << "Density "             << conf.density         << "\n";
    file << "SpecificHeat "        << conf.specific_heat   << "\n";
    file << "Nodes_number "        << nodes.size()         << "\n";
    file << "Elements_number "     << mesh.size()          << "\n";

    file << "*Node\n";
    for(size_t i=0;i<nodes.size();++i)
        file << (i+1) << ", " << nodes[i].pos.x << ", " << nodes[i].pos.y << "\n";

    file << "*Element\n";
    for(size_t i=0;i<mesh.size();++i){

        std::vector<go::Node> v = { mesh[i].vertices[0], mesh[i].vertices[1], mesh[i].vertices[2], mesh[i].vertices[3] };
        v = order_ccw_nodes(v);         
        if(signed_area(v) < 0) std::reverse(v.begin(), v.end()); 

        int ids[4];
        for(int k=0;k<4;k++){
            ids[k] = find_node_id(v[k], nodes);

            if(ids[k] < 0){
                std::cerr << "[to_FEM] Element " << (i+1) << " ma nieznany węzeł; sprawdź tolerancję dopasowania.\n";
                ids[k] = 0; 
            }
        }
        file << (i+1) << ", " << (ids[0]+1) << ", " << (ids[1]+1) << ", " << (ids[2]+1) << ", " << (ids[3]+1) << "\n";

    }

    file << "*BC\n";
    for(auto& b : bc_nodes){
        int bc_id = find_node_id(b, nodes);
        if(bc_id >= 0){
            float neumann = 0.f;
            float alfa = 300.f;
            float t_ext = 1200.f;

            file << (bc_id+1) << ", " << neumann << ", " << alfa << ", " << t_ext << "\n";
        }
    }

    std::cout << "File data written\n";
    file.close();
}

*/