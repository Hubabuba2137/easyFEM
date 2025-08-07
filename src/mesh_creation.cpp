#include "imgui.h"
#include "raylib.h"
#include <vector>

#include "../include/geometry.h"

#include "mes/"

void add_new_node(Vector2 worldPos, go::Vertex &mesh){
    bool can_add_node = true;
    for(auto&it:mesh.vertices){
        if(it.pos == worldPos){
            can_add_node = false;
            break;
        }
    }

    if(can_add_node){
        go::Node node_1(worldPos);
        mesh.add_vertex(node_1);
    }
}

void pop_node(go::Vertex &mesh){
    if(!mesh.vertices.empty()){
        mesh.vertices.pop_back();
        mesh.create_edges();
    }
}

bool check_if_params_correct(Fem::GlobalData conf){
    if(conf.)
}
