#include "imgui.h"
#include "raylib.h"
#include <vector>

#include "../include/geometry.h"
#include "mes/structs.h"

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
    if(conf.total_time<=0){
        std::cout<<"total_time<=0\n";
        return false;
    }
    else if(conf.time_step<=0){
        std::cout<<"time_step<=0\n";
        return false;
    }
    else if(conf.conductivity<=0){
        std::cout<<"conductivity<=0 \n";
        return false;}
    else if(conf.alfa<=0){
        std::cout<<"alfa<=0 \n";
        return false;}
    else if(conf.tot<-273.0){
        std::cout<<"tot<-273.0 \n";
        return false;}
    else if(conf.init_temperature<-273.0){
        std::cout<<"init_temperature<273.0 \n";
        return false;}
    else if(conf.density<0){
        std::cout<<"density<0 \n";
        return false;}
    else if(conf.specific_heat<0){
        std::cout<<"specific_heat<0\n";
        return false;}
    else if(conf.node_number<=0){
        std::cout<<"node_number<=0\n";
        return false;}
    else if(conf.elem_number<=0){
        std::cout<<"elem_number<=0 \n";
        return false;}
    else{
        return true;
    }
}
