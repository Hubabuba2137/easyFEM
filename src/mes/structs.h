#pragma once

#include <vector>
#include "matrix.cpp"


namespace Fem{
    float dN1dxi(float eta){return -0.25*(1-eta);};
    float dN2dxi(float eta){return 0.25*(1-eta);};
    float dN3dxi(float eta){return 0.25*(1+eta);};
    float dN4dxi(float eta){return -0.25*(1+eta);};

    float dN1deta(float xi){return -0.25*(1-xi);};
    float dN2deta(float xi){return -0.25*(1+xi);};
    float dN3deta(float xi){return 0.25*(1+xi);};
    float dN4deta(float xi){return 0.25*(1-xi);};

    float N1(float xi, float eta){return 0.25*((1-xi)*(1-eta));}
    float N2(float xi, float eta){return 0.25*((1+xi)*(1-eta));}
    float N3(float xi, float eta){return 0.25*((1+xi)*(1+eta));}
    float N4(float xi, float eta){return 0.25*((1-xi)*(1+eta));}

    std::vector<double> pc_xi = {-1.0/std::sqrt(3), 1.0/std::sqrt(3), 1.0/std::sqrt(3), -1.0/std::sqrt(3)};
    std::vector<double> pc_eta = {-1.0/std::sqrt(3), -1.0/std::sqrt(3), 1.0/std::sqrt(3), 1.0/std::sqrt(3)};
    std::vector<double> pc_weights = {1,1,1,1};

    std::vector<double> bc_xi = {-1.0/std::sqrt(3), 1.0/std::sqrt(3), 1,1, 1.0/std::sqrt(3), -1.0/std::sqrt(3), -1, -1};
    std::vector<double> bc_eta = {-1,-1, -1.0/std::sqrt(3), 1.0/std::sqrt(3), 1, 1, 1.0/std::sqrt(3), -1.0/std::sqrt(3)};
    std::vector<double> bc_weights = {1,1,1,1,1,1,1,1};

    struct Robin_bc{
        float alfa;
        float t_ext;
    };
    struct Neumann_bc{
        float q_n;
    };

    struct Bc_type{
        Robin_bc robin;
        Neumann_bc neumann;
    };

    struct Node{
        float x;
        float y;

        bool bc;

        Node(){
        this-> x = 0; 
        this-> y = 0; 
        }

        Node(float x, float y){
            this->x = x;
            this->y = y;
            this->bc = false;
        }
    };

    struct Element{
        int id;
        int node_ids[4];

        Matrix H_local;
        Matrix H_bc;
        Matrix P;
        Matrix C;
        
        Element(): H_local(4,4), H_bc(4,4), P(4,1), C(4,4){
          
        }


        Element(int id, int n1, int n2, int n3, int n4) : H_local(4, 4), H_bc(4,4), P(4,1), C(4,4) {
            this->id = id;
            this->node_ids[0] = n1;
            this->node_ids[1] = n2;
            this->node_ids[2] = n3;
            this->node_ids[3] = n4;
        }
    };


    struct Triangle{
        int id;
        int node_ids[3];

        Triangle(int id, int n1, int n2, int n3){
            this->node_ids[0] = n1;
            this->node_ids[1] = n2;
            this->node_ids[2] = n3;
            this->id = id;
        }
    };

    struct GlobalData{
        float total_time;
        float time_step;
        float conductivity;
        float alfa;
        float tot;
        float init_temperature;
        float density;
        float specific_heat;
        int node_number;
        int elem_number;
        
        std::vector<Node> nodes;
        std::vector<Element> elements;

        GlobalData(){
            this->total_time=500;
            this->time_step=50;
            this->conductivity=25;
            this->alfa=300;
            this->tot=1200;
            this->init_temperature=100;
            this->density=7800;
            this->specific_heat=700;
            this->node_number=0;
            this->elem_number=0;
        }


    };

    struct Solution{
        Matrix Global_H;
        Matrix Global_C;
        Matrix Global_P;
        GlobalData conf;

        std::vector<Node> nodes;
        std::vector<Element> elements;
        
        Solution(): Global_H(4,4), Global_C(4,4), Global_P(4,4){

        }

        Solution(Matrix &Global_H, Matrix &Global_C, Matrix &Global_P, GlobalData &conf, std::vector<Node> &nodes, std::vector<Element> &elements): Global_H(4,4), Global_C(4,4), Global_P(4,4){
            this->Global_H = Global_H;
            this->Global_C = Global_C;
            this->Global_P = Global_P;
            this->conf = conf;
            this->nodes = nodes;
            this->elements = elements;
        }

        Solution(GlobalData data): Global_H(4,4), Global_C(4,4), Global_P(4,4){
            this->conf=data;
            this->nodes = data.nodes;
            

            this->elements = data.elements;
        }
    };
}
