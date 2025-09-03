#pragma once

/*

#include <vector>
#include "matrix/matrix.h"

namespace Fem{
    extern float dN1dxi(float eta);
    extern float dN2dxi(float eta);
    extern float dN3dxi(float eta);
    extern float dN4dxi(float eta);

    extern float dN1deta(float xi);
    extern float dN2deta(float xi);
    extern float dN3deta(float xi);
    extern float dN4deta(float xi);

    extern float N1(float xi, float eta);
    extern float N2(float xi, float eta);
    extern float N3(float xi, float eta);
    extern float N4(float xi, float eta);

    extern std::vector<double> pc_xi;
    extern std::vector<double> pc_eta;
    extern std::vector<double> pc_weights;

    extern std::vector<double> bc_xi;
    extern std::vector<double> bc_eta;
    extern std::vector<double> bc_weights;

    void showProgress(int current, int max);

    struct BC_node{
        int id;

        bool exist=false;

        float flux;
        float alfa;
        float t_ext;

        BC_node(){
            this->id=0;
            this->flux=0.0;
            this->alfa=0.0;
            this->t_ext=0.0;
        }

        BC_node(int id, float flux, float alfa, float t_ext){
            this->exist=true;
            this->id=id;
            this->flux=flux;
            this->alfa=alfa;
            this->t_ext=t_ext;
        }
    };
    
    struct Node{
        float x;
        float y;
        
        BC_node bc;
        
        Node(){
            this-> x = 0; 
            this-> y = 0; 
        }
        
        Node(float x, float y){
            this->x = x;
            this->y = y;
        }
    };
    
    struct Element{
        int id;
        int node_ids[4];
        
        Matrix H_local;
        Matrix H_bc;
        Matrix P;
        Matrix C;
        
        Element(int id, int n1, int n2, int n3, int n4) : H_local(3, 3), H_bc(3,3), P(3,1), C(3,3) {
            this->id = id;
            this->node_ids[0] = n1;
            this->node_ids[1] = n2;
            this->node_ids[2] = n3;
            this->node_ids[3] = n4;
        }
    };

    std::vector<Fem::Node> load_nodes(std::string file_name);
    std::vector<Fem::Element> load_quad_elements(std::string file_name);
    std::vector<Fem::Node> set_bc(std::string file_name, std::vector<Fem::Node>& nodes);

    float distance(Fem::Node A, Fem::Node B);
    Fem::Matrix jacobian_mat(Fem::Element &element, std::vector<Fem::Node> &nodes, float pc_xi, float pc_eta);
    float det_jacobian(Fem::Matrix jacobian_mat);
    Fem::Matrix inv_jacobian_mat(Fem::Matrix jacobian_mat);
    
    Fem::Matrix calc_local_H(Fem::Element &element, std::vector<Fem::Node> &nodes, float conductivity);
    Fem::Matrix calc_local_Hbc(Fem::Element &element, std::vector<Fem::Node> &nodes);
    Fem::Matrix calc_P(Fem::Element &element, std::vector<Fem::Node> &nodes);
    Fem::Matrix calc_local_C(Fem::Element &element, std::vector<Fem::Node> &nodes, float rho, float c);
    
    void aggregate(Fem::Matrix &Global, Fem::Element element, Fem::Matrix &Local);
    void aggregate_p_vec(Fem::Matrix &P_vec, Fem::Element element, Fem::Matrix &Local);
    
    void write_to_vtu_file(int step, const std::vector<Fem::Node>& nodes, const std::vector<double>& temp, const std::vector<Fem::Element>& elements);
    
    struct GlobalData{
        float total_time;
        float time_step;
        float conductivity;
        //float alfa;
        //float tot;
        float init_temperature;
        float density;
        float specific_heat;
        int node_number;
        int elem_number;
        
        GlobalData(){
            this->total_time=0;
            this->time_step=0;
            this->conductivity=25;
            //this->alfa=0;
            //this->tot=0;
            this->init_temperature=0;
            this->density=7800;
            this->specific_heat=700;
            this->node_number=0;
            this->elem_number=0;
        }
    };
    
    GlobalData load_configuration(std::string file_name);
    void print_config(GlobalData configuration);

    struct Solution{
        Matrix Global_H;
        Matrix Global_C;
        Matrix Global_P;
        GlobalData conf;
        std::vector<Fem::Node> nodes;
        std::vector<Fem::Element> elements;

        Solution(Matrix &Global_H, Matrix &Global_C, Matrix &Global_P, GlobalData &conf, std::vector<Node> &nodes, std::vector<Element> &elements): Global_H(4,4), Global_C(4,4), Global_P(4,4){
            this->Global_H = Global_H;
            this->Global_C = Global_C;
            this->Global_P = Global_P;
            this->conf = conf;
            this->nodes = nodes;
            this->elements = elements;
        }

        Solution(std::string filename): Global_H(4,4), Global_C(4,4), Global_P(4,1){
            this->nodes = load_nodes(filename);

            this->nodes = set_bc(filename, this->nodes);
            
            this->elements = load_quad_elements(filename);
            this->conf = load_configuration(filename);

            this->conf.node_number = this->nodes.size();
            this->conf.elem_number= this->elements.size();

            print_config(this->conf);

            this->Global_H = Fem::Matrix(conf.node_number,conf.node_number);
            this->Global_C = Fem::Matrix(conf.node_number,conf.node_number);

            this->Global_P = Fem::Matrix(conf.node_number,1);

            //tworzenie macierzy dla każdego elementu
            std::cout<<"Assembling elemental matrices...\n";
            int i =0;
            int max_iter = this->elements.size();
            for(Fem::Element &element: this->elements){
                
                element.H_local = calc_local_H(element, this->nodes, this->conf.conductivity);
                element.H_bc = calc_local_Hbc(element, this->nodes);
                element.P = calc_P(element, this->nodes);
                element.C = calc_local_C(element, this->nodes, this->conf.density, this->conf.specific_heat);
                
                //sumowanie H_l i H_bc
                for(int row=0; row<4; row++){
                    for(int col=0; col<4; col++){
                        element.H_local[row][col] += element.H_bc[row][col];
                    }
                }
                
                i++;
                showProgress(i, max_iter);
            }

            //agreagacja
            std::cout<<"Agregating matrices...\n";
            i=0;
            for(Fem::Element &it:elements){
                aggregate(this->Global_H, it, it.H_local);
                aggregate(this->Global_C, it, it.C);
                aggregate_p_vec(this->Global_P, it, it.P);
                i++;
                showProgress(i, max_iter);
            }
        }
            
        }

        void solve(bool write_to_vtu=false, bool write_temp_in_time = true);
    };

    struct temp_data{
        std::string temp_time;
        std::string temp_value;

        temp_data(){
            this->temp_time = "";
            this->temp_value = "";
        }

        temp_data(std::string time, std::string val){
            this->temp_time = time;
            this->temp_value = val;
        }
    };
}*/