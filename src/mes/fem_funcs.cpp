#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <chrono>

#include "structs.h"
#include "matrix.h"

#include "load_function.cpp"

float distance(Fem::Node A, Fem::Node B){
    return std::sqrt(pow(A.x - B.x,2) + pow(A.y - B.y,2));
}

Fem::Matrix jacobian_mat(Fem::Element &element, std::vector<Fem::Node> &nodes, float pc_xi, float pc_eta){
    Fem::Matrix jacobian(2,2);

    jacobian[0][0] = Fem::dN1dxi(pc_eta)*nodes[element.node_ids[0]-1].x+Fem::dN2dxi(pc_eta)*nodes[element.node_ids[1]-1].x+Fem::dN3dxi(pc_eta)*nodes[element.node_ids[2]-1].x+Fem::dN4dxi(pc_eta)*nodes[element.node_ids[3]-1].x; //dxdxi
    jacobian[0][1] = Fem::dN1dxi(pc_eta)*nodes[element.node_ids[0]-1].y+Fem::dN2dxi(pc_eta)*nodes[element.node_ids[1]-1].y+Fem::dN3dxi(pc_eta)*nodes[element.node_ids[2]-1].y+Fem::dN4dxi(pc_eta)*nodes[element.node_ids[3]-1].y; //dydxi
    jacobian[1][0] = Fem::dN1deta(pc_xi)*nodes[element.node_ids[0]-1].x+Fem::dN2deta(pc_xi)*nodes[element.node_ids[1]-1].x+Fem::dN3deta(pc_xi)*nodes[element.node_ids[2]-1].x+Fem::dN4deta(pc_xi)*nodes[element.node_ids[3]-1].x; //dxdeta
    jacobian[1][1] = Fem::dN1deta(pc_xi)*nodes[element.node_ids[0]-1].y+Fem::dN2deta(pc_xi)*nodes[element.node_ids[1]-1].y+Fem::dN3deta(pc_xi)*nodes[element.node_ids[2]-1].y+Fem::dN4deta(pc_xi)*nodes[element.node_ids[3]-1].y; //dydeta


    
    //std::cout<<Fem::dN1deta(pc_xi)<<"*"<<nodes[element.node_ids[0]-1].x<<"+"<<Fem::dN1deta(pc_xi)<<"*"<<nodes[element.node_ids[1]-1].x<<"+"<<Fem::dN1deta(pc_xi)<<"*"<<nodes[element.node_ids[2]-1].x<<"+"<<Fem::dN1deta(pc_xi)<<"*"<<nodes[element.node_ids[3]-1].x<<"\n";

    return jacobian;
}

float det_jacobian(Fem::Matrix jacobian_mat){
    return jacobian_mat[0][0]*jacobian_mat[1][1] - jacobian_mat[1][0]*jacobian_mat[0][1];
}

Fem::Matrix inv_jacobian_mat(Fem::Matrix jacobian_mat){
    Fem::Matrix inv_jacobian(2,2);
    float det_J = det_jacobian(jacobian_mat);
    inv_jacobian[0][0] = jacobian_mat[1][1]/det_J;
    inv_jacobian[0][1] = -jacobian_mat[0][1]/det_J;
    inv_jacobian[1][0] = -jacobian_mat[1][0]/det_J;
    inv_jacobian[1][1] = jacobian_mat[0][0]/det_J;

    return inv_jacobian;
}

Fem::Matrix calc_local_H(Fem::Element &element, std::vector<Fem::Node> &nodes, float conductivity){
    Fem::Matrix H(4,4);

    for(int i=0; i<4; i++){
        Fem::Matrix H_local(4,4);
        Fem::Matrix jacob = jacobian_mat(element, nodes, Fem::pc_xi[i], Fem::pc_eta[i]);
        float det_J = det_jacobian(jacob);

        Fem::Matrix inv_jacob = inv_jacobian_mat(jacob);

        Fem::Matrix dNdx(4,1);
        Fem::Matrix dNdy(4,1);

        dNdx[0][0] = inv_jacob[0][0]*Fem::dN1dxi(Fem::pc_eta[i]) + inv_jacob[0][1]*Fem::dN1deta(Fem::pc_xi[i]);
        dNdx[1][0] = inv_jacob[0][0]*Fem::dN2dxi(Fem::pc_eta[i]) + inv_jacob[0][1]*Fem::dN2deta(Fem::pc_xi[i]);
        dNdx[2][0] = inv_jacob[0][0]*Fem::dN3dxi(Fem::pc_eta[i]) + inv_jacob[0][1]*Fem::dN3deta(Fem::pc_xi[i]);
        dNdx[3][0] = inv_jacob[0][0]*Fem::dN4dxi(Fem::pc_eta[i]) + inv_jacob[0][1]*Fem::dN4deta(Fem::pc_xi[i]);
        
        dNdy[0][0] = inv_jacob[1][0]*Fem::dN1dxi(Fem::pc_eta[i]) + inv_jacob[1][1]*Fem::dN1deta(Fem::pc_xi[i]);
        dNdy[1][0] = inv_jacob[1][0]*Fem::dN2dxi(Fem::pc_eta[i]) + inv_jacob[1][1]*Fem::dN2deta(Fem::pc_xi[i]);
        dNdy[2][0] = inv_jacob[1][0]*Fem::dN3dxi(Fem::pc_eta[i]) + inv_jacob[1][1]*Fem::dN3deta(Fem::pc_xi[i]);
        dNdy[3][0] = inv_jacob[1][0]*Fem::dN4dxi(Fem::pc_eta[i]) + inv_jacob[1][1]*Fem::dN4deta(Fem::pc_xi[i]);

        Fem::Matrix Bx = dNdx*dNdx.transpose();
        Fem::Matrix By = dNdy*dNdy.transpose();

        for(int row=0; row<4; row++){
            for(int col=0; col<4; col++){
                H_local[row][col] = (Bx[row][col] + By[row][col])*conductivity*det_J;
            }
        }

        for(int row=0; row<4; row++){
            for(int col=0; col<4; col++){
                H[row][col] += H_local[row][col]*Fem::pc_weights[i];
            }
        }
    }

    return H;
}     

Fem::Matrix calc_local_Hbc(Fem::Element &element, std::vector<Fem::Node> &nodes, float alfa){
    Fem::Matrix H_bc(4,4);

    for(int edge=0; edge<4; edge++){
        Fem::Matrix Hbc_edge(4,4);
        int id1 = element.node_ids[edge]-1;
        int id2 = element.node_ids[(edge+1)%4]-1;

        
        if(!(nodes[id1].bc && nodes[id2].bc)){
            continue;
        }

        float det_J = 0.5*distance(nodes[id1], nodes[id2]);

        for(int pc=0; pc<2; pc++){
            int pc_on_egde = edge*2+pc;
            Fem::Matrix Ns(4,1);

            Ns[0][0] = Fem::N1(Fem::bc_xi[pc_on_egde], Fem::bc_eta[pc_on_egde]);
            Ns[1][0] = Fem::N2(Fem::bc_xi[pc_on_egde], Fem::bc_eta[pc_on_egde]);
            Ns[2][0] = Fem::N3(Fem::bc_xi[pc_on_egde], Fem::bc_eta[pc_on_egde]);
            Ns[3][0] = Fem::N4(Fem::bc_xi[pc_on_egde], Fem::bc_eta[pc_on_egde]);

            Fem::Matrix BTB = Ns*Ns.transpose();

            for(int i=0; i<4;i++){
                for(int j=0; j<4; j++){
                    Hbc_edge[i][j] += BTB[i][j]*Fem::bc_weights[pc_on_egde];
                }
            }
        }

        for(int i=0; i<4;i++){
            for(int j=0; j<4; j++){
                H_bc[i][j] += Hbc_edge[i][j]*det_J*alfa;
            }
        }
    }

    return H_bc;
}

Fem::Matrix calc_P(Fem::Element &element, std::vector<Fem::Node> &nodes, float alfa, float temperature){
    Fem::Matrix p_vec(4,1);

    for(int edge=0; edge<4; edge++){
        Fem::Matrix P_edge(4,1);
        int id1 = element.node_ids[edge]-1;
        int id2 = element.node_ids[(edge+1)%4]-1;

        
        if(!(nodes[id1].bc && nodes[id2].bc)){
            continue;
        }

        float det_J = 0.5*distance(nodes[id1], nodes[id2]);

        for(int pc=0; pc<2; pc++){
            int pc_on_egde = edge*2+pc;
            Fem::Matrix Ns(4,1);

            Ns[0][0] = Fem::N1(Fem::bc_xi[pc_on_egde], Fem::bc_eta[pc_on_egde]);
            Ns[1][0] = Fem::N2(Fem::bc_xi[pc_on_egde], Fem::bc_eta[pc_on_egde]);
            Ns[2][0] = Fem::N3(Fem::bc_xi[pc_on_egde], Fem::bc_eta[pc_on_egde]);
            Ns[3][0] = Fem::N4(Fem::bc_xi[pc_on_egde], Fem::bc_eta[pc_on_egde]);

            for(int i=0; i<4;i++){
                P_edge[i][0] += Ns[i][0]*Fem::bc_weights[pc_on_egde]*temperature;
            }
        }

        for(int i=0; i<4;i++){
            p_vec[i][0] += P_edge[i][0]*det_J*alfa;
        }
    }

    return p_vec;
}

Fem::Matrix calc_local_C(Fem::Element &element, std::vector<Fem::Node> &nodes, float rho, float c){
    Fem::Matrix C(4,4);

    for(int i=0; i<4; i++){
        Fem::Matrix C_local(4,4);
        Fem::Matrix jacob = jacobian_mat(element, nodes, Fem::pc_xi[i], Fem::pc_eta[i]);
        float det_J = det_jacobian(jacob);

        Fem::Matrix inv_jacob = inv_jacobian_mat(jacob);

        Fem::Matrix base_functions_values(4,1);

        base_functions_values[0][0] = Fem::N1(Fem::pc_xi[i], Fem::pc_eta[i]);
        base_functions_values[1][0] = Fem::N2(Fem::pc_xi[i], Fem::pc_eta[i]);
        base_functions_values[2][0] = Fem::N3(Fem::pc_xi[i], Fem::pc_eta[i]);
        base_functions_values[3][0] = Fem::N4(Fem::pc_xi[i], Fem::pc_eta[i]);

        Fem::Matrix B_mat = base_functions_values*base_functions_values.transpose();

        for(int row=0; row<4; row++){
            for(int col=0; col<4; col++){
                C_local[row][col] = B_mat[row][col]*rho*c*det_J;
            }
        }

        for(int row=0; row<4; row++){
            for(int col=0; col<4; col++){
                C[row][col] += C_local[row][col]*Fem::pc_weights[i];
            }
        }
    }

    return C;
} 

void aggregate(Fem::Matrix &Global, Fem::Element element, Fem::Matrix &Local){
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            int glob_i = element.node_ids[i]-1;
            int glob_j = element.node_ids[j]-1;

            Global[glob_i][glob_j] += Local[i][j];
        }
    }
}

void aggregate_p_vec(Fem::Matrix &P_vec, Fem::Element element, Fem::Matrix &Local){
    for(int i=0; i<4; i++){
        int glob_i = element.node_ids[i]-1;
        P_vec[glob_i][0] += Local[i][0];
    }
}

Fem::Solution prepare_sol(std::string filename){
    std::chrono::steady_clock::time_point prep_beg = std::chrono::steady_clock::now();
    
    Fem::GlobalData conf = init_solution(filename);
    std::cout<<"Data loaded form file "<<filename<<"\n";

    std::vector<Fem::Node> nodes = conf.nodes;
    std::vector<Fem::Element> elements = conf.elements;

    conf.node_number = nodes.size();
    conf.elem_number = elements.size();

    print_config(conf);

    Fem::Matrix Global_H(conf.node_number,conf.node_number);
    Fem::Matrix Global_C(conf.node_number,conf.node_number);
    Fem::Matrix Global_P(conf.node_number,1);

    //tworzenie macierzy dla każdego elementu
    for(Fem::Element &it: elements){
        it.H_local = calc_local_H(it, nodes, conf.conductivity);
        it.H_bc = calc_local_Hbc(it, nodes, conf.alfa);
        it.P = calc_P(it, nodes, conf.alfa, conf.tot);
        it.C = calc_local_C(it, nodes, conf.density, conf.specific_heat);

        //std::cout<<it.C<<"\n";
        //sumowanie H_l i H_bc
        for(int row=0; row<4; row++){
            for(int col=0; col<4; col++){
                it.H_local[row][col] += it.H_bc[row][col];
            }
        }
    }

    std::cout<<"Created elemt matrices...\n";

    //agreagacja
    for(Fem::Element &it:elements){
        aggregate(Global_H, it, it.H_local);
        aggregate(Global_C, it, it.C);
        aggregate_p_vec(Global_P, it, it.P);
    }

    //std::cout<<Global_P<<"\n";

    std::cout<<"Agregated global matrices...\n";
    std::chrono::steady_clock::time_point prep_end = std::chrono::steady_clock::now();
    std::cout << "Prepping time = " << std::chrono::duration_cast<std::chrono::milliseconds>(prep_end - prep_beg).count() << "[ms]" << std::endl;

    return Fem::Solution(Global_H, Global_C, Global_P, conf, nodes, elements);
}
