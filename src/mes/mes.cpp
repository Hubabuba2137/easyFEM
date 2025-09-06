#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <algorithm>

#include "mes.h"
#include "../gauss/gauss.h"

namespace Fem {

float dN1dxi(float eta) { return -0.25 * (1 - eta); }
float dN2dxi(float eta) { return 0.25 * (1 - eta); }
float dN3dxi(float eta) { return 0.25 * (1 + eta); }
float dN4dxi(float eta) { return -0.25 * (1 + eta); }

float dN1deta(float xi) { return -0.25 * (1 - xi); }
float dN2deta(float xi) { return -0.25 * (1 + xi); }
float dN3deta(float xi) { return 0.25 * (1 + xi); }
float dN4deta(float xi) { return 0.25 * (1 - xi); }

float N1(float xi, float eta) { return 0.25 * ((1 - xi) * (1 - eta)); }
float N2(float xi, float eta) { return 0.25 * ((1 + xi) * (1 - eta)); }
float N3(float xi, float eta) { return 0.25 * ((1 + xi) * (1 + eta)); }
float N4(float xi, float eta) { return 0.25 * ((1 - xi) * (1 + eta)); }

std::vector<double> pc_xi = {-1.0/std::sqrt(3), 1.0/std::sqrt(3), 1.0/std::sqrt(3), -1.0/std::sqrt(3)};
std::vector<double> pc_eta = {-1.0/std::sqrt(3), -1.0/std::sqrt(3), 1.0/std::sqrt(3), 1.0/std::sqrt(3)};
std::vector<double> pc_weights = {1,1,1,1};

std::vector<double> bc_xi = {-1.0/std::sqrt(3), 1.0/std::sqrt(3), 1,1, 1.0/std::sqrt(3), -1.0/std::sqrt(3), -1, -1};
std::vector<double> bc_eta = {-1,-1, -1.0/std::sqrt(3), 1.0/std::sqrt(3), 1, 1, 1.0/std::sqrt(3), -1.0/std::sqrt(3)};
std::vector<double> bc_weights = {1,1,1,1,1,1,1,1};

void showProgress(int current, int max)
{
    if (max <= 0) return;
    
    float percentage = (float)current / max * 100.0f;
    
    const int barWidth = 60;
    int filledWidth = (int)(percentage / 100.0f * barWidth);
    
    std::string bar = "[";
    for (int i = 0; i < barWidth; i++) {
        if (i < filledWidth) {
            bar += "="; 
        } else {
            bar += " ";  
        }
    }
    bar += "]";
    
    std::cout << "\r" << bar << " " 
              << std::fixed << std::setprecision(1) << percentage << "% "
              << "(" << current<<"/"<<max << ")" << std::flush;
    

    if (current >= max) {
        std::cout << std::endl;
    }
}

void showProgress_time(int current, int max, float elapsedTime) {
    if (max <= 0 || current <= 0) return;
    
    float percentage = (float)current / max * 100.0f;
    
    const int barWidth = 60;
    int filledWidth = (int)(percentage / 100.0f * barWidth);
    
    std::string bar = "[";
    for (int i = 0; i < barWidth; i++) {
        if (i < filledWidth) {
            bar += "=";
        } else {
            bar += " ";
        }
    }
    bar += "]";
    
    float avgTimePerOperation = elapsedTime / current;
    float remainingOperations = max - current;
    float estimatedTimeLeft = avgTimePerOperation * remainingOperations;
    
    std::cout << "\r" << bar << " "
              << std::fixed << std::setprecision(1) << percentage << "% "
              << "(ETA " << std::setprecision(0) << estimatedTimeLeft << "ms)" << std::flush;
    
    if (current >= max) {
        std::cout << std::endl;
    }
}


std::vector<Node> load_nodes(std::string file_name)
{
    std::fstream file;
    std::vector<Node> result;

    file.open(file_name);
    if(!file.good()){
        std::cout << "Couldn't load text file\n";
        return result;
    }

    std::string line;
    bool node_selection = false;

    while (std::getline(file, line)) {
        if (line == "*Node") {
            node_selection = true;
            continue;
        }
        
        if (line == "*Element" || line == "*BC") {
            node_selection = false;
        }

        if (node_selection) {
            int id;
            float x, y;
            char comma;
            std::istringstream iss(line);
            if (iss >> id >> comma >> x >> comma >> y) {
                result.emplace_back(x, y);
            }
        }
    }

    file.close();
    return result;
}

std::vector<Element> load_quad_elements(std::string file_name)
{
    std::fstream file;
    std::vector<Fem::Element> result;

    file.open(file_name);
    if(!file.good()){
        std::cout << "Couldn't load text file\n";
        return result;
    }

    std::string line;
    bool element_selection = false;

    while (std::getline(file, line)) {
        if (line == "*Element") {
            element_selection = true;
            continue;
        }
        
        if (line == "*Node" || line == "*BC") {
            element_selection = false;
        }

        if (element_selection) {
            int id;
            int node_id1, node_id2, node_id3, node_id4;
            char comma;
            std::istringstream iss(line);
            if (iss >> id >> comma >> node_id1 >> comma >> node_id2 >> comma >> node_id3 >> comma >> node_id4) {
                result.emplace_back(id, node_id1, node_id2, node_id3, node_id4);
            }
        }
    }

    file.close();

    return result;
}

GlobalData load_configuration(std::string file_name)
{
    std::fstream file;
    GlobalData data;

    file.open(file_name);
    if(!file.good()){
        std::cout << "Couldn't load text file\n";
        return data;
    }

    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        float value;

        if (iss >> key >> value) {
            if (key == "SimulationTime") data.total_time = value;
            else if (key == "SimulationStepTime") data.time_step = value;
            else if (key == "Conductivity") data.conductivity = value;
            //else if (key == "Alfa") data.alfa = value;
            //else if (key == "Tot") data.tot = value;
            else if (key == "InitialTemp") data.init_temperature = value;
            else if (key == "Density") data.density = value;
            else if (key == "SpecificHeat") data.specific_heat = value;
            else if (key == "Nodes_number") data.node_number = value;
            else if (key == "Elements_number") data.elem_number = value;
        }
    }

    file.close();
    return data;
}

std::vector<Node> set_bc(std::string file_name, std::vector<Node> &nodes)
{
    std::fstream file;

    file.open(file_name);
    if(!file.good()){
        std::cout << "Couldn't load text file: " << file_name << "\n";
        return nodes;
    }

    std::string line;
    bool bc_selection = false;

    while (std::getline(file, line)) {
        if (line == "*BC") {
            bc_selection = true;
            continue;
        }

        if (line == "*Node" || line == "*Element") {
            bc_selection = false;
        }

        if (bc_selection) {
            int node_id;
            float flux, alfa, t_ext;
            char comma;
            std::istringstream iss(line);
            if (iss >> node_id >> comma >> flux >> comma >> alfa >> comma >> t_ext) {
                //std::cout << "Parsed BC: Node ID=" << node_id << ", Flux=" << flux << ", Alfa=" << alfa << ", T_ext=" << t_ext << "\n"; // Debug output
                BC_node temp_bc(node_id, flux, alfa, t_ext);
                nodes[node_id-1].bc = temp_bc;
            } 
            /*
            else {
                std::cout << "Failed to parse line: " << line << "\n"; // Debug output
            }*/
        }
    }

    file.close();
    return nodes;
}

void print_config(GlobalData configuration)
{
    std::cout<<"-----------------Configuration-----------------\n";
    std::cout<<"Total time = "<<configuration.total_time<<"\n";
    std::cout<<"Time step = "<<configuration.time_step<<"\n";
    std::cout<<"Conductivity = "<<configuration.conductivity<<"\n";
    //std::cout<<"Alfa = "<<configuration.alfa<<"\n";
    //std::cout<<"Tot = "<<configuration.tot<<"\n";
    std::cout<<"Initial temperature = "<<configuration.init_temperature<<"\n";
    std::cout<<"Density = "<<configuration.density<<"\n";
    std::cout<<"Specific heat = "<<configuration.specific_heat<<"\n";
    std::cout<<"Number of nodes = "<<configuration.node_number<<"\n";
    std::cout<<"Number of elements = "<<configuration.elem_number<<"\n";
    std::cout<<"-----------------Configuration-----------------\n";
}

float distance(Node A, Node B)
{
    return std::sqrt(pow(A.x - B.x,2) + pow(A.y - B.y,2));
}

Matrix jacobian_mat(Element &element, std::vector<Node> &nodes, float pc_xi, float pc_eta)
{
    Matrix jacobian(2,2);

    jacobian[0][0] = dN1dxi(pc_eta)*nodes[element.node_ids[0]-1].x+Fem::dN2dxi(pc_eta)*nodes[element.node_ids[1]-1].x+Fem::dN3dxi(pc_eta)*nodes[element.node_ids[2]-1].x+Fem::dN4dxi(pc_eta)*nodes[element.node_ids[3]-1].x; //dxdxi
    jacobian[0][1] = dN1dxi(pc_eta)*nodes[element.node_ids[0]-1].y+Fem::dN2dxi(pc_eta)*nodes[element.node_ids[1]-1].y+Fem::dN3dxi(pc_eta)*nodes[element.node_ids[2]-1].y+Fem::dN4dxi(pc_eta)*nodes[element.node_ids[3]-1].y; //dydxi
    jacobian[1][0] = dN1deta(pc_xi)*nodes[element.node_ids[0]-1].x+Fem::dN2deta(pc_xi)*nodes[element.node_ids[1]-1].x+Fem::dN3deta(pc_xi)*nodes[element.node_ids[2]-1].x+Fem::dN4deta(pc_xi)*nodes[element.node_ids[3]-1].x; //dxdeta
    jacobian[1][1] = dN1deta(pc_xi)*nodes[element.node_ids[0]-1].y+Fem::dN2deta(pc_xi)*nodes[element.node_ids[1]-1].y+Fem::dN3deta(pc_xi)*nodes[element.node_ids[2]-1].y+Fem::dN4deta(pc_xi)*nodes[element.node_ids[3]-1].y; //dydeta
    
    //std::cout<<Fem::dN1deta(pc_xi)<<"*"<<nodes[element.node_ids[0]-1].x<<"+"<<Fem::dN1deta(pc_xi)<<"*"<<nodes[element.node_ids[1]-1].x<<"+"<<Fem::dN1deta(pc_xi)<<"*"<<nodes[element.node_ids[2]-1].x<<"+"<<Fem::dN1deta(pc_xi)<<"*"<<nodes[element.node_ids[3]-1].x<<"\n";

    return jacobian;
}

float det_jacobian(Fem::Matrix jacobian_mat)
{
    return jacobian_mat[0][0]*jacobian_mat[1][1] - jacobian_mat[1][0]*jacobian_mat[0][1];
}

Matrix inv_jacobian_mat(Fem::Matrix jacobian_mat)
{
    Fem::Matrix inv_jacobian(2,2);
    float det_J = det_jacobian(jacobian_mat);
    inv_jacobian[0][0] = jacobian_mat[1][1]/det_J;
    inv_jacobian[0][1] = -jacobian_mat[0][1]/det_J;
    inv_jacobian[1][0] = -jacobian_mat[1][0]/det_J;
    inv_jacobian[1][1] = jacobian_mat[0][0]/det_J;

    return inv_jacobian;
}

Matrix calc_local_H(Fem::Element &element, std::vector<Fem::Node> &nodes, float conductivity)
{
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

Matrix calc_local_Hbc(Fem::Element &element, std::vector<Fem::Node> &nodes)
{
    Fem::Matrix H_bc(4,4);

    for(int edge=0; edge<4; edge++){
        Fem::Matrix Hbc_edge(4,4);
        int id1 = element.node_ids[edge]-1;
        int id2 = element.node_ids[(edge+1)%4]-1;
        
        if(!(nodes[id1].bc.exist && nodes[id2].bc.exist)){
            continue;
        }

        float alfa = nodes[id1].bc.alfa; //alfa z W.B. Robina

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

Matrix calc_P(Fem::Element &element, std::vector<Fem::Node> &nodes)
{
    Fem::Matrix p_vec(4,1);

    for(int edge=0; edge<4; edge++){
        Fem::Matrix P_edge(4,1);
        int id1 = element.node_ids[edge]-1;
        int id2 = element.node_ids[(edge+1)%4]-1;

        if(!(nodes[id1].bc.exist && nodes[id2].bc.exist)){
            continue;
        }

        //W.B. Robina
        float temperature = nodes[id1].bc.t_ext;
        float alfa = nodes[id1].bc.alfa;

        //W.B. Neumanna;
        float flux = nodes[id1].bc.flux;

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

        //W.B. Robina
        for(int i=0; i<4;i++){
            p_vec[i][0] += P_edge[i][0]*det_J*alfa;
        }

        
        //W.B. Neumanna
        for(int i=0; i<4;i++){
            p_vec[i][0] += P_edge[i][0]*det_J*flux;
        }
    }

    return p_vec;
}

Matrix calc_local_C(Fem::Element &element, std::vector<Fem::Node> &nodes, float rho, float c)
{
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

void aggregate(Fem::Matrix &Global, Fem::Element element, Fem::Matrix &Local)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            int glob_i = element.node_ids[i]-1;
            int glob_j = element.node_ids[j]-1;

            Global[glob_i][glob_j] += Local[i][j];
        }
    }
}

void aggregate_p_vec(Fem::Matrix &P_vec, Fem::Element element, Fem::Matrix &Local)
{
    for(int i=0; i<4; i++){
        int glob_i = element.node_ids[i]-1;
        P_vec[glob_i][0] += Local[i][0];
    }
}

void write_to_vtu_file(int step, const std::vector<Fem::Node> &nodes, const std::vector<double> &temp, const std::vector<Fem::Element> &elements)
{

    std::ostringstream fname;
    fname << "sol_" << step << ".vtu";
    const std::string path = "../Data/"+fname.str();

    std::ofstream out(path);
    if (!out.is_open()) {
        std::perror(("Nie można otworzyć pliku " + path).c_str());
        return;
    }

    out << R"(<?xml version="1.0"?>)"
        << "\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           "  <UnstructuredGrid>\n"
           "    <Piece NumberOfPoints=\"" << nodes.size()
        << "\" NumberOfCells=\"" << elements.size() << "\">\n";

    // --- 4) Points ---
    out << "      <Points>\n"
           "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << std::setprecision(6);
    for (auto& n : nodes) {
        out << "          " << n.x << " " << n.y << " 0\n";
    }
    out << "        </DataArray>\n"
           "      </Points>\n";

    // --- 5) Cells ---
    out << "      <Cells>\n"
           "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (auto& e : elements) {
        out << "          "
            << (e.node_ids[0] - 1) << " "
            << (e.node_ids[1] - 1) << " "
            << (e.node_ids[2] - 1) << " "
            << (e.node_ids[3] - 1) << "\n";
    }
    out << "        </DataArray>\n"
           "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (size_t i = 0; i < elements.size(); ++i) {
        offset += 4; 
        out << "          " << offset << "\n";
    }
    out << "        </DataArray>\n"
           "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < elements.size(); ++i) {
        out << "          9\n";
    }
    out << "        </DataArray>\n"
           "      </Cells>\n";

    // --- 6) PointData: temperatura ---
    out << "      <PointData Scalars=\"Temperature\">\n"
           "        <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n";
    for (double t : temp) {
        out << "          " << t << "\n";
    }
    out << "        </DataArray>\n"
           "      </PointData>\n";

    out << "    </Piece>\n"
           "  </UnstructuredGrid>\n"
           "</VTKFile>\n";

    out.close();
}

void Fem::Solution::solve(bool write_to_vtu, bool write_temp_in_time)
{
    std::vector<double> t0(conf.node_number, conf.init_temperature);
    std::vector<double> t1 = t0;

    Fem::Matrix H_C(conf.node_number, conf.node_number);
    for(int i=0; i<conf.node_number; i++){
        for(int j=0; j<conf.node_number;j++){
            H_C[i][j] = this->Global_H[i][j]+this->Global_C[i][j]/conf.time_step;
        }
    }

    std::chrono::steady_clock::time_point for_beg = std::chrono::steady_clock::now();

    std::vector<temp_data> temperature_data;

    std::cout<<"Beginning time integration...\n";
    for(double time = conf.time_step; time <= conf.total_time; time+=conf.time_step){

        //showProgress(time, conf.total_time);
        for(int row=0; row<conf.node_number; row++){
            t1[row] = 0;
            for(int col=0; col<conf.node_number; col++){
                t1[row] += (this->Global_C[row][col] / conf.time_step) * t0[col];
            }
            t1[row] += this->Global_P[row][0];
        }

        t1 = Gauss(H_C, t1);

        std::string time_string = "Temperature at time " + std::to_string(time) + "s\n";
        std::string val_string = "\tMIN: " + std::to_string(*std::min_element(t1.begin(), t1.end())) + " MAX: " + std::to_string(*std::max_element(t1.begin(), t1.end())) + "\n";

        temperature_data.emplace_back(temp_data(time_string, val_string));

        if(write_to_vtu){
            write_to_vtu_file((int)time, this->nodes, t1, this->elements);
        }
        
        t0=t1;

        std::chrono::steady_clock::time_point progress_bar = std::chrono::steady_clock::now();
        float elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(progress_bar - for_beg).count();

        showProgress_time(time, conf.total_time, elapsed_time);
    }

    if(write_temp_in_time){
        std::cout<<"\n";
        for(auto&it:temperature_data){
            std::cout<<it.temp_time;
            std::cout<<it.temp_value;
        }
    }

    std::chrono::steady_clock::time_point for_end = std::chrono::steady_clock::now();
    std::cout << "Integration time = " << std::chrono::duration_cast<std::chrono::milliseconds>(for_end - for_beg).count() << "[ms]" << std::endl;
}
}
