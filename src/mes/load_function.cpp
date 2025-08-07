#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm> // For std::remove_if
#include "structs.h"
#include "load_function.h"

// Helper function to trim whitespace and carriage returns from strings
std::string trim_string(const std::string &str) {
    std::string result = str;
    // Remove carriage returns and whitespace from both ends
    result.erase(result.begin(), std::find_if(result.begin(), result.end(), 
        [](int ch) { return !std::isspace(ch) && ch != '\r'; }));
    result.erase(std::find_if(result.rbegin(), result.rend(), 
        [](int ch) { return !std::isspace(ch) && ch != '\r'; }).base(), result.end());
    return result;
}

std::vector<Fem::Node> load_nodes(std::string file_name) {
    std::fstream file;
    std::vector<Fem::Node> result;
    file.open(file_name);
    if (!file.good()) {
        std::cout << "Couldn't load text file\n";
        return result;
    }
    std::string line;
    bool node_selection = false;
    while (std::getline(file, line)) {
        // Trim the line to remove any carriage returns or whitespace
        line = trim_string(line);
        
        if (line == "*Node") {
            node_selection = true;
            continue;
        }
        
        if (line == "*Element" || line == "*BC") {
            node_selection = false;
        }
        if (node_selection && !line.empty()) {
            int id;
            float x, y;
            char comma;
            std::istringstream iss(line);
            // Replace commas with spaces for easier parsing
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss2(line);
            if (iss2 >> id >> x >> y) {
                result.emplace_back(x, y);
            }
        }
    }
    file.close();
    return result;
}

std::vector<Fem::Element> load_quad_elements(std::string file_name) {
    std::fstream file;
    std::vector<Fem::Element> result;
    file.open(file_name);
    if (!file.good()) {
        std::cout << "Couldn't load text file\n";
        return result;
    }
    std::string line;
    bool element_selection = false;
    while (std::getline(file, line)) {
        // Trim the line to remove any carriage returns or whitespace
        line = trim_string(line);
        
        if (line == "*Element") {
            element_selection = true;
            continue;
        }
        
        if (line == "*Node" || line == "*BC") {
            element_selection = false;
        }
        if (element_selection && !line.empty()) {
            int id;
            int node_id1, node_id2, node_id3, node_id4;
            // Replace commas with spaces for easier parsing
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);
            if (iss >> id >> node_id1 >> node_id2 >> node_id3 >> node_id4) {
                result.emplace_back(id, node_id1, node_id2, node_id3, node_id4);
            }
        }
    }
    file.close();
    return result;
}



std::vector<Fem::Node> load_bc(std::string file_name, std::vector<Fem::Node>& nodes) {
    std::fstream file;
    file.open(file_name);
    if (!file.good()) {
        std::cout << "Couldn't load text file\n";
        return nodes;
    }
    std::string line;
    bool bc_selection = false;
    while (std::getline(file, line)) {
        // Trim the line to remove any carriage returns or whitespace
        line = trim_string(line);
        
        if (line == "*BC") {
            bc_selection = true;
            continue;
        }
        if (line == "*Node" || line == "*Element") {
            bc_selection = false;
        }
        if (bc_selection && !line.empty()) {
            int node_id;
            std::istringstream iss(line);
            if (iss >> node_id) {
                if (node_id > 0 && node_id <= nodes.size()) {
                    nodes[node_id-1].bc = true;
                }
            }
        }
    }
    file.close();
    return nodes;
}

Fem::GlobalData load_configuration(std::string file_name) {
    std::fstream file;
    Fem::GlobalData data;
    file.open(file_name);
    if (!file.good()) {
        std::cout << "Couldn't load text file\n";
        return data;
    }
    std::string line;
    while (std::getline(file, line)) {
        // Trim the line to remove any carriage returns or whitespace
        line = trim_string(line);
        
        std::istringstream iss(line);
        std::string key;
        float value;
        if (iss >> key >> value) {
            if (key == "SimulationTime") data.total_time = value;
            else if (key == "SimulationStepTime") data.time_step = value;
            else if (key == "Conductivity") data.conductivity = value;
            else if (key == "Alfa") data.alfa = value;
            else if (key == "Tot") data.tot = value;
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

Fem::GlobalData init_solution(std::string filename){
    Fem::GlobalData configuration = load_configuration(filename);
    std::vector<Fem::Node> init_nodes = load_nodes(filename);
    configuration.elements = load_quad_elements(filename);
    
    configuration.nodes = load_bc(filename, init_nodes);

    return configuration;
}

// Prints current configuration from file
void print_config(Fem::GlobalData configuration) {
    std::cout << "-----------------Configuration-----------------\n";
    std::cout << "Total time = " << configuration.total_time << "\n";
    std::cout << "Time step = " << configuration.time_step << "\n";
    std::cout << "Conductivity = " << configuration.conductivity << "\n";
    std::cout << "Alfa = " << configuration.alfa << "\n";
    std::cout << "Tot = " << configuration.tot << "\n";
    std::cout << "Initial temperature = " << configuration.init_temperature << "\n";
    std::cout << "Density = " << configuration.density << "\n";
    std::cout << "Specific heat = " << configuration.specific_heat << "\n";
    std::cout << "Number of nodes = " << configuration.node_number << "\n";
    std::cout << "Number of elements = " << configuration.elem_number << "\n";
    std::cout << "-----------------Configuration-----------------\n";
}
