#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm> // For std::remove_if
#include "structs.h"


// Helper function to trim whitespace and carriage returns from strings
std::string trim_string(const std::string &str);
std::vector<Fem::Node> load_nodes(std::string file_name);
std::vector<Fem::Element> load_quad_elements(std::string file_name);
std::vector<Fem::Node> load_bc(std::string file_name, std::vector<Fem::Node>& nodes);
Fem::GlobalData load_configuration(std::string file_name);
Fem::GlobalData init_solution(std::string filename);
void print_config(Fem::GlobalData configuration);
