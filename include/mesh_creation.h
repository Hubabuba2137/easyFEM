#include "imgui.h"
#include "raylib.h"

#include "../include/geometry.h"
#include "mes/structs.h"

void add_new_node(Vector2 worldPos, go::Vertex &mesh);
void pop_node(go::Vertex &mesh);
bool check_if_params_correct(Fem::GlobalData conf);
void print_nodes(std::vector<Fem::Node> &nodes);
