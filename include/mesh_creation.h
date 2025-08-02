#include "imgui.h"
#include "raylib.h"

#include "../include/geometry.h"

void add_new_node(Vector2 worldPos, go::Vertex &mesh);
void pop_node(go::Vertex &mesh);