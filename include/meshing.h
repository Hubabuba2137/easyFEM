#pragma once

#include <raylib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <limits>

#include "geometry.h"

std::vector<go::Node> add_boundary_nodes_on_edge(go::Segment seg, int N);
void remove_duplicate_nodes(std::vector<go::Node> &nodes);

bool is_node_inside_trian(go::Vertex vert, go::Node node);
float cross(const Vector2& a, const Vector2& b, const Vector2& c);
bool is_convex(const Vector2& prev, const Vector2& curr, const Vector2& next);
bool is_ear(const std::vector<go::Node>& poly, int i);
std::vector<go::Vertex> ear_cut_triangulation(const go::Vertex& polygon);


bool is_node_same(go::Node n1, go::Node n2);
bool have_same_side(go::Vertex v1, go::Vertex v2);
std::vector<go::Vertex> make_quads(std::vector<go::Vertex> &init_triangles);
std::vector<go::Vertex> meshing(go::Vertex polygon, float spacing);