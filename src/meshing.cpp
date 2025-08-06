#include <raylib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <limits>
#include <unordered_set>

#include "../include/geometry.h"
#include "../include/meshing.h"

std::vector<go::Node> add_boundary_nodes_on_edge(go::Segment seg, int N){
    go::Node A = seg.tab[0];
    go::Node B = seg.tab[1];

    auto x_u = [&](float u) { return A.pos.x + (B.pos.x - A.pos.x) * u; };
    auto y_u = [&](float u) { return A.pos.y + (B.pos.y - A.pos.y) * u; };

    std::vector<float> u(N+2);
    u[0] = 0;
    for(float i =1; i<=N; i++){
        u[i] = i/(float)N;
    }

    std::vector<go::Node> nodes;
    for(float it:u){
        go::Node temp(x_u(it),y_u(it));

        nodes.emplace_back(temp);
    }

    return nodes;
}

struct NodeHash {
    std::size_t operator()(const std::pair<float, float>& pos) const {
        return std::hash<float>()(pos.first) ^ (std::hash<float>()(pos.second) << 1);
    }
};

void remove_duplicate_nodes(std::vector<go::Node> &nodes) {
    std::unordered_set<std::pair<float, float>, NodeHash> unique_positions;
    auto it = nodes.begin();
    while (it != nodes.end()) {
        auto pos = std::make_pair(it->pos.x, it->pos.y);
        if (unique_positions.find(pos) != unique_positions.end()) {
            it = nodes.erase(it);
        } else {
            unique_positions.insert(pos);
            ++it;
        }
    }
}

bool is_node_inside_trian(go::Vertex vert, go::Node node) {
    go::Node a = vert.vertices[0];
    go::Node b = vert.vertices[1];
    go::Node c = vert.vertices[2];
    
    Vector2 p = node.pos;
    
    float cross1 = (b.pos.x - a.pos.x) * (p.y - a.pos.y) - (b.pos.y - a.pos.y) * (p.x - a.pos.x);
    float cross2 = (c.pos.x - b.pos.x) * (p.y - b.pos.y) - (c.pos.y - b.pos.y) * (p.x - b.pos.x);
    float cross3 = (a.pos.x - c.pos.x) * (p.y - c.pos.y) - (a.pos.y - c.pos.y) * (p.x - c.pos.x);
    
    bool has_same_sign = (cross1 >= 0 && cross2 >= 0 && cross3 >= 0) || (cross1 <= 0 && cross2 <= 0 && cross3 <= 0);
    return has_same_sign;
}

float cross(const Vector2& a, const Vector2& b, const Vector2& c) {
    //std::cout<<(b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)<<"\n";
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool is_convex(const Vector2& prev, const Vector2& curr, const Vector2& next) {
    return cross(prev, curr, next) < 0; 
}


bool is_ear(std::vector<go::Node>& poly, int i) {
    int prev = (i - 1 + poly.size()) % poly.size();
    int next = (i + 1) % poly.size();

    if (!is_convex(poly[prev].pos, poly[i].pos, poly[next].pos)){
        //std::cout<<"not convex\n";
        return false;
    }

    go::Vertex triangle({ poly[prev], poly[i], poly[next] });

    for (int j = 0; j < poly.size(); ++j) {
        if (j == prev || j == i || j == next){
            continue;
        }
        if (is_node_inside_trian(triangle, poly[j])){
            return false;
        }
    }

    return true;
}

std::vector<go::Vertex> ear_cut_triangulation(go::Vertex& polygon) {
    std::vector<go::Node> verts = polygon.vertices;
    std::vector<go::Vertex> triangles;

    while (verts.size() > 3) {
        bool ear_found = false;
        for (size_t i = 0; i < verts.size(); ++i) {
            if (is_ear(verts, i)) {
                //std::cout<<"Found ear\n";
                int prev = (i - 1 + verts.size()) % verts.size();
                int next = (i + 1) % verts.size();

                std::vector<go::Node> tri = { verts[prev], verts[i], verts[next] };
                triangles.push_back(go::Vertex(tri));
                verts.erase(verts.begin() + i);
                ear_found = true;
                break;
            }
        }

        if (!ear_found) {
            break;
        }
    }

    if (verts.size() == 3) {
        triangles.push_back(go::Vertex(verts));
    }

    //std::cout<<"Triangles size "<<triangles.size()<<"\n";
    return triangles;
}

bool is_node_same(go::Node n1, go::Node n2){
    float tolerance = 0.1;

    if(std::abs(n1.pos.x - n2.pos.x) < tolerance && std::abs(n1.pos.y - n2.pos.y) < tolerance){
        return true;
    }
    else{
        return false;
    }
}

bool have_same_side(go::Vertex v1, go::Vertex v2){

    for(auto&it:v1.edges){
        for(auto&that:v2.edges){
            if(is_node_same(it.tab[0], that.tab[0])||
            is_node_same(it.tab[0], that.tab[1])||
            is_node_same(it.tab[1], that.tab[0])||
            is_node_same(it.tab[1], that.tab[1])){
                return true;
            }
        }
    }

    return false;
}

std::vector<go::Vertex> make_quads(std::vector<go::Vertex> &init_triangles){
    std::vector<go::Vertex> triangles = init_triangles;
    std::vector<go::Vertex> quads;

    int max_iter = 100;
    int iter = 0;

    if(triangles.size()%2==0){
        //łączenie trójkątów w czworokąty
        while(!triangles.empty() && iter < max_iter){
            iter++;
            go::Vertex tr1 = triangles[0];
            go::Vertex tr2 = triangles[1];
    
            if(have_same_side(tr1, tr2)){
                std::vector<go::Node> temp_nodes;
                for(auto&node_1: tr1.vertices){
                    temp_nodes.push_back(node_1);
                }
                for(auto&node_2: tr2.vertices){
                    temp_nodes.push_back(node_2);
                }
    
                remove_duplicate_nodes(temp_nodes);
    
                if(temp_nodes.size() == 4){
                    go::Vertex temp_quad(temp_nodes);
                    quads.push_back(temp_quad);
    
                    triangles.erase(triangles.begin());
                    triangles.erase(triangles.begin());
                }
            }
        }
    }
    else{
        for(auto&triangle:triangles){
            go::Node n1 = triangle.vertices[0];
            go::Node n2 = triangle.vertices[1];
            go::Node n3 = triangle.vertices[2];
    
            go::Node center((n1.pos.x+ n2.pos.x+n3.pos.x)/3,(n1.pos.y+ n2.pos.y+n3.pos.y)/3);
    
            go::Node n4((n1.pos.x+n2.pos.x)/2,(n1.pos.y+n2.pos.y)/2);
            go::Node n5((n2.pos.x+n3.pos.x)/2,(n2.pos.y+n3.pos.y)/2);
            go::Node n6((n3.pos.x+n1.pos.x)/2,(n3.pos.y+n1.pos.y)/2);
    
            std::vector<go::Node> vert1_ns = {n1, n4, center, n6};
            std::vector<go::Node> vert2_ns = {n4, n2, n5, center};
            std::vector<go::Node> vert3_ns = {n5, n3, n6, center};
    
            go::Vertex vert1(vert1_ns);
            go::Vertex vert2(vert2_ns);
            go::Vertex vert3(vert3_ns);
    
            quads.push_back(vert1);
            quads.push_back(vert2);
            quads.push_back(vert3);
        }
    }


    return quads;
}

std::vector<go::Vertex> meshing(go::Vertex polygon, float spacing){
    std::vector<go::Vertex> triangles = ear_cut_triangulation(polygon);
    std::vector<go::Vertex> quads = make_quads(triangles);

    //std::cout<<"Triangle vector size "<<triangles.size()<<"\n";
    //std::cout<<"Quads vector size "<<quads.size()<<"\n";

    float L=0;
    for(auto&edge:polygon.edges){
        L+=edge.len();
    }

    int N = (int)L/spacing;

    std::vector<go::Vertex> result;
    for(auto &quad: quads){
        //creating internal points
        std::vector<go::Node> edge_0_nodes = add_boundary_nodes_on_edge(quad.edges[0], N);
        std::vector<go::Node> edge_2_nodes = add_boundary_nodes_on_edge(quad.edges[2], N);

        std::vector<std::vector<go::Node>> quad_nodes;
        std::vector<go::Node> temp_nodes;
        quad_nodes.reserve(edge_0_nodes.size());

        for(size_t i = 0; i < edge_0_nodes.size(); i++) {
            remove_duplicate_nodes(edge_0_nodes);
            remove_duplicate_nodes(edge_2_nodes);
            
            go::Node A = edge_0_nodes[i];
            go::Node B = edge_2_nodes[edge_0_nodes.size()-1-i];

            go::Segment seg(A, B);
        
            temp_nodes = add_boundary_nodes_on_edge(seg, N);
            //std::cout<<"Temp nodes size "<<temp_nodes.size()<<"\n";
            quad_nodes.push_back(temp_nodes);
        }

        //creating quad elements
        for(size_t row =0; row<quad_nodes.size()-1; row++){
            for( size_t col = 0; col<quad_nodes[row].size()-1; col++){
                std::vector<go::Node> quad_vertices = {
                    quad_nodes[row][col],
                    quad_nodes[row][col+1],
                    quad_nodes[row+1][col+1],
                    quad_nodes[row+1][col]
                };

                go::Vertex temp_quad(quad_vertices);
                result.push_back(temp_quad);
            }
        }
    }

    //std::cout<<"Result size "<<result.size()<<"\n";
    return result;
}