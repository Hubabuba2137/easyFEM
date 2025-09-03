#include "geometry.h"

#include <algorithm>
#include <string>
#include <cmath>
#include <limits>
#include <iostream>

namespace go {
    int Node::next_id = 0; // Initialize static ID counter

    inline bool operator==(const Triangle& t1, const Triangle& t2) {
        for (int i = 0; i < 3; ++i) {
            bool found = false;
            for (int j = 0; j < 3; ++j) {
                if(t1.points[i] == t2.points[j]) {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        return true;
    }

    //----------- Node -------------
    Node::Node(float x_in, float y_in) : id(next_id++) {
        //pos = (Vector2){x_in, SCREEN_HEIGHT-y_in}; //normalny układ współrzędnych
        Vector2 pos_vec{x_in, y_in};
        pos = pos_vec;
    }

    Node::Node(Vector2 pos_in) : id(next_id++) {
        this->pos = pos_in;
    }

    Node::Node() : id(next_id++) {
        Vector2 pos_vec{0, 0};
        this->pos = pos_vec;
    }

    void Node::draw(){
        DrawCircleV(pos, radius, node_color);
    }

    void Node::draw(float scale){
        DrawCircleV(pos, radius*scale, node_color);
    }

    void Node::move(Vector2 vec){
        this->pos.x += vec.x;
        this->pos.y += vec.y;
    }

    void Node::change_color(Color new_color){
        this->node_color = new_color;
    }

    //----------- Segment -------------
    Segment::Segment(Node node_start, Node node_end){
        tab[0] = node_start;
        tab[1] = node_end;
    }

    Segment::Segment(){}

    void Segment::draw(){
    DrawLineV(tab[0].pos, tab[1].pos, WHITE);
    tab[0].draw();
    tab[1].draw();
}

    void Segment::draw(float scale){
        DrawLineV(tab[0].pos, tab[1].pos, WHITE);
        tab[0].draw(scale);
        tab[1].draw(scale);
    }
    
    void Segment::move(Vector2 vec){
        this->tab[0].pos.x+=vec.x;
        this->tab[0].pos.y+=vec.y;
        this->tab[1].pos.x+=vec.x;
        this->tab[1].pos.y+=vec.y;
    }

    float Segment::len(){
        float dx = tab[1].pos.x - tab[0].pos.x;
        float dy = tab[1].pos.y - tab[0].pos.y;
        return sqrt(dx * dx + dy * dy);
    }
    
    //----------- Triangle -------------
    Triangle::Triangle(Node a, Node b, Node c){
        points[0] = a;
        points[1] = b;
        points[2] = c;

        edges[0] = Segment(a,b);
        edges[1] = Segment(b,c);
        edges[2] = Segment(a,c);
    }

    Triangle::Triangle()
    {
    }
    void Triangle::draw()
    {
        for(Segment it: edges){
            it.draw();
        }
    }
    void Triangle::draw(float scale)
    {
        for(Segment it: edges){
            it.draw(scale);
        }
    }

    //----------- Vertex -------------
    Vertex::Vertex(std::vector<Node> nodes) {
        float centroidX = 0.0f, centroidY = 0.0f;
        for (const auto& node : nodes) {
            centroidX += node.pos.x;
            centroidY += node.pos.y;
        }
        centroidX /= nodes.size();
        centroidY /= nodes.size();

        std::sort(nodes.begin(), nodes.end(), [centroidX, centroidY](const Node& a, const Node& b) {
            float angleA = atan2(a.pos.y - centroidY, a.pos.x - centroidX);
            float angleB = atan2(b.pos.y - centroidY, b.pos.x - centroidX);
            return angleA < angleB;
        });

        vertices = nodes;
        //sort_vertices_by_position();
        create_edges();
    }

    Vertex::Vertex() {
        vertices = {};
        edges = {};
    }

    void Vertex::sort_vertices_by_position() {
        std::sort(vertices.begin(), vertices.end(), [](const Node& a, const Node& b) {
            if (a.pos.x != b.pos.x) {
                return a.pos.x < b.pos.x;
            }
            return a.pos.y < b.pos.y;
        });
        create_edges();
    }
    
    void Vertex::create_edges(){
        std::vector<Segment> temp;
        std::vector<Node> vert = this->vertices;
        
        for(size_t i=0; i<vertices.size();i++){
    
            if(i<vertices.size()-1){
                Segment temp_line(vert[i], vert[i+1]);
                temp.push_back(temp_line);
            }
            else{
                Segment temp_line(vert[i], vert[0]);
                temp.push_back(temp_line);
            }
        }
    
        this->edges = temp;
    }
    void Vertex::draw_nodes(){
        for(size_t i =0; i<vertices.size();i++){
            vertices[i].draw();
        }
    }
    
    void Vertex::draw(){
        for(size_t i =0; i<edges.size();i++){
            edges[i].draw();
        }
    
        for(size_t i =0; i<vertices.size();i++){
            vertices[i].draw();
        }
    }

    void Vertex::draw(float scale){
        for(size_t i =0; i<edges.size();i++){
            edges[i].draw(scale);
        }
    
        for(size_t i =0; i<vertices.size();i++){
            vertices[i].draw(scale);
        }
    }
    
    void Vertex::add_vertex(Node node){
        vertices.push_back(node);
        //sort_vertices_by_position();
        create_edges();
    }

    bool Vertex::is_node_inside(Node &point) {
        bool inside = false;

        std::vector<go::Node> polygon = this->vertices;
        int n = polygon.size();
        for (int i = 0; i < n; i++){
            int j = (i + 1) % n;
            if (((polygon[i].pos.y > point.pos.y) != (polygon[j].pos.y > point.pos.y)) &&
                (point.pos.x < (polygon[j].pos.x - polygon[i].pos.x) *
                (point.pos.y - polygon[i].pos.y) / (polygon[j].pos.y - polygon[i].pos.y) + polygon[i].pos.x)){
                    inside = !inside;
            }
        }
        return inside;
    }
}
