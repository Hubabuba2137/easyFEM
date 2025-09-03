#pragma once

#include <raylib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <limits>

//resusing old library made for computational geometry class

namespace go{
    //struct declarations

    struct Node{
        //Node is a single point representem by a small circle
        
        Vector2 pos;
        bool bc = false;
        float radius = 3.0;
        Color node_color = GREEN;

        static int next_id; // Static counter for unique IDs
        int id; // Unique ID for each Node

        Node(float x_in, float y_in);
        Node(Vector2 pos_in);
        Node();
        void draw();
        void draw(float scale);

        void move(Vector2 vec);
        void change_color(Color new_color);
    };
    
    inline bool operator==(const Node& n1, const Node& n2) {
        return (n1.pos.x == n2.pos.x && n1.pos.y == n2.pos.y);
    }
    inline bool operator!=(const Node& n1, const Node& n2) {
        return (n1.pos.x != n2.pos.x && n1.pos.y != n2.pos.y);
    }
    
    struct Segment{
        //segement is made up of two points
        Node tab[2];
    
        Segment(Node node_start, Node node_end);
        Segment();
        void draw();
        void draw(float scale);
    
        void move(Vector2);
    
        float len();
        
        bool operator<(const Segment &other) const {
            if (tab[0].pos.x != other.tab[0].pos.x)
                return tab[0].pos.x < other.tab[0].pos.x;
            if (tab[0].pos.y != other.tab[0].pos.y)
                return tab[0].pos.y < other.tab[0].pos.y;
            if (tab[1].pos.x != other.tab[1].pos.x)
                return tab[1].pos.x < other.tab[1].pos.x;
            return tab[1].pos.y < other.tab[1].pos.y;
        }
    };

    struct Vertex{
        //vertex is just a shape made of multiple points (nodes)
        std::vector<Node> vertices;
        std::vector<Segment> edges;
    
        Vertex(std::vector<Node> nodes);
        Vertex();

        void create_edges();
        void draw();
        void draw(float scale);
        void draw_nodes();
        void add_vertex(Node node);

        void move(int dx, int dy);
        void set_pos(int x, int y);

        bool is_node_inside(Node &point);
        void sort_vertices_by_position();
    };

    //jaki trójkąt jest każdy widzi
    struct Triangle{
        Node points[3];
        Segment edges[3];

        Triangle(Node a, Node b, Node c);
        Triangle();

        void draw();
        void draw(float scale);
    };
}


