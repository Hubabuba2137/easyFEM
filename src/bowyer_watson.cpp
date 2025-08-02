#include<iostream>
#include <vector>
#include <map>

#include "../include/geometry.h"

/*
float dist(go::Node A, go::Node B){
    return sqrt(pow(A.pos.x - B.pos.x,2)+ pow(A.pos.y - B.pos.y,2));
}

int min_x(std::vector<go::Node> &nodes){
    int min = nodes[0].pos.x;
    for(auto& it: nodes){
        if(it.pos.x < min){
            min = it.pos.x;
        }
    }
    return min;
}

int min_y(std::vector<go::Node> &nodes){
    int min = nodes[0].pos.y;
    for(auto& it: nodes){
        if(it.pos.y < min){
            min = it.pos.y;
        }
    }
    return min;
}

int max_x(std::vector<go::Node> &nodes){
    int max = nodes[0].pos.x;
    for(auto& it: nodes){
        if(it.pos.x >= max){
            max = it.pos.x;
        }
    }
    return max;
}

int max_y(std::vector<go::Node> &nodes){
    int max = nodes[0].pos.y;
    for(auto& it: nodes){
        if(it.pos.y >= max){
            max = it.pos.y;
        }
    }
    return max;
}

go::Triangle super_triangle(std::vector<go::Node> &nodes){

    int x_min = min_x(nodes);
    int y_min = min_y(nodes);

    int x_max = max_x(nodes);
    int y_max = max_y(nodes);


    int x_a, y_a, x_b, y_b, x_c, y_c;

    x_a = x_min - 20;
    y_a = y_min - 20;

    x_b = (-y_min + y_max + x_max )+ 100;
    y_b = y_a;

    x_c = x_a;
    y_c = -x_min + y_max + x_max + 50;

    return go::Triangle(go::Node(x_a, y_a),go::Node(x_b, y_b), go::Node(x_c, y_c)); 
}

bool is_inside_circumcircle(go::Triangle trian, go::Node point){
    //tworzymy okrąg opisany na trójkącie
    float x1 = trian.points[0].pos.x;
    float y1 = trian.points[0].pos.y;

    float x2 = trian.points[1].pos.x;
    float y2 = trian.points[1].pos.y;

    float x3 = trian.points[2].pos.x;
    float y3 = trian.points[2].pos.y;

    float D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

    float Ux = ((x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2)) / D;
    float Uy = ((x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1)) / D;

    float R = sqrt((Ux - x1) * (Ux - x1) + (Uy - y1) * (Uy - y1));

    go::Node center(Ux, Uy);

    //i sprawdzamy czy punkt jest wewnątrz niego
    if(dist(center, point) < R){
        return true;
    }
    else{
        return false;
    }
}

int find_triangle(const std::vector<go::Triangle>& triangulation, const go::Triangle& tri) {
    for (size_t i = 0; i < triangulation.size(); ++i) {
        if (triangulation[i] = tri) {
            return i;
        }
    }
    return -1; 
}

std::vector<go::Triangle> bowyer_watson(std::vector<go::Node> &nodes){
    std::vector<go::Triangle> triangulation;

    //dodajemy super trójkąt    
    go::Triangle super = super_triangle(nodes);
    triangulation.push_back(super);

    for(auto& point:nodes){

        //znajdujemt trójkąty, których okrąg opisany zawiera nowy punkt
        std::vector<go::Triangle> bad_trian;
        for(auto& triangle: triangulation){
            if(is_inside_circumcircle(triangle, point)){
                bad_trian.push_back(triangle);
            }
        }

        //znajdujemy krawędzie graniczne
        std::vector<go::Segment> polygon;
        std::map<go::Segment, int> edgeCount;

        for(auto& bad_tr: bad_trian){
            std::vector<go::Segment> edges = {bad_tr.edges[0], bad_tr.edges[1], bad_tr.edges[2]};

            for(auto& e: edges){
                edgeCount[e]++;
            }
        }

        for (auto& kv : edgeCount) {
            if (kv.second == 1) {
                polygon.push_back(kv.first);
            }
        }

        //usuwamy złe trójkąty
        std::vector<int> indicesToRemove;
		indicesToRemove.reserve(bad_trian.size());

		for (auto& tri : bad_trian) {
			int idx = find_triangle(triangulation, tri);
			if (idx != -1) {
				indicesToRemove.push_back(idx);
			}
		}
		std::sort(indicesToRemove.begin(), indicesToRemove.end(), std::greater<int>());

		for (int idx : indicesToRemove) {
			if (idx >= 0 && idx < (int)triangulation.size()) {
				triangulation.erase(triangulation.begin() + idx);
			}
		}

        //tworzymy nowe trójkąty
        for(auto& edge: polygon){
            triangulation.push_back(go::Triangle(edge.tab[0], edge.tab[1], point));
        }

    }
    std::vector<go::Triangle> result;

    for (auto& tri : triangulation) {
        if (tri.points[0] != super.points[0] && tri.points[0] != super.points[1] && tri.points[0] != super.points[2] &&
            tri.points[1] != super.points[0] && tri.points[1] != super.points[1] && tri.points[1] != super.points[2] &&
            tri.points[2] != super.points[0] && tri.points[2] != super.points[1] && tri.points[2] != super.points[2]) {
            result.push_back(tri);
        }
    }

    return result;
}
*/