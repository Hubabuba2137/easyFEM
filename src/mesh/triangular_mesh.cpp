#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <fstream>

#include "triangular_mesh.h"

void add_new_node(Vector2 worldPos, go::Vertex &polygon){
    go::Node temp_node(worldPos);

    polygon.vertices.push_back(temp_node);
    polygon.create_edges();
}

void pop_node(go::Vertex &polygon){
    if(!polygon.vertices.empty()){
        polygon.vertices.pop_back();
        polygon.create_edges();
    }
}

namespace msh{

    bool is_node_same(const go::Node &n1, const go::Node &n2){
        float tolerance = 0.1;

        if(std::abs(n1.pos.x - n2.pos.x) < tolerance && std::abs(n1.pos.y - n2.pos.y) < tolerance){
            return true;
        }
        else{
            return false;
        }
    }

    void write_node_ids(std::vector<go::Node> &nodes){
        for(int i=0; i<nodes.size(); i++){
            nodes[i].id = i;
        }
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

    std::vector<go::Node> add_boundary_nodes_on_edge(const go::Segment &seg, int N)
    {
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

    std::vector<go::Node> add_boundary_nodes_on_vertex(const go::Vertex shape, float spacing){
        std::vector<go::Node> nodes;
        float L=0;
        for(auto it:shape.edges){
            L+=it.len();
        }

        for(auto edge:shape.edges){
            int N = edge.len()/spacing;
            std::vector<go::Node> temp_nodes = add_boundary_nodes_on_edge(edge, N);

            for(auto it:temp_nodes){
                nodes.emplace_back(it);
            }
        }

        remove_duplicate_nodes(nodes);
        
        write_node_ids(nodes);

        return nodes;
    }

    go::Triangle super_trian(const std::vector<go::Node> &node_list){
        float max_x, max_y, min_x, min_y;
        
        max_x = node_list[0].pos.x;
        max_y = node_list[0].pos.y;
        min_x = node_list[0].pos.x;
        min_y = node_list[0].pos.y;
        for(auto&it:node_list){
            if(it.pos.x > max_x){
                max_x = it.pos.x;
            }

            if(it.pos.x < min_x){
                min_x = it.pos.x;
            }

            if(it.pos.y > max_y){
                max_y = it.pos.y;
            }

            if(it.pos.y < min_y){
                min_y = it.pos.y;
            }
        }

        float w = max_x - min_x;
        float h = max_y - min_y;

        go::Node p1(min_x - w*0.1f, min_y-h);
        go::Node p2(min_x+w*1.7f, min_y+h*0.5f);
        go::Node p3(min_x-w*0.1f, min_y+h*2.0f);
        
        go::Triangle super_triangle(p1, p2, p3);
        return super_triangle;
    }

    bool inside_circumcircle(const go::Triangle &triangle, const go::Node &point){
        Vector2 A = triangle.points[0].pos;
        Vector2 B = triangle.points[1].pos;
        Vector2 C = triangle.points[2].pos;
        Vector2 D = point.pos;
        
        double adx = A.x - D.x;
        double ady = A.y - D.y;
        double bdx = B.x - D.x;
        double bdy = B.y - D.y;
        double cdx = C.x - D.x;
        double cdy = C.y - D.y;
        
        double bcdet = bdx * cdy - bdy * cdx;
        double cadet = cdx * ady - cdy * adx;
        double abdet = adx * bdy - ady * bdx;
        
        double alift = adx * adx + ady * ady;
        double blift = bdx * bdx + bdy * bdy;
        double clift = cdx * cdx + cdy * cdy;
        
        double det = alift * bcdet + blift * cadet + clift * abdet;
        
        double area_ABC = (B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y);
        
        if (area_ABC == 0.0) {
            return false;
        }
        
        return (det * area_ABC > 0);
    }
    bool same_triangle(const go::Triangle tr1, const go::Triangle tr2) {
        int matches = 0;
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (is_node_same(tr1.points[i], tr2.points[j])) {
                    matches++;
                    break;
                }
            }
        }
        
        return matches == 3;
    }

    inline bool same_edge(const go::Segment &e1, const go::Segment &e2) {
        return (is_node_same(e1.tab[0], e2.tab[0]) && is_node_same(e1.tab[1], e2.tab[1])) ||
            (is_node_same(e1.tab[0], e2.tab[1]) && is_node_same(e1.tab[1], e2.tab[0]));
    }

    bool is_boundary_edge(const go::Segment &edge, const std::vector<go::Triangle> &bad_triangles) {
        int count = 0;
        for (const auto& triangle : bad_triangles) {
            for (const auto& tri_edge : triangle.edges) {
                if (same_edge(edge, tri_edge)) {
                    count++;
                    if (count > 1) return false;
                }
            }
        }
        return count == 1;
    }

    void filter_triangles(std::vector<go::Triangle> &triangles, go::Vertex &polygon){
        std::vector<go::Triangle> temp_triangles;
        temp_triangles.reserve(triangles.size());

        for(auto &tr:triangles){
            float cent_x = (tr.points[0].pos.x + tr.points[1].pos.x + tr.points[2].pos.x)/3;
            float cent_y = (tr.points[0].pos.y + tr.points[1].pos.y + tr.points[2].pos.y)/3;

            go::Node ccentroid(cent_x, cent_y);
            if(polygon.is_node_inside(ccentroid)){
                temp_triangles.push_back(tr);
            }
        }

        triangles = temp_triangles;
    }

    std::vector<go::Triangle> bowyer_watson(std::vector<go::Node>& node_list) {
        if (node_list.empty()) {
            return {};
        }

        std::vector<go::Triangle> triangulation;
        triangulation.reserve(2 * node_list.size());

        go::Triangle super = super_trian(node_list);
        triangulation.push_back(super);

        std::vector<go::Triangle> bad_triangles;
        std::vector<go::Segment> polygon;
        bad_triangles.reserve(triangulation.capacity() / 4);
        polygon.reserve(triangulation.capacity() / 4);

        for (go::Node& node : node_list) {
            bad_triangles.clear();
            
            for (go::Triangle& triangle : triangulation) {
                if (inside_circumcircle(triangle, node)) {
                    bad_triangles.push_back(triangle);
                }
            }
            
            if (bad_triangles.empty()) {
                continue;
            }

            polygon.clear();
            for (const go::Triangle& bad_tr : bad_triangles) {
                for (const go::Segment& edge : bad_tr.edges) {
                    if (is_boundary_edge(edge, bad_triangles)) {
                        polygon.push_back(edge);
                    }
                }
            }

            auto new_end = std::remove_if(triangulation.begin(), triangulation.end(),
                [&bad_triangles](const go::Triangle& tr) {
                    for (const go::Triangle& bad_tr : bad_triangles) {
                        if (same_triangle(tr, bad_tr)) {
                            return true;
                        }
                    }
                    return false;
                });

            triangulation.erase(new_end, triangulation.end());

            triangulation.reserve(triangulation.size() + polygon.size());
            for (const go::Segment& edge : polygon) {
                triangulation.emplace_back(edge.tab[0], edge.tab[1], node);
            }
        }

        std::vector<go::Triangle> final_triangulation;
        final_triangulation.reserve(triangulation.size());
        
        for (const go::Triangle& tr : triangulation) {
            bool shares_super_vertex = false;
            
            for (int i = 0; i < 3 && !shares_super_vertex; i++) {
                for (int j = 0; j < 3; j++) {
                    if (is_node_same(tr.points[i], super.points[j])) {
                        shares_super_vertex = true;
                        break;
                    }
                }
            }
            
            if (!shares_super_vertex) {
                final_triangulation.push_back(tr);
            }
        }

        return final_triangulation;
    }

    float len(go::Node A, go::Node B){
        return sqrt(pow(B.pos.x-A.pos.x,2)+pow(B.pos.y-A.pos.y,2));
    }

    float tr_size(go::Triangle &tr){
        float a = len(tr.points[0], tr.points[1]);
        float b = len(tr.points[1], tr.points[2]);
        float c = len(tr.points[2], tr.points[0]);

        return 0.25*(std::sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)));
    }


    //algorithm taken from:  Owen, Steven. (2000). A Survey of Unstructured Mesh Generation Technology. 7th International Meshing Roundtable. 3.
    //section 2.2.1

    //todo: dodać error handling że jak nie da się dodać punktów bo spacing za duży
    //to coś zrobić
    void triangulate_mesh(go::Vertex polygon, float spacing, std::vector<go::Triangle> &triangles, std::vector<go::Node> &nodes){

        //1. interpolujemy punkty na brzegach
        nodes = add_boundary_nodes_on_vertex(polygon, spacing);
        //std::cout<< int_nodes.size()<<"\n'";
        
        //2. inicjalizacja siatki
        triangles = bowyer_watson(nodes);
        float mean_size = 0.0f; 
        for(go::Triangle tr:triangles){
            mean_size += tr_size(tr);
        }
        mean_size = mean_size/triangles.size();

        constexpr float divider = (4.0f * 1.732f) / 3.0f;

        int max_iter = 10;
        int current_iter = 0;
        while((std::sqrt(divider*mean_size) > spacing) && (current_iter < max_iter)){
            current_iter++;
            //std::cout<<current_iter<<"\n";
            //3. liczymy średnią wielkość trójkątów i średnią odległość między circumcenter
            for(go::Triangle tr:triangles){
                mean_size += tr_size(tr);
            }
            mean_size = mean_size/static_cast<int>(triangles.size());
            
            //std::cout<<std::sqrt(divider*mean_size)<<" <-> " << spacing<<"\n";
            
            //4. dodajemy nowe punkty wewnątrz każdego trójkąta
            //int_nodes.resize(triangles.size()-1);
            for(go::Triangle tr:triangles){
                //for now only the center point of a triangle
                
                //4.2 dodajemy trójkąt jeżeli jego wielkość jest większa od średniej
                if(tr_size(tr) > mean_size){
                    
                    float x_p = (tr.points[0].pos.x+tr.points[1].pos.x+tr.points[2].pos.x)/3;
                    float y_p = (tr.points[0].pos.y+tr.points[1].pos.y+tr.points[2].pos.y)/3;
                    go::Node center(x_p, y_p);
                    
                    if(polygon.is_node_inside(center)){
                        nodes.push_back(center);
                    }
                }
            }

            //5. triangulacja siatki
            triangles = bowyer_watson(nodes);
        }

        filter_triangles(triangles, polygon);

        if(triangles.empty()){
            nodes = polygon.vertices;
            triangles = bowyer_watson(nodes);
        }
    }
}


namespace to_fem{
    float FLOAT_TOLERANCE = 0.01f;

    struct Node_key {
        float x;
        float y;

        bool operator==(const Node_key& other) const {
            return x == other.x && y == other.y;
        }
    };

    struct Node_hashing {
        std::size_t operator()(const Node_key& k) const {
            return std::hash<float>()(k.x) ^ (std::hash<float>()(k.y) << 1);
        }
    };

    
    std::vector<Triangle_ref> convert_to_fem(const std::vector<go::Node> &nodes, const std::vector<go::Triangle> &triangles){
        std::unordered_map<Node_key, int, Node_hashing> node_to_id_map;
        std::vector<Triangle_ref> triangle_refs;


        for (int i = 0; i < nodes.size(); ++i) {
            node_to_id_map[{nodes[i].pos.x, nodes[i].pos.y}] = i;
        }
        

        for (const auto& triangle : triangles) {

            Triangle_ref ref_element;
            for (int i = 0; i < 3; ++i) {
                Node_key key = {triangle.points[i].pos.x, triangle.points[i].pos.y};
                auto it = node_to_id_map.find(key);
                if (it != node_to_id_map.end()) {
                    ref_element.node_ids[i] = it->second;
                } else {
                    ref_element.node_ids[i] = -1; // Or throw an exception
                }
            }
            triangle_refs.push_back(ref_element);
        }

        return triangle_refs;
    }

    void write_to_FEM(const std::vector<go::Node> &nodes, 
                const std::vector<go::Node> &bc_nodes,
                const Fem::GlobalData &conf,
                const std::vector<go::Triangle> mesh)
    {
        std::string file_name = "../Data/fem_data.txt";
        std::ofstream file(file_name, std::ios::out);

        if(!file.is_open()){
            std::cerr << "Error: Could not open file " << file_name << " for writing.\n";
            return;
        }

        std::vector<Triangle_ref> triangle_refs = convert_to_fem(nodes, mesh);

        file << "SimulationTime "      << conf.total_time      << "\n";
        file << "SimulationStepTime "  << conf.time_step       << "\n";
        file << "Conductivity "        << conf.conductivity    << "\n";
        file << "InitialTemp "         << conf.init_temperature<< "\n";
        file << "Density "             << conf.density         << "\n";
        file << "SpecificHeat "        << conf.specific_heat   << "\n";
        file << "Nodes_number "        << nodes.size()         << "\n";
        file << "Elements_number "     << triangle_refs.size()          << "\n";

        file << "*Nodes\n";
        for(size_t i=0;i<nodes.size();i++)
            file << (i+1) << ", " << nodes[i].pos.x << ", " << nodes[i].pos.y << "\n";

        file << "*Elements\n";
        for(size_t i=0;i<triangle_refs.size();i++){
            file << i+1 <<", " << 
            triangle_refs[i].node_ids[0]+1 <<", " << 
            triangle_refs[i].node_ids[1]+1 <<", " << 
            triangle_refs[i].node_ids[2]+1 <<"\n";
        }

        file << "*BC\n";

        std::unordered_map<Node_key, int, Node_hashing> bc_nodes_to_id_map;

        for (int i = 0; i < nodes.size(); ++i) {
            bc_nodes_to_id_map[{nodes[i].pos.x, nodes[i].pos.y}] = i;
        }

        for(auto& b : bc_nodes){
            Node_key key = {b.pos.x, b.pos.y};
            auto it = bc_nodes_to_id_map.find(key);


            if(it != bc_nodes_to_id_map.end()){
                float neumann = 0.f;
                float alfa = 300.f;
                float t_ext = 1200.f;

                int bc_id = it->second;

                file << (bc_id+1) << ", " << neumann << ", " << alfa << ", " << t_ext << "\n";
            }
        }

        std::cout << "File data written\n";
        file.close();
    }

}