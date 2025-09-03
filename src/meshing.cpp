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

std::vector<go::Node> add_boundary_nodes_on_vertex(go::Vertex shape, float spacing){
    std::vector<go::Node> nodes;
    float L=0;
    for(auto it:shape.edges){
        L+=it.len();
    }

    int N_total = (int)(L/spacing);
    int N = N_total/shape.edges.size();

    for(auto edge:shape.edges){

        std::vector<go::Node> temp_nodes = add_boundary_nodes_on_edge(edge, (int)N);

        for(auto it:temp_nodes){
            nodes.emplace_back(it);
        }
    }

    remove_duplicate_nodes(nodes);
    
    write_node_ids(nodes);

    return nodes;
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
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool is_convex(const Vector2& prev, const Vector2& curr, const Vector2& next) {
    return cross(prev, curr, next) > 0; 
}


bool is_ear(const std::vector<go::Node>& poly, int i) {
    int prev = (i - 1 + poly.size()) % poly.size();
    int next = (i + 1) % poly.size();

    if (!is_convex(poly[prev].pos, poly[i].pos, poly[next].pos)){
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

    if(verts.size()<3){
        std::cerr<<"Polygon too small to triangulate\n";
    }

    while (verts.size() > 3) {
        bool ear_found = false;

        for (size_t i = 0; i < verts.size(); ++i) {
            if (is_ear(verts, i)) {
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
            if((is_node_same(it.tab[0], that.tab[0]) && is_node_same(it.tab[1], that.tab[1])) ||
                (is_node_same(it.tab[0], that.tab[1]) && is_node_same(it.tab[1], that.tab[0]))) {
                return true;
            }
        }
    }

    return false;
}

bool have_same_side(go::Triangle v1, go::Triangle v2){

    for(auto&it:v1.edges){
        for(auto&that:v2.edges){
            if((is_node_same(it.tab[0], that.tab[0]) && is_node_same(it.tab[1], that.tab[1])) ||
                (is_node_same(it.tab[0], that.tab[1]) && is_node_same(it.tab[1], that.tab[0]))) {
                return true;
            }
        }
    }

    return false;
}

namespace {
    using Pair = std::pair<int,int>;

    struct Candidate {
        int a;
        int b;
        int priority; // lower = better (sum of degrees)
    };

    // --- Angle/shape helpers -------------------------------------------------
    static inline float dot2D(const go::Node& a, const go::Node& b){ return a.pos.x*b.pos.x + a.pos.y*b.pos.y; }
    static inline go::Node sub2D(const go::Node& a, const go::Node& b){ return go::Node(a.pos.x-b.pos.x, a.pos.y-b.pos.y); }
    static inline float len2(const go::Node& v){ return v.pos.x*v.pos.x + v.pos.y*v.pos.y; }

    static inline float angle_between_deg(const go::Node& u, const go::Node& v){
        const float uu = std::max(1e-12f, len2(u));
        const float vv = std::max(1e-12f, len2(v));
        float c = (u.pos.x*v.pos.x + u.pos.y*v.pos.y) / std::sqrt(uu*vv);
        if(c> 1.f) c = 1.f; if(c<-1.f) c = -1.f;
        return std::acos(c) * 180.0f / 3.14159265358979323846f;
    }

    static std::vector<go::Node> order_ccw(const std::vector<go::Node>& nodes){
        go::Node c(0.f,0.f);
        for(const auto& n: nodes){ c.pos.x+=n.pos.x; c.pos.y+=n.pos.y; }
        c.pos.x/=nodes.size(); c.pos.y/=nodes.size();
        std::vector<std::pair<float,go::Node>> tmp; tmp.reserve(nodes.size());
        for(const auto& n: nodes){
            float ang = std::atan2(n.pos.y - c.pos.y, n.pos.x - c.pos.x);
            tmp.emplace_back(ang, n);
        }
        std::sort(tmp.begin(), tmp.end(), [](auto& A, auto& B){ return A.first < B.first; });
        std::vector<go::Node> out; out.reserve(nodes.size());
        for(auto& p: tmp) out.push_back(p.second);
        return out;
    }

    // lenient shape check: reject quads with nearly collinear triplets (triangle-ish)
    // Defaults chosen to be permissive: only filter extreme cases
    static bool quad_shape_ok(const std::vector<go::Node>& nodes,
                              float min_angle_deg = 6.0f,
                              float max_angle_deg = 174.0f){
        if(nodes.size()!=4) return false;
        auto P = order_ccw(nodes);
        for(int i=0;i<4;i++){
            const go::Node& pi = P[i];
            const go::Node& pp = P[(i+3)&3]; // previous
            const go::Node& pn = P[(i+1)&3]; // next
            go::Node a = sub2D(pp, pi);
            go::Node b = sub2D(pn, pi);
            float ang = angle_between_deg(a,b);
            if(!(ang > min_angle_deg && ang < max_angle_deg)) return false;
        }
        return true;
    }

    // --- Pairing utilities ---------------------------------------------------
    // Find the shared edge of two triangles; also return the third vertex of each triangle
    static bool shared_edge_and_tips(const go::Triangle& t1, const go::Triangle& t2,
                                     go::Segment& shared,
                                     go::Node& tip1, go::Node& tip2){
        // locate shared edge
        for(const auto& e1: t1.edges){
            for(const auto& e2: t2.edges){
                if(same_edge(e1, e2)){
                    shared = e1; // arbitrary orientation

                    // tip1 is vertex of t1 not on shared
                    for(const auto& n: t1.points){
                        if(!is_node_same(n, e1.tab[0]) && !is_node_same(n, e1.tab[1])){
                            tip1 = n; break;
                        }
                    }
                    // tip2 is vertex of t2 not on shared
                    for(const auto& n: t2.points){
                        if(!is_node_same(n, e2.tab[0]) && !is_node_same(n, e2.tab[1])){
                            tip2 = n; break;
                        }
                    }
                    return true;
                }
            }
        }
        return false;
    }

    static std::vector<std::vector<int>> build_adjacency(const std::vector<go::Triangle>& tris){
        const int n = (int)tris.size();
        std::vector<std::vector<int>> adj(n);
        for(int i=0;i<n;i++){
            for(int j=i+1;j<n;j++){
                if(have_same_side(tris[i], tris[j])){
                    adj[i].push_back(j);
                    adj[j].push_back(i);
                }
            }
        }
        return adj;
    }

    static std::vector<Candidate> build_candidates(const std::vector<std::vector<int>>& adj){
        const int n = (int)adj.size();
        std::vector<int> deg(n);
        for(int i=0;i<n;i++) deg[i] = (int)adj[i].size();

        std::vector<Candidate> cands;
        for(int i=0;i<n;i++){
            for(int j: adj[i]) if(j>i){
                cands.push_back({i,j, deg[i]+deg[j]});
            }
        }
        std::sort(cands.begin(), cands.end(), [](const Candidate& A, const Candidate& B){
            if(A.priority != B.priority) return A.priority < B.priority; // prefer leaves
            // tie-breaker: lower indices stable
            if(A.a != B.a) return A.a < B.a;
            return A.b < B.b;
        });
        return cands;
    }

    static void triangle_into_three_quads(const go::Triangle& tri, std::vector<go::Vertex>& out){
        const go::Node n1 = tri.points[0];
        const go::Node n2 = tri.points[1];
        const go::Node n3 = tri.points[2];

        go::Node center((n1.pos.x + n2.pos.x + n3.pos.x)/3.0f,
                        (n1.pos.y + n2.pos.y + n3.pos.y)/3.0f);

        go::Node n4((n1.pos.x+n2.pos.x)/2.0f, (n1.pos.y+n2.pos.y)/2.0f);
        go::Node n5((n2.pos.x+n3.pos.x)/2.0f, (n2.pos.y+n3.pos.y)/2.0f);
        go::Node n6((n3.pos.x+n1.pos.x)/2.0f, (n3.pos.y+n1.pos.y)/2.0f);

        out.emplace_back(go::Vertex(std::vector<go::Node>{n1, n4, center, n6}));
        out.emplace_back(go::Vertex(std::vector<go::Node>{n4, n2, n5, center}));
        out.emplace_back(go::Vertex(std::vector<go::Node>{n5, n3, n6, center}));
    }

    static void vertex_triangle_into_three_quads(const go::Vertex& tri, std::vector<go::Vertex>& out){
        const go::Node n1 = tri.vertices[0];
        const go::Node n2 = tri.vertices[1];
        const go::Node n3 = tri.vertices[2];

        go::Node center((n1.pos.x + n2.pos.x + n3.pos.x)/3.0f,
                        (n1.pos.y + n2.pos.y + n3.pos.y)/3.0f);

        go::Node n4((n1.pos.x+n2.pos.x)/2.0f, (n1.pos.y+n2.pos.y)/2.0f);
        go::Node n5((n2.pos.x+n3.pos.x)/2.0f, (n2.pos.y+n3.pos.y)/2.0f);
        go::Node n6((n3.pos.x+n1.pos.x)/2.0f, (n3.pos.y+n1.pos.y)/2.0f);

        out.emplace_back(go::Vertex(std::vector<go::Node>{n1, n4, center, n6}));
        out.emplace_back(go::Vertex(std::vector<go::Node>{n4, n2, n5, center}));
        out.emplace_back(go::Vertex(std::vector<go::Node>{n5, n3, n6, center}));
    }
}


std::vector<go::Vertex> make_quads(std::vector<go::Triangle>& input_tris){
    const int n = (int)input_tris.size();
    if(n == 0) return {};

    // 1) Build dual-graph adjacency and candidate list
    auto adj = build_adjacency(input_tris);
    auto cands = build_candidates(adj);

    // 2) Greedy maximal matching over triangle graph
    std::vector<char> used(n, 0);
    std::vector<go::Vertex> quads; quads.reserve(n/2 + n); // generous

    for(const auto& c: cands){
        if(used[c.a] || used[c.b]) continue;
        go::Segment shared; go::Node ta, tb;
        if(!shared_edge_and_tips(input_tris[c.a], input_tris[c.b], shared, ta, tb)) continue;

        std::vector<go::Node> nodes = {ta, shared.tab[0], tb, shared.tab[1]};
        remove_duplicate_nodes(nodes);
        if(nodes.size() != 4) continue; // degenerate, skip

        // Vertex ctor sorts nodes around centroid and constructs edges for us.
        // Only accept if quad is not triangle-like
        if(!quad_shape_ok(nodes)) continue;
        quads.emplace_back(go::Vertex(order_ccw(nodes)));
        used[c.a] = used[c.b] = 1;
    }

    // 3) Split any remaining single triangles into three quads each
    for(int i=0;i<n;i++) if(!used[i]){
        triangle_into_three_quads(input_tris[i], quads);
    }

    return quads;
}

std::vector<go::Vertex> make_quads(std::vector<go::Vertex>& input_tris){
    const int n = (int)input_tris.size();
    if(n == 0) return {};

    // Build adjacency by checking shared sides (Vertex overload exists)
    std::vector<std::vector<int>> adj(n);
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            if(have_same_side(input_tris[i], input_tris[j])){
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }

    // Candidates (sum of degrees priority)
    std::vector<int> deg(n); for(int i=0;i<n;i++) deg[i] = (int)adj[i].size();
    std::vector<Candidate> cands;
    for(int i=0;i<n;i++) for(int j:adj[i]) if(j>i) cands.push_back({i,j,deg[i]+deg[j]});
    std::sort(cands.begin(), cands.end(), [](const Candidate&a,const Candidate&b){
        if(a.priority!=b.priority) return a.priority<b.priority;
        if(a.a!=b.a) return a.a<b.a; return a.b<b.b;
    });

    std::vector<char> used(n,0);
    std::vector<go::Vertex> quads; quads.reserve(n/2 + n);

    for(const auto& c: cands){
        if(used[c.a] || used[c.b]) continue;
        // Build quad nodes: union of the 2 triangles' vertices
        std::vector<go::Node> nodes; nodes.reserve(6);
        for(const auto& v: input_tris[c.a].vertices) nodes.push_back(v);
        for(const auto& v: input_tris[c.b].vertices) nodes.push_back(v);
        remove_duplicate_nodes(nodes);
        if(nodes.size()!=4) continue; // skip degenerate
        if(!quad_shape_ok(nodes)) continue; // avoid triangle-like quads
        quads.emplace_back(go::Vertex(order_ccw(nodes)));
        used[c.a]=used[c.b]=1;
    }

    for(int i=0;i<n;i++) if(!used[i]){
        vertex_triangle_into_three_quads(input_tris[i], quads);
    }

    return quads;
}


std::vector<go::Node> creating_nodes(go::Vertex polygon, float spacing){
    //std::vector<go::Vertex> triangles = ear_cut_triangulation(polygon);
    
    //std::vector<go::Vertex> quads = make_quads(triangles);

    std::vector<go::Vertex> quads;
    quads.push_back(polygon);
    
    float L=0;
    for(auto&edge:polygon.edges){
        L+=edge.len();
    }

    int N = (int)L/spacing;

    std::vector<go::Node> result;
    for(auto &quad: quads){
        //creating internal points
        std::vector<go::Node> edge_0_nodes = add_boundary_nodes_on_edge(quad.edges[0], N);
        std::vector<go::Node> edge_2_nodes = add_boundary_nodes_on_edge(quad.edges[2], N);

        std::vector<go::Node> temp_nodes;
        
        for(size_t i = 0; i < edge_0_nodes.size(); i++) {
            remove_duplicate_nodes(edge_0_nodes);
            remove_duplicate_nodes(edge_2_nodes);
            
            go::Node A = edge_0_nodes[i];
            go::Node B = edge_2_nodes[edge_0_nodes.size()-1-i];

            go::Segment seg(A, B);
        
            temp_nodes = add_boundary_nodes_on_edge(seg, N);
            
            for(auto&node:temp_nodes){
                result.push_back(node);
            }
        }
    }
    return result;
}

go::Triangle super_trian(std::vector<go::Node> &node_list){
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

bool inside_circumcircle(go::Triangle &triangle, go::Node &point){
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
bool same_triangle(go::Triangle tr1, go::Triangle tr2) {
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

inline bool same_edge(const go::Segment& e1, const go::Segment& e2) {
    return (is_node_same(e1.tab[0], e2.tab[0]) && is_node_same(e1.tab[1], e2.tab[1])) ||
           (is_node_same(e1.tab[0], e2.tab[1]) && is_node_same(e1.tab[1], e2.tab[0]));
}

bool is_boundary_edge(const go::Segment& edge, const std::vector<go::Triangle>& bad_triangles) {
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
    polygon.reserve(32);

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

std::vector<go::Node> create_bcs(go::Vertex polygon, float spacing){
    float L=0;
    for(auto&edge:polygon.edges){
        L+=edge.len();
    }

    int N = (int)L/spacing;

    std::vector<go::Node> bcs;

    for(go::Segment &edge:polygon.edges){
        std::vector<go::Node> temp_nodes = add_boundary_nodes_on_edge(edge, N);

        for(go::Node &node:temp_nodes){
            bcs.push_back(node);
        }
    }

    return bcs;
}

void create_mesh(go::Vertex polygon, float spacing, std::vector<go::Vertex> &quad_mesh, std::vector<go::Node> &glob_nodes){
    //interpoalting nodes
    //std::cout<<"creating mesh...\n";

    glob_nodes = creating_nodes(polygon, spacing);
    //std::cout<<"\tnodes size: "<<glob_nodes.size()<<"\n";
    
    //creating triangle mesh with delonay triangulation
    std::vector<go::Triangle> triangles = bowyer_watson(glob_nodes);
    //std::cout<<"\ttriangles size: "<<triangles.size()<<"\n";

    //filtering triangles outside polygon
    filter_triangles(triangles, polygon);

    //connecting triangles into quads
    quad_mesh = make_quads(triangles);
    //std::cout<<"\tquad mesh size: "<<quad_mesh.size()<<"\n";
}