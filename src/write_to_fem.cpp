#include "mes/structs.h"
#include "../include/geometry.h"
#include <unordered_map>
#include <cmath>
#include <utility>

// Hash function for node positions to help with duplicate detection
struct NodePositionHash {
    std::size_t operator()(const std::pair<float, float>& pos) const {
        // Round to avoid floating point precision issues
        float rounded_x = std::round(pos.first * 1000.0f) / 1000.0f;
        float rounded_y = std::round(pos.second * 1000.0f) / 1000.0f;
        return std::hash<float>()(rounded_x) ^ (std::hash<float>()(rounded_y) << 1);
    }
};

// Check if two nodes are at the same position (within tolerance)
bool nodes_equal(const go::Node& n1, const go::Node& n2, float tolerance = 0.001f) {
    return std::abs(n1.pos.x - n2.pos.x) < tolerance && 
           std::abs(n1.pos.y - n2.pos.y) < tolerance;
}

// Convert go::Node vector to Fem::Node vector and remove duplicates
// Returns pair: first = fem_nodes, second = node_mapping
std::pair<std::vector<Fem::Node>, std::vector<int> > create_fem_nodes(const std::vector<go::Node>& go_nodes) {
    std::vector<Fem::Node> fem_nodes;
    std::vector<int> node_mapping; // Maps original index to new index after deduplication
    std::unordered_map<std::pair<float, float>, int, NodePositionHash> position_to_index;
    
    node_mapping.reserve(go_nodes.size());
    
    for (size_t i = 0; i < go_nodes.size(); ++i) {
        float x = go_nodes[i].pos.x;
        float y = go_nodes[i].pos.y;
        
        // Round positions to avoid floating point precision issues
        float rounded_x = std::round(x * 1000.0f) / 1000.0f;
        float rounded_y = std::round(y * 1000.0f) / 1000.0f;
        std::pair<float, float> pos = std::make_pair(rounded_x, rounded_y);
        
        std::unordered_map<std::pair<float, float>, int, NodePositionHash>::iterator it = position_to_index.find(pos);
        if (it != position_to_index.end()) {
            // Node already exists, map to existing index
            node_mapping.push_back(it->second);
        } else {
            // New unique node
            int new_index = fem_nodes.size();
            fem_nodes.push_back(Fem::Node(rounded_x, rounded_y));
            position_to_index[pos] = new_index;
            node_mapping.push_back(new_index);
        }
    }
    
    return std::make_pair(fem_nodes, node_mapping);
}

// Convert go::Vertex vector to Fem::Element vector using the node mapping
std::vector<Fem::Element> create_fem_elements(const std::vector<go::Vertex>& go_elements, 
                                             const std::vector<int>& node_mapping,
                                             int node_offset = 0) {
    std::vector<Fem::Element> fem_elements;
    fem_elements.reserve(go_elements.size());
    
    int current_node_index = node_offset;
    
    for (size_t elem_id = 0; elem_id < go_elements.size(); ++elem_id) {
        const go::Vertex& go_element = go_elements[elem_id];
        
        // Each go::Vertex should have 4 nodes for quadrilateral elements
        if (go_element.vertices.size() != 4) {
            std::cerr << "Warning: Element " << elem_id << " does not have 4 vertices. Skipping.\n";
            current_node_index += go_element.vertices.size();
            continue;
        }
        
        // Map the local node indices to global deduplicated indices
        int n1 = node_mapping[current_node_index];
        int n2 = node_mapping[current_node_index + 1];
        int n3 = node_mapping[current_node_index + 2];
        int n4 = node_mapping[current_node_index + 3];
        
        fem_elements.push_back(Fem::Element(elem_id, n1, n2, n3, n4));
        current_node_index += 4;
    }
    
    return fem_elements;
}

// Complete conversion function that handles the entire process
std::pair<std::vector<Fem::Node>, std::vector<Fem::Element> >
convert_mesh_to_fem(const std::vector<go::Vertex>& go_elements) {
    // First, collect all nodes from all elements
    std::vector<go::Node> all_go_nodes;
    for (std::vector<go::Vertex>::const_iterator it = go_elements.begin(); it != go_elements.end(); ++it) {
        const go::Vertex& element = *it;
        for (std::vector<go::Node>::const_iterator node_it = element.vertices.begin(); 
             node_it != element.vertices.end(); ++node_it) {
            all_go_nodes.push_back(*node_it);
        }
    }
    
    // Remove duplicates and get mapping
    std::pair<std::vector<Fem::Node>, std::vector<int> > node_result = create_fem_nodes(all_go_nodes);
    std::vector<Fem::Node> fem_nodes = node_result.first;
    std::vector<int> node_mapping = node_result.second;
    
    // Create elements using the mapping
    std::vector<Fem::Element> fem_elements = create_fem_elements(go_elements, node_mapping);
    
    return std::make_pair(fem_nodes, fem_elements);
}

// Alternative version if you already have separate node and element vectors
std::pair<std::vector<Fem::Node>, std::vector<Fem::Element> >
convert_mesh_to_fem(const std::vector<go::Node>& go_nodes, 
                   const std::vector<go::Vertex>& go_elements) {
    // Remove duplicates and get mapping
    std::pair<std::vector<Fem::Node>, std::vector<int> > node_result = create_fem_nodes(go_nodes);
    std::vector<Fem::Node> fem_nodes = node_result.first;
    std::vector<int> node_mapping = node_result.second;
    
    // Create elements using the mapping
    std::vector<Fem::Element> fem_elements = create_fem_elements(go_elements, node_mapping);
    
    return std::make_pair(fem_nodes, fem_elements);
}
