#include <iostream>
#include <vector>
#include <limits>

#include "../include/geometry.h"

go::Vertex simple_mesh(go::Vertex border, float number_of_rows, int number_of_columns){
    go::Vertex result;

    // Initialize min and max values
    int min_x = std::numeric_limits<int>::max();
    int min_y = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int max_y = std::numeric_limits<int>::min();

    for(const auto& vertex : border.vertices){
        if(vertex.pos.x < min_x){
            min_x = vertex.pos.x;
        }
        if(vertex.pos.y < min_y){
            min_y = vertex.pos.y;
        }
        if(vertex.pos.x > max_x){
            max_x = vertex.pos.x;
        }
        if(vertex.pos.y > max_y){
            max_y = vertex.pos.y;
        }
    }

    float row_separator = static_cast<float>(max_y - min_y) / number_of_rows;
    float column_separator = static_cast<float>(max_x - min_x) / number_of_columns;

    /*
    for(auto &it:border.vertices){
        result.add_vertex(it);
    }*/

    for(size_t i = 0; i < number_of_rows; i++){
        for(size_t j = 0; j < number_of_columns; j++){
            float temp_x = min_x + j * column_separator;
            float temp_y = min_y + i * row_separator;

            go::Node temp(temp_x, temp_y);

            if(border.is_node_inside(temp)){
                result.add_vertex(temp);
            }
        }
    }

    return result;
}