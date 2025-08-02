#include "imgui.h"

#include "raylib.h"
#include "../include/rlgl.h"
#include "../include/raymath.h"
#include "../include/axis.h"

#include "rlImGui.h"
#include "rlImGuiColors.h"

#define LIGHTGRAY_dim CLITERAL(Color){ 200, 200, 200, 100 } 
#define GRAY_dim CLITERAL(Color){ 130, 130, 130, 100 } 

void show_grid(Camera2D camera, float grid_size = 1.0f, Color grid_color = LIGHTGRAY, Color axis_color = GRAY, int screenWidth=1280, int screenHeight=720)
{
    Vector2 top_left = GetScreenToWorld2D({0, 0}, camera);
    Vector2 bottom_right = GetScreenToWorld2D({(float)screenWidth, (float)screenHeight}, camera);
    
    float start_x = floorf(top_left.x / grid_size) * grid_size - grid_size;
    float end_x = ceilf(bottom_right.x / grid_size) * grid_size + grid_size;
    float start_y = floorf(top_left.y / grid_size) * grid_size - grid_size;
    float end_y = ceilf(bottom_right.y / grid_size) * grid_size + grid_size;
    
    for (float x = start_x; x <= end_x; x += grid_size)
    {
        Color line_color = (fabsf(x) < 0.001f) ? axis_color : grid_color; // Oś Y
        DrawLineV({x, start_y}, {x, end_y}, line_color);
    }
    
    for (float y = start_y; y <= end_y; y += grid_size)
    {
        Color line_color = (fabsf(y) < 0.001f) ? axis_color : grid_color; // Oś X
        DrawLineV({start_x, y}, {end_x, y}, line_color);
    }
}

void show_grid_dim(Camera2D camera, float grid_size = 1.0f, int screenWidth=1280, int screenHeight=720)
{
    Vector2 top_left = GetScreenToWorld2D({0, 0}, camera);
    Vector2 bottom_right = GetScreenToWorld2D({(float)screenWidth, (float)screenHeight}, camera);
    
    float start_x = floorf(top_left.x / grid_size) * grid_size - grid_size;
    float end_x = ceilf(bottom_right.x / grid_size) * grid_size + grid_size;
    float start_y = floorf(top_left.y / grid_size) * grid_size - grid_size;
    float end_y = ceilf(bottom_right.y / grid_size) * grid_size + grid_size;
    
    for (float x = start_x; x <= end_x; x += grid_size)
    {
        Color line_color = GRAY_dim; // Oś Y
        DrawLineV({x, start_y}, {x, end_y}, line_color);
    }
    
    for (float y = start_y; y <= end_y; y += grid_size)
    {
        Color line_color = GRAY_dim; // Oś X
        DrawLineV({start_x, y}, {end_x, y}, line_color);
    }
}

void show_adaptive_grid(Camera2D camera, Color grid_color = LIGHTGRAY, Color axis_color = GRAY)
{
    // Automatyczne dostosowanie rozmiaru siatki do poziomu zoom
    float base_grid_size = 1.0f;
    float zoom_factor = camera.zoom / 80.0f; // 80.0f to twój camera_base_zoom
    
    // Wybierz odpowiedni rozmiar siatki
    float grid_size = base_grid_size;
    if (zoom_factor < 0.1f) {
        grid_size = 10.0f; // Duża siatka dla małego zoom
    } else if (zoom_factor < 0.5f) {
        grid_size = 5.0f;  // Średnia siatka
    } else if (zoom_factor < 2.0f) {
        grid_size = 1.0f;  // Normalna siatka
    } else {
        grid_size = 0.5f;  // Drobna siatka
    } 
    show_grid(camera, grid_size, grid_color, axis_color);
}

void show_adaptive_grid_dim(Camera2D camera){
    float base_grid_size = 1.0f;
    float zoom_factor = camera.zoom / 80.0f; // 80.0f to twój camera_base_zoom
    
    // Wybierz odpowiedni rozmiar siatki
    float grid_size = base_grid_size;
    if (zoom_factor < 0.1f) {
        grid_size = 10.0f; // Duża siatka dla małego zoom
    } else if (zoom_factor < 0.5f) {
        grid_size = 5.0f;  // Średnia siatka
    } else if (zoom_factor < 2.0f) {
        grid_size = 1.0f;  // Normalna siatka
    } else {
        grid_size = 0.5f;  // Drobna siatka
    } 
    show_grid_dim(camera, grid_size);
}