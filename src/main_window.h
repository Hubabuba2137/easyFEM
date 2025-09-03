#include "raylib.h"
#include <vector>
#include "mes/mes.h"
#include "mesh/geometry.h"

struct MainWindow{
    const int screenWidth = 1280;
    const int screenHeight = 720;
    float toolbarWidth = 250.0f;

    // camera
    Camera2D camera = { 0 };

    bool showDemoWindow = false;
    bool showXY = true;
    bool creatingMesh = false;
    bool show_grid_bool = true;

    //mesh options
    float spacing = 1.0f;
    bool mesh_created = false;
    
    //main options
    float mouseSensitivity = 8.0f;
    float camera_base_zoom = 80.0f;
    bool resetSensitivity = false;

    //solver options
    bool can_solve=false;
    Fem::GlobalData configuration;

    std::vector<go::Node> glob_nodes;
    std::vector<go::Node> bc_nodes;

    MainWindow();
};