#include "imgui.h"

#include "rlImGui.h"
#include "rlImGuiColors.h"

#include "raylib.h"
#include "include/rlgl.h"
#include "include/raymath.h"

#include "include/geometry.h"
#include "include/meshing.h"
#include "src/axis.cpp"
#include "src/grid.cpp"
#include "src/mesh_creation.cpp"

#include "src/mes/structs.h"

#include "src/write_to_fem.cpp"
#include "src/mes/fem_solve.cpp"

//#include <cstdlib>

/*//first time:
mkdir build
cd build
cmake ..
cmake --build .
./Debug/easyFEM.exe

//building:
cmake --build .
./Debug/easyFEM.exe
*/
 
static const int screenWidth = 1280;
static const int screenHeight = 720;
float toolbarWidth = 250.0f;

// camera
Camera2D camera = { 0 };

bool showDemoWindow = false;
bool showXY = true;
bool creatingMesh = false;
bool show_grid_bool = true;

//mesh options
float spacing = 4;
bool mesh_created = false;
 
//main options
float mouseSensitivity = 8.0f;
float camera_base_zoom = 80.0f;
bool resetSensitivity = false;

Fem::GlobalData configuration;

int main()
{
    InitWindow(screenWidth, screenHeight, "easyFEM");
    SetTargetFPS(60);
    rlImGuiSetup(true);

    camera.offset = {screenWidth/2.0f + 0.5f*toolbarWidth, screenHeight/2.0f};
    camera.target = {0.0f, 0.0f};
    camera.rotation = 0.0f;
    camera.zoom = camera_base_zoom;


    go::Vertex polygon;
    std::vector<go::Vertex> mesh;

    while (!WindowShouldClose())
    {
        // DRAW START
        BeginDrawing();
        ClearBackground(RAYWHITE);

        
        // ------------------ IMGUI LAYOUT ------------------
        rlImGuiBegin();

        ImGui::SetNextWindowPos({0,0}, ImGuiCond_Always);
        ImGui::SetNextWindowSize({(float)screenWidth,(float)screenHeight}, ImGuiCond_Always);
        ImGui::Begin("Root", nullptr,
                     ImGuiWindowFlags_NoDecoration |
                     ImGuiWindowFlags_NoBringToFrontOnFocus);

        // Toolbar
        ImGui::BeginChild("Toolbar", {toolbarWidth,0}, true);
        ImGui::Text("Tools");
        ImGui::Separator();
        
        if (ImGui::Button("Reset Camera")){
            camera.zoom   = camera_base_zoom;
            camera.target = { 0.0f, 0.0f };
        }

        if (ImGui::CollapsingHeader("Mesh Creation"), 1){
            ImGui::BeginChild("##child1", ImVec2(0, 200), true);
            ImGui::Checkbox("Create mesh", &creatingMesh);
            ImGui::Text("Press E to place point");
            ImGui::Text("Press Q to pop points");
            ImGui::Text("Press Enter to create mesh");

            if (ImGui::Button("Reset nodes")){
                polygon.vertices.clear();
            }

            if (ImGui::CollapsingHeader("Mesh options"))
            {

                ImGui::BeginChild("##child3", ImVec2(0, 150), true);
                ImGui::Text("Node spacing");
                ImGui::SliderFloat("##", &spacing, 1, 30);
                
                
                ImGui::EndChild();
            }

            ImGui::EndChild();
        }

        if (ImGui::CollapsingHeader("Solver options"))
        {
            ImGui::BeginChild("##child2", ImVec2(0, 200), true);
            ImGui::Text("Input solver data below");
            
            
            ImGui::Text("Total time");
            ImGui::InputFloat("s", &configuration.total_time);
            ImGui::Text("Time step");
            ImGui::InputFloat("s deltaT", &configuration.time_step);
            ImGui::Text("Conductivity");
            ImGui::InputFloat("W/mk", &configuration.conductivity);
            ImGui::Text("Afla");
            ImGui::InputFloat("Alfa", &configuration.alfa);
            ImGui::Text("External temperature");
            ImGui::InputFloat("°C T_ext", &configuration.tot);
            ImGui::Text("Initial temperature");
            ImGui::InputFloat("°C T0", &configuration.init_temperature);
            ImGui::Text("Density");
            ImGui::InputFloat("kg/m^3", &configuration.density);
            ImGui::Text("Speific hest");
            ImGui::InputFloat("J/kgK", &configuration.specific_heat);
            
            if (ImGui::Button("Solve")){
                std::cout<<"Solving clicked\n";
                if(mesh_created){
                    configuration.node_number = configuration.nodes.size();
                    std::cout<<"configuration nodes number = "<<configuration.node_number<<"\n";
                    configuration.elem_number = configuration.elements.size();
                    std::cout<<"configuration elements number = "<<configuration.elem_number<<"\n"; //-<incorrect

                    if(check_if_params_correct(configuration)){
                        Fem::Solution solution(configuration);
                        solve(solution, false, true);
                    }
                }
            }           

            ImGui::EndChild();
        }

        if (ImGui::CollapsingHeader("Options"))
        {

            ImGui::BeginChild("##child4", ImVec2(0, 200), true);
            ImGui::Checkbox("Show XY", &showXY);
            ImGui::Checkbox("Show grid", &show_grid_bool);

            ImGui::Text("Mouse sensitivity");
            ImGui::SliderFloat("##", &mouseSensitivity, 0.1, 20);
            if(ImGui::Button("Reset sensitivity")){
                mouseSensitivity = 8.0f;
            }
            ImGui::EndChild();
        }

        ImGui::EndChild();
        
        ImGui::SameLine();
        
        ImGui::BeginChild("MainPanel", {0,0}, false);
        ImVec2 panelMin = ImGui::GetCursorScreenPos();
        ImVec2 panelSize = ImGui::GetContentRegionAvail();
        bool panelHover = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByPopup);
        ImGui::EndChild();

        ImGui::End();

        rlImGuiEnd();


        // ------------------ CAMERA INTERACTION ------------------
        // Only pan/zoom when hovering MainPanel
        if (panelHover)
        {
            // Pan
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
            {
                Vector2 delta = GetMouseDelta();
                camera.target.x -= delta.x / camera.zoom;
                camera.target.y -= delta.y / camera.zoom;
            }
            // Zoom
            float wheel = GetMouseWheelMove();
            if (wheel != 0.0f)
            {
                camera.zoom += wheel * mouseSensitivity;
                if (camera.zoom < 10.0f) camera.zoom = 10.0f;
            }
            
        }

        // ------------------ MAIN WINDOW ------------------
        // Creating scissors to cut main view part for the raylib
        int sx = (int)panelMin.x;
        int sy = screenHeight - (int)(panelMin.y + panelSize.y);
        int sw = (int)panelSize.x;
        int sh = (int)panelSize.y;

        rlEnableScissorTest();
        rlScissor(sx, sy, sw, sh);

        BeginMode2D(camera);
            // Whole rendering is being done here

            if(show_grid_bool){
                //show_adaptive_grid(camera, BLUE, RED);
                show_adaptive_grid_dim(camera);
            }

            polygon.draw((1/camera.zoom));

            for(auto &quad:mesh){
                quad.draw((1/camera.zoom));
            }

        EndMode2D();

        if(creatingMesh){
            Vector2 worldPos = GetScreenToWorld2D(GetMousePosition(), camera);
            // Snap to grid
            worldPos.x = round(worldPos.x / 0.1f) * 0.1f;
            worldPos.y = round(worldPos.y / 0.1f) * 0.1f;

            // Convert snapped world position back to screen position for drawing
            Vector2 screenPos = GetWorldToScreen2D(worldPos, camera);
            DrawCircleV(screenPos, 3, RED);

            Vector2 pos_vec({-44, -24 });
            DrawTextEx(GetFontDefault(), TextFormat("[%0.1f, %0.1f]", (float)worldPos.x, (float)worldPos.y), Vector2Add(screenPos, pos_vec), 20, 2, BLACK);

            if(IsKeyPressed(KEY_E)){
                add_new_node(worldPos, polygon);
                if(mesh_created){
                    mesh.clear();
                }
            }

            if(IsKeyPressed(KEY_Q)){
                pop_node(polygon);
                if(mesh_created){
                    mesh.clear();
                }
                
                if(polygon.vertices.size() == 0 && mesh_created){
                    mesh_created = false;
                }
            }

            if(IsKeyPressed(KEY_ENTER)){
                //creating simple mesh  
                mesh = meshing(polygon, spacing);
                
                //creating configuration and mesh object
                
                std::pair<std::vector<Fem::Node>, std::vector<Fem::Element> > pair_of_nodes = convert_mesh_to_fem(mesh);
                
                configuration.nodes = pair_of_nodes.first;
                configuration.elements = pair_of_nodes.second;

                mesh_created = true;
                std::cout<<"creating mesh\n";
            }
        }

        rlDisableScissorTest();

        if(showXY){
            DrawCoordinateAxes(camera, panelMin, panelSize);
        }

        EndDrawing();
        // DRAW END
    }

    rlImGuiShutdown();
    CloseWindow();
    return 0;
}
