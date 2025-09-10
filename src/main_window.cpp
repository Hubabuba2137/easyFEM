#include "../../_deps/imgui-src/imgui.h"
#include "../../_deps/rlimgui-src/rlImGui.h"
#include "../../_deps/rlimgui-src/rlImGuiColors.h"

#include "raylib.h"
#include "../include/rlgl.h"
#include "../include/raymath.h"

#include <vector>
#include <future>

#include "app_ui/axis.h"
#include "app_ui/grid.h"

#include "mes/mes.h"
#include "mesh/geometry.h"
#include "mesh/mesh.h"
#include "main_window.h"

MainWindow::MainWindow()
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

    std::thread worker_thread; //for fem_solver

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
                polygon.create_edges();
                mesh.clear();
            }

            if (ImGui::CollapsingHeader("Mesh options"))
            {

                ImGui::BeginChild("##child3", ImVec2(0, 150), true);
                ImGui::Text("Node spacing");
                ImGui::SliderFloat("##", &spacing, 0.5f, 5.0f);

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
            ImGui::Text("Initial temperature");
            ImGui::InputFloat("°C T0", &configuration.init_temperature);
            ImGui::Text("Density");
            ImGui::InputFloat("kg/m^3", &configuration.density);
            ImGui::Text("Speific hest");
            ImGui::InputFloat("J/kgK", &configuration.specific_heat);

            if (ImGui::Button("Save")){
                std::cout<<"Save clicked\n";
                if(mesh_created){
                    to_fem::write_to_FEM(glob_nodes, bc_nodes, configuration, mesh);
                }
            }

            if (ImGui::Button("Solve")){
                Fem::Solution solution(DATA_DIR "/fem_data.txt");
                solution.solve(false, true);
                
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
                show_grid_dim(camera);
            }

            polygon.draw((1/camera.zoom));

            for(auto &quad:mesh){
                quad.draw((1/camera.zoom));
            }

            for(auto&it:bc_nodes){
                it.change_color(RED);
                it.draw((1/camera.zoom));
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
                bc_nodes.clear();
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

                mesh = msh::create_quad_mesh(polygon, spacing, glob_nodes);
                bc_nodes = msh::create_bcs(glob_nodes, polygon);

                mesh_created = true;
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
}