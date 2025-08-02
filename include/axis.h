#ifndef AXIS_H
#define AXIS_H

#include "imgui.h"
#include "raylib.h"


void DrawCoordinateAxes(Camera2D camera, ImVec2 panelMin, ImVec2 panelSize);
void DrawAxisX(float worldLeft, float worldRight, ImVec2 panelMin, ImVec2 panelSize, float axisHeight);
void DrawAxisY(float worldTop, float worldBottom, ImVec2 panelMin, ImVec2 panelSize, float axisWidth, float axisHeight, Camera2D camera);
float CalculateAxisStep(float range);


#endif
