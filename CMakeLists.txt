cmake_minimum_required(VERSION 3.24)
project(easyFEM)

# Workaround for CLang and RayLib's compound literals
# See: https://github.com/raysan5/raylib/issues/1343
set(CMAKE_CXX_STANDARD 11)

# ===== Dependencies =====
include(FetchContent)

include(cmake/raylib.cmake)
include(cmake/imgui.cmake)
include(cmake/rlimgui.cmake)


# ===== Our Project =====
set(SOURCE_FILES
    main.cpp
    src/triangular_mesh.cpp
    src/geometry.cpp
    src/mes/mes.cpp
    src/mes/matrix/matrix.cpp
    src/to_fem.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} rlimgui)

# Checks if OSX and links appropriate frameworks (Only required on MacOS)
if (APPLE)
    target_link_libraries(${PROJECT_NAME} "-framework IOKit")
    target_link_libraries(${PROJECT_NAME} "-framework Cocoa")
    target_link_libraries(${PROJECT_NAME} "-framework OpenGL")
endif()


# Copy the data files
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/Data/
    $<TARGET_FILE_DIR:${PROJECT_NAME}>/Data/)
