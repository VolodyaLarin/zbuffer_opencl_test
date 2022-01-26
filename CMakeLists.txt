cmake_minimum_required(VERSION 3.21)
project(zbuffer_test)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCL REQUIRED)


add_executable(zbuffer_test main.cpp ObjModelLoader.cpp
        imgui/imgui.cpp
        imgui/imgui_demo.cpp
        imgui/imgui_draw.cpp
        imgui/imgui_tables.cpp
        imgui/imgui_widgets.cpp
        imgui/backends/imgui_impl_glut.cpp
        imgui/backends/imgui_impl_opengl3.cpp
        )
target_include_directories(zbuffer_test PUBLIC imgui .)


target_link_libraries(zbuffer_test PUBLIC OpenCL freeglut glew32 opengl32)