cmake_minimum_required(VERSION 3.21)
project(zbuffer_test)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCL REQUIRED)


add_executable(zbuffer_test main.cpp ObjModelLoader.cpp)
target_link_libraries(zbuffer_test PUBLIC OpenCL freeglut glew32 opengl32)