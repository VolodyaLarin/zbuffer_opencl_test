#include <iostream>
#include "ObjModelLoader.h"

#include <CL/opencl.hpp>

#include <GL/glut.h>

#include <sstream>

#include "glm/matrix.hpp"
#include "imgui/imgui.h"
#include "imgui/backends/imgui_impl_glut.h"
#include "imgui/backends/imgui_impl_opengl3.h"
#include <glm/gtc/matrix_transform.hpp>

size_t w = 1800, h = 1000;

float speed_k = 1;

cl::Kernel kernel;
cl::Device device;
cl::Context ctx;
cl::Program program;

cl_uchar3 *colorBuffer = nullptr;
cl_float2 *velocityBuffer = nullptr;
cl_float *depthBuffer = nullptr;


cl::Buffer cl_image;
cl::Buffer cl_depth;
cl::Buffer cl_velocity;
cl_uint2 wh = {static_cast<cl_uint>(w), static_cast<cl_uint>(h)};
cl::Buffer cl_wh;


struct TransformationParams {
    glm::vec3 rotate;
    glm::vec3 scale;
    glm::vec3 translate;
};

TransformationParams trans_camera_old = {glm::vec3(0,0,0), glm::vec3(1,1,1), glm::vec3(0,0,0)};
TransformationParams trans_camera = {glm::vec3(0,0,0), glm::vec3(1,1,1), glm::vec3(0,0,0)};

TransformationParams trans_obj_old = {glm::vec3(0,0,0), glm::vec3(1,1,1), glm::vec3(0,0,0)};
TransformationParams trans_obj = {glm::vec3(0,0,0), glm::vec3(1,1,1), glm::vec3(0,100,0)};


std::vector<Model> models = {};


int mode = 1;


void timeFunc(void) {
    glutPostRedisplay();
}

typedef struct __attribute__ ((packed)) {
    cl_uchar4 color;
    cl_float4 points[3];
    cl_float16 mvp;
    cl_float16 mvp_old;
} cl_polygon;



glm::mat4 genetate_srt(glm::mat4 M, glm::vec3 rot, glm::vec3 scale, glm::vec3 transl) {
    M = glm::scale(M, scale);
    M = glm::rotate(M, rot.x, glm::vec3(1.f, 0, 0));
    M = glm::rotate(M, rot.y, glm::vec3(0, 1, 0));
    M = glm::rotate(M, rot.z, glm::vec3(0, 0, 1));
    M = glm::translate(M, transl);
    return M;
}

glm::mat4 generate_mvp(glm::vec3 rot_cam, glm::vec3 scale_cam, glm::vec3 transl_cam, glm::vec3 rot, glm::vec3 scale, glm::vec3 transl){
    glm::mat4 P = glm::mat4 (1);
    P[2][3] = 0.001;

    glm::mat4 V = genetate_srt(glm::mat4(1.f), rot_cam, scale_cam, transl_cam);
    glm::mat4 M = genetate_srt(glm::mat4(1.f), rot, scale, transl);
    glm::mat4 MVP = P * V * M;

    return MVP;
}



void displayMe(void) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGLUT_NewFrame();

    glClear(GL_COLOR_BUFFER_BIT);


    std::vector<cl_polygon> polygons = {};
    for (auto &model: models) {
        for (auto &polygon: model.triangles) {
            float k = 100;
            float c = 0;
            auto color = glm::vec3(1, 0, 0);
            cl_polygon polygon1 = {
                    {cl_uchar(color.x * 255), cl_uchar(color.y * 255), cl_uchar(color.z * 255), 255},
                    {
                     {polygon.t[0].x * k + c, polygon.t[0].y * k + c, polygon.t[0].z * k + c, 1},
                                              {polygon.t[1].x * k + c, polygon.t[1].y * k + c, polygon.t[1].z * k + c, 1},
                                                                       {polygon.t[2].x * k + c, polygon.t[2].y * k + c,polygon.t[2].z * k + c, 1},
                    }

            };


            glm::mat4 MVP = generate_mvp(
                    trans_camera.rotate, trans_camera.scale, trans_camera.translate,
                    trans_obj.rotate, trans_obj.scale, trans_obj.translate
            );
            glm::mat4 MVP_old = generate_mvp(
                    trans_camera_old.rotate, trans_camera_old.scale, trans_camera_old.translate,
                    trans_obj_old.rotate, trans_obj_old.scale, trans_obj_old.translate
            );

            polygon1.mvp = *(cl_float16 *)(&MVP[0][0]);
            polygon1.mvp_old = *(cl_float16 *)(&MVP_old[0][0]);

            polygons.emplace_back(polygon1);
        }
    }

    auto *cl_polygons_buffer = polygons.data();
    auto cl_models_count = polygons.size();

    cl::Buffer cl_models(
            ctx,
            CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR | CL_MEM_COPY_HOST_PTR,
            sizeof(cl_polygon) * cl_models_count,
            cl_polygons_buffer
    );

    kernel.setArg(0, cl_image);
    kernel.setArg(1, cl_depth);
    kernel.setArg(2, cl_wh);
    kernel.setArg(3, cl_models);
    kernel.setArg(4, cl_velocity);


    cl::CommandQueue queue(ctx, device);
    int err_q = 0;
    err_q = queue.enqueueFillBuffer<cl_uchar4>(cl_image, {255, 255, 255, 255}, 0, w * h * sizeof (cl_uchar4));
    err_q = queue.enqueueFillBuffer<cl_float2>(cl_velocity, {0,0}, 0, w * h * sizeof (cl_float2));
    err_q = queue.enqueueFillBuffer<cl_float>(cl_depth, -CL_INFINITY, 0, w * h * sizeof (cl_float));
    err_q = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(cl_models_count));


    err_q = queue.enqueueReadBuffer(cl_image, CL_TRUE, 0, w * h * sizeof(cl_uchar3), colorBuffer);
    err_q = queue.enqueueReadBuffer(cl_depth, CL_TRUE, 0, w * h * sizeof(cl_float), depthBuffer);
    err_q = queue.enqueueReadBuffer(cl_velocity, CL_TRUE, 0, w * h * sizeof(cl_float2), velocityBuffer);

    if (!mode) {
        for (int x = 0; x < w; x++) {
            for (int y = 0; y < h; y++) {
                int inx = x + y * w;
                colorBuffer[inx].x = velocityBuffer[inx].x * speed_k;
                colorBuffer[inx].y = velocityBuffer[inx].y * speed_k;
                colorBuffer[inx].z = 0;
                colorBuffer[inx].w = 255;
            }
        }
    }

    glDrawPixels(w, h, GL_RGBA, GL_UNSIGNED_BYTE, colorBuffer);
    glFlush();

    ImGui::ShowDemoWindow();

    {
        static float f = 0.0f;
        static int counter = 0;

        ImGui::Begin("!");                          // Create a window called "Hello, world!" and append into it.

        ImGui::InputInt("Mode", &mode);
        ImGui::InputFloat3("Camera translate", &trans_camera.translate[0]);
        ImGui::InputFloat3("Camera rotate", &trans_camera.rotate[0]);
//        ImGui::InputFloat3("Camera scale", &trans_camera.scale[0]);

        ImGui::InputFloat3("Object translate", &trans_obj.translate[0]);
        ImGui::InputFloat3("Object rotate", &trans_obj.rotate[0]);
        ImGui::InputFloat3("Object scale", &trans_obj.scale[0]);

        ImGui::InputFloat3("Camera old translate", &trans_camera_old.translate[0]);
        ImGui::InputFloat3("Camera old rotate", &trans_camera_old.rotate[0]);

        ImGui::InputFloat3("Object old translate", &trans_obj_old.translate[0]);
        ImGui::InputFloat3("Object old rotate", &trans_obj_old.rotate[0]);
        ImGui::InputFloat3("Object old scale", &trans_obj_old.scale[0]);
        ImGui::SliderFloat("speed_k", &speed_k, 0, 10);

        ImGui::End();
    }

    ImGui::Render();

    ImGuiIO &io = ImGui::GetIO();
    glViewport(0, 0, (GLsizei) io.DisplaySize.x, (GLsizei) io.DisplaySize.y);

    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    ImGui::EndFrame();

    glutSwapBuffers();
    glutPostRedisplay();

}



static void initGui() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    (void) io;

    ImFontConfig font_config;
    font_config.OversampleH = 1; //or 2 is the same
    font_config.OversampleV = 1;
    font_config.PixelSnapH = 1;

    static const ImWchar ranges[] =
            {
                    0x0020, 0x00FF, // Basic Latin + Latin Supplement
                    0x0400, 0x044F, // Cyrillic
                    0,
            };
    io.Fonts->AddFontFromFileTTF("../fonts/Tahoma.ttf", 16.0f, &font_config, ranges);
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    auto &style = ImGui::GetStyle();
    style.FrameRounding = 10;
    style.WindowBorderSize = 0;
    style.FrameBorderSize = 1;
    ImGui::StyleColorsDark(&style);


    // Setup Platform/Renderer backends
    ImGui_ImplGLUT_Init();


    ImGui_ImplGLUT_InstallFuncs();
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    ObjModelLoader loader;



    auto m = loader.load("../ozyx.obj");
    models.push_back(m);
    models.push_back(m);
    models.push_back(m);
    models.push_back(m);
    models.push_back(m);
    models.push_back(m);
    models.push_back(m);



    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if (platforms.empty()) {
        std::cerr << "No platforms found!" << std::endl;
        exit(1);
    }

    auto platform = platforms.front();

    std::vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);

    if (devices.empty()) std::cerr << "No devices found!" << std::endl;

    device = devices.front();


    std::ifstream kernel_file("../test.cl");
    std::string src(std::istreambuf_iterator<char>(kernel_file), (std::istreambuf_iterator<char>()));

    cl::Program::Sources sources = {src};

    ctx = cl::Context(device);
    program = cl::Program(ctx, sources);

    auto err = program.build();
    if (err != CL_BUILD_SUCCESS) {
        std::cerr << "Error!\nBuild Status: " << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(device)
                  << "\nBuild Log:\t " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
        exit(1);
    }
    std::cerr << "Program was built!" << std::endl;

    kernel = cl::Kernel(program, "render", nullptr);


    colorBuffer = new cl_uchar3[w * h];
    depthBuffer = new cl_float[w * h];
    velocityBuffer = new cl_float2[w * h];


    cl_image = cl::Buffer(
            ctx,
            CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY,
            sizeof(cl_uchar3) * w * h
    );

    cl_depth = cl::Buffer(
            ctx,
            CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY,
            sizeof(cl_float) * w * h
    );

    cl_velocity = cl::Buffer(
            ctx,
            CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY,
            sizeof(cl_float2) * w * h
    );

    cl_wh = cl::Buffer(ctx, CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_USE_HOST_PTR, sizeof(cl_uint2), &wh);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_MULTISAMPLE);

    glutInitWindowSize(w, h);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("SAMPLE TEST");
    glutDisplayFunc(displayMe);
    glutIdleFunc(timeFunc);

    initGui();

    ImGui_ImplOpenGL3_Init();

    glutMainLoop();

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGLUT_Shutdown();
    ImGui::DestroyContext();
    return 0;
}

