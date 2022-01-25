#include <iostream>
#include "ObjModelLoader.h"

#include <CL/opencl.hpp>

#include <GL/glut.h>

#include <sstream>


size_t w = 800, h = 600;
cl_uchar3 *colorBuffer = nullptr;

void displayMe(void) {
    for (int i = 0; i < w * h; ++i) {
        colorBuffer[i].w = 255;

    }
//    int x = rand()% (h - 2) * w;
//    for (int i = x; i < w  * 1 + x; ++i) {
//        colorBuffer[i].x += 5;
//        colorBuffer[i].y += 50;
//    }
    glClear(GL_COLOR_BUFFER_BIT);
    if (!colorBuffer) {

        glBegin(GL_POLYGON);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.5, 0.0, 0.0);
        glVertex3f(0.5, 0.5, 0.0);
        glVertex3f(0.0, 0.5, 0.0);
        glEnd();
    } else {
        glDrawPixels(w, h, GL_RGBA, GL_BYTE, colorBuffer);
    }
    glFlush();
}


void timeFunc(void) {
    glutPostRedisplay();
}

typedef struct __attribute__ ((packed)) {
    cl_uchar4 color;
    cl_float3 points[3];
} cl_polygon;


int main(int argc, char *argv[]) {
    srand(time(NULL));

    ObjModelLoader loader;

    std::vector<Model> models = {};


//    auto m = loader.load("../model.obj");
//    models.push_back(m);

    std::stringstream ss(
            "v 1.000000 1.000000 -1.000000\n"
            "v 1.000000 -1.000000 -1.000000\n"
            "v 1.000000 1.000000 1.000000\n"
            "v 1.000000 -1.000000 1.000000\n"
            "v -1.000000 1.000000 -1.000000\n"
            "v -1.000000 -1.000000 -1.000000\n"
            "v -1.000000 1.000000 1.000000\n"
            "v -1.000000 -1.000000 1.000000\n"
            "f 1/1/1 5/2/1 7/3/1 3/4/1\n"
            "f 4/5/2 3/4/2 7/6/2 8/7/2\n"
            "f 8/8/3 7/9/3 5/10/3 6/11/3\n"
            "f 6/12/4 2/13/4 4/5/4 8/14/4\n"
            "f 2/13/5 1/1/5 3/4/5 4/5/5\n"
            "f 6/11/6 5/10/6 1/1/6 2/13/6\n"
    );

    models.push_back(loader.loadStream(ss));


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

    auto device = devices.front();


    std::ifstream kernel_file("../test.cl");
    std::string src(std::istreambuf_iterator<char>(kernel_file), (std::istreambuf_iterator<char>()));

    cl::Program::Sources sources = {src};

    auto ctx = cl::Context(device);
    auto program = cl::Program(ctx, sources);

    auto err = program.build();
    if (err != CL_BUILD_SUCCESS) {
        std::cerr << "Error!\nBuild Status: " << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(device)
                  << "\nBuild Log:\t " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
        exit(1);
    }
    std::cerr << "Program was built!" << std::endl;


    std::vector<cl_polygon> polygons = {};
    for (auto &model: models) {
        for (auto &polygon: model.triangles) {
            float k = 300;
            auto color = glm::vec3(1, 0, 0);
            cl_polygon polygon1 = {
                    {cl_uchar(color.x * 255), cl_uchar(color.y * 255), cl_uchar(color.z * 255), 255},
                    {
                     {polygon.t[0].x * k, polygon.t[0].y * k, polygon.t[0].z * k},
                                              {polygon.t[1].x * k, polygon.t[1].y * k, polygon.t[1].z * k},
                                                                       {polygon.t[2].x * k, polygon.t[2].y * k,
                                                                               polygon.t[2].z * k},
                    }

            };

            std::cout << polygon1.points[0].x << " " << polygon1.points[0].y << " " << polygon1.points[0].z: << std::endl;
            polygons.emplace_back(polygon1);
        }
    }

    auto *cl_models_buffer = new cl_polygon[polygons.size()];
    unsigned cl_models_count = 0;
    for (auto &p: polygons) {
        cl_models_buffer[cl_models_count] = p;
        cl_models_count++;
    }

    cl::Kernel kernel(program, "render", nullptr);

    cl::Buffer cl_image(
            ctx,
            CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY,
            sizeof(cl_uchar3) * w * h
    );

    cl::Buffer cl_depth(
            ctx,
            CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY,
            sizeof(cl_float) * w * h
    );

    cl::Buffer cl_models(
            ctx,
            CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY,
            sizeof(*cl_models_buffer) * cl_models_count,
            cl_models_buffer
    );

    cl_uint2 wh = {static_cast<cl_uint>(w), static_cast<cl_uint>(h)};
    cl::Buffer cl_wh(ctx, CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_USE_HOST_PTR, sizeof(cl_uint2), &wh);


    kernel.setArg(0, cl_image);
    kernel.setArg(1, cl_depth);
    kernel.setArg(2, cl_wh);
    kernel.setArg(3, cl_models);


    cl::CommandQueue queue(ctx, device);
    auto err_q = queue.enqueueFillBuffer<cl_uchar4>(cl_image, {255, 255, 255, 255}, 0, w * h);

    err_q = queue.enqueueFillBuffer<float>(cl_depth, -CL_INFINITY, 0, w * h);

    err_q = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(cl_models_count));


    colorBuffer = new cl_uchar3[w * h];
    auto depthBuffer = new cl_float[w * h];

    err_q = queue.enqueueReadBuffer(cl_image, CL_TRUE, 0, w * h * sizeof(cl_uchar3), colorBuffer);
    err_q = queue.enqueueReadBuffer(cl_depth, CL_TRUE, 0, w * h * sizeof(cl_float), depthBuffer);


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(w, h);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("SAMPLE TEST");
    glutDisplayFunc(displayMe);
    glutIdleFunc(timeFunc);
    glutMainLoop();
    return 0;
}

