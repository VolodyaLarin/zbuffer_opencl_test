//
// Created by volodya on 24.01.2022.
//

#ifndef ZBUFFER_TEST_MODEL_H
#define ZBUFFER_TEST_MODEL_H


#include <vector>
#include <glm/vec3.hpp>
#include <memory>

struct triangle {
    glm::vec3 t[3];
};
class Model {
public:
    std::vector<triangle> triangles = {};
};


#endif //ZBUFFER_TEST_MODEL_H
