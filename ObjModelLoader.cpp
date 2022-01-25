//
// Created by volodya on 29.08.2021.
//
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

#include "ObjModelLoader.h"

Model ObjModelLoader::load(const std::string &filename) {
    std::ifstream fin;
    fin.open(filename);

    if (!fin.is_open()) {
        throw std::exception();
    }

    return loadStream(fin);
}



Model ObjModelLoader::loadStream(std::istream &fin) {
    std::vector<glm::vec3> coords = {};
    std::vector<triangle> polygons = {};

    while (!fin.eof()) {
        std::string lineHeader;
        std::getline(fin, lineHeader);

        std::stringstream str(lineHeader);

        std::string cmd;

        str >> cmd;

        if (cmd == "v") {
            glm::vec3 xyz;
            str >> xyz.x >> xyz.y >> xyz.z;
            coords.emplace_back(xyz);
        } else if (cmd == "f") {
            triangle polygon;
            std::string vertexStr;
            int i = 0;
            while (std::getline(str, vertexStr, ' ')) {
                if (!vertexStr.length()) continue;

                if (i >= 3) continue;

                std::stringstream vertexStream(vertexStr);
                int x;
                vertexStream >> x;

                x--;

                polygon.t[i] = coords[x];
                i++;
            }
            polygons.push_back(polygon);
        }

    }

    Model model;
    model.triangles = std::move(polygons);
    return  std::move(model);
}
