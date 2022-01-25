//
// Created by volodya on 29.08.2021.
//

#ifndef PROJECT_OBJMODELLOADER_H
#define PROJECT_OBJMODELLOADER_H

#include <vector>
#include <fstream>
#include "model.h"
class ObjModelLoader {
public:
    Model load(const std::string& filename) ;
    Model loadStream(std::istream &fin);
};


#endif //PROJECT_OBJMODELLOADER_H
