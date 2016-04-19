//
// Created by Eric Fang on 4/10/16.
//

#ifndef PBF_15462_PARTICLE_H
#define PBF_15462_PARTICLE_H


#include <glm/glm.hpp>
#include <vector>

class Particle {
public:
    Particle(double x, double y, double z, size_t maxNeighbors);
    glm::dvec3 x, v, x_next;
    double lambda, density;
    std::vector<Particle *> neighbors;
    bool boundary, surface;
};


#endif //PBF_15462_PARTICLE_H
