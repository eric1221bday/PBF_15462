//
// Created by Eric Fang on 4/10/16.
//

#include "Particle.h"


Particle::Particle(double x, double y, double z, size_t maxNeighbors) {
    this->x = glm::dvec3(x,y,z);
    v = glm::dvec3(0,0,0);
    x_next = glm::dvec3(0,0,0);
    lambda = 0.0f;
    neighbors.reserve(maxNeighbors);
    surface = false;
    density = 0.0f;
}