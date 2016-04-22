//
// Created by Eric Fang on 4/10/16.
//
#include <glm/vec3.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/ext.hpp>
#include "Particle.h"
#include <vector>
#include <unordered_map>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "shader.hpp"
#include "controls.hpp"
#ifndef PBF_15462_PARTICLESYSTEM_H
#define PBF_15462_PARTICLESYSTEM_H

struct KeyHash {


    std::size_t operator()(const std::tuple<size_t,size_t,size_t>& k) const
    {
        return ((std::get<0>(k)*73856093)+(std::get<0>(k)*19349663)+(std::get<0>(k)*83492791))%200003;
    }
};

struct KeyEqual {
    bool operator()(const std::tuple<size_t,size_t,size_t>& lhs, const std::tuple<size_t,size_t,size_t>& rhs) const
    {
        return (std::get<0>(lhs)==std::get<0>(rhs) && std::get<1>(lhs)==std::get<1>(rhs) &&
                std::get<2>(lhs)==std::get<2>(rhs));
    }
};

typedef std::unordered_multimap<std::tuple<size_t,size_t,size_t>,Particle *,KeyHash,KeyEqual> hash_map;

class ParticleSystem {
public:
    ParticleSystem();
    ~ParticleSystem();
    void step();
    void add_particle(double x, double y, double z);
    size_t maxNeighbors;
    size_t maxParticles;
    glm::dvec3 gravity;
    glm::dvec3 bounds_min;
    glm::dvec3 bounds_max;
    int iterations;
    double dt;
    double h;
    double rest_density;
    double epsilon;
    double k;
    double delta_q;
    double dist_from_bound;
    double c;
    double poly6_const;
    double spiky_const;
    size_t imax,jmax,kmax;
    std::vector<Particle *> particles;
    std::vector<glm::vec3> positions;
    std::vector<double> scalar_field;
protected:
    void apply_forces();
    void find_neighbors();
    void get_lambda();
    glm::dvec3 get_delta_pos(Particle *i);
    void apply_pressure();
    void collision_check(Particle* i);
    glm::dvec3 get_viscosity(Particle *i);
    double poly6(glm::dvec3 r);
    glm::dvec3 spiky_prime(glm::dvec3 r);
    hash_map neighbor_hash;
    double calc_scalar(size_t i, size_t j, size_t k);
    double calc_cell_density(size_t i, size_t j, size_t k, glm::dvec3 grid_vertex);

};


#endif //PBF_15462_PARTICLESYSTEM_H
