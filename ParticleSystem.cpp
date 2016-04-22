//
// Created by Eric Fang on 4/10/16.
//

#include "ParticleSystem.h"
#include <iostream>

ParticleSystem::ParticleSystem() {
    maxNeighbors = 50;
    maxParticles = 100000;
    gravity = glm::dvec3(0.0f,0.0f,-9.8f);
    bounds_min = glm::dvec3(0.f,0.f,0.f);
    bounds_max = glm::dvec3(40.f,40.f,30.f);
    iterations = 3  ;
    dt = 0.05f;
    h = 1.5f;
    rest_density = 2000;
    epsilon = 0.01f;
    k = 0.1f;
    delta_q = 0.2f*h;
    dist_from_bound = 0.0001f;
    c = 0.01f;
    poly6_const = (315.f/(64.f*glm::pi<double>()*h*h*h*h*h*h*h*h*h));
    spiky_const = (45.f/(glm::pi<double>()*h*h*h*h*h*h));
    imax = size_t(ceil((bounds_max.x-bounds_min.x)/h));
    jmax = size_t(ceil((bounds_max.y-bounds_min.y)/h));
    kmax = size_t(ceil((bounds_max.z-bounds_min.z)/h));
    //particles.reserve(maxParticles);
}

ParticleSystem::~ParticleSystem() {
    for (auto i : particles) {
        delete(i);
    }
}

void ParticleSystem::add_particle(double x, double y, double z) {
    if ((x<=bounds_min.x || x>=bounds_max.x) || (y<=bounds_min.y || y>=bounds_max.y) ||
            (z<=bounds_min.z || z>=bounds_max.z)) {
        std::cout << "particle out of bounds: " << x <<"," << y <<"," << z << std::endl;
        return;}
    if (particles.size()==maxParticles) {return;}
    Particle *particle = new Particle(x,y,z,maxNeighbors);
    particles.push_back(particle);
    positions.push_back(particle->x);

}

double ParticleSystem::poly6(glm::dvec3 r) {
    double norm_coeff = (h*h-glm::length2(r));
    if (norm_coeff<=0) {return double(0);}
    if (r.x==0.0f && r.y==0.0f && r.z==0.0f) {return double(0);}
    return poly6_const*norm_coeff*norm_coeff*norm_coeff;
}

glm::dvec3 ParticleSystem::spiky_prime(glm::dvec3 r) {
    glm::dvec3 r_norm = glm::normalize(r);
    double norm_coeff = (h-glm::l2Norm(r));
    if (norm_coeff<=0) {return glm::dvec3(0.0f);}
    if (r.x==0.0f && r.y==0.0f && r.z==0.0f) {return glm::dvec3(0.0f);}
    return spiky_const*norm_coeff*norm_coeff*r_norm;

}

void ParticleSystem::apply_forces() {
    for (auto i : particles) {
        i->v = i->v + dt*gravity;
        i->x_next = i->x + dt*i->v;
        i->boundary=false;
        //collision_check(i);
    }
}

void ParticleSystem::find_neighbors() {

    neighbor_hash.clear();
    for (auto i : particles) {
        neighbor_hash.emplace(std::make_tuple(floor(i->x_next[0]/h),floor(i->x_next[1]/h),floor(i->x_next[2]/h)),i);
    }
    for (auto i : particles) {
        i->neighbors.clear();
        glm::dvec3 BB_min = i->x_next-glm::dvec3(h,h,h);
        glm::dvec3 BB_max = i->x_next+glm::dvec3(h,h,h)+glm::dvec3(h,h,h);
        for (double x=BB_min.x;x<BB_max.x;x+=h) {
            for (double y=BB_min.y;y<BB_max.y;y+=h) {
                for (double z=BB_min.z;z<BB_max.z;z+=h) {
                    //std::cout << x<<y<<z<<std::endl;
                    auto range = neighbor_hash.equal_range(std::make_tuple(floor(x/h),floor(y/h),floor(z/h)));
                    if (range.first==range.second) { continue;}
                    for(auto it=range.first; it != range.second; ++it) {
                        Particle *j = it->second;
                        if (j != i) {
                            double length = glm::l2Norm(i->x_next,j->x_next);
                            if (length < h) {i->neighbors.push_back(j);}
                        }
                    }
                }
            }
        }
    }
}

double ParticleSystem::calc_cell_density(size_t i, size_t j, size_t k, glm::dvec3 grid_vertex) {
    double scalar=0.0f;
    auto range = neighbor_hash.equal_range(std::make_tuple(i,j,k));
    if (range.first==range.second) { return 0.0f;}
    for(auto it=range.first; it != range.second; ++it) {
        Particle *p = it->second;
        double length = glm::l2Norm(grid_vertex,p->x_next);
        if (length < h) {
            //std::cout << scalar << std::endl;
            scalar+=poly6(grid_vertex-p->x_next);
        }
    }
    return scalar;
}

double ParticleSystem::calc_scalar(size_t i, size_t j, size_t k) {
    double scalar=0.0f;
    glm::dvec3 grid_vertex(i*h,j*h,k*h);
    scalar+=calc_cell_density(i,j,k,grid_vertex);
    scalar+=calc_cell_density(i-1,j,k,grid_vertex);
    scalar+=calc_cell_density(i,j-1,k,grid_vertex);
    scalar+=calc_cell_density(i-1,j-1,k,grid_vertex);
    scalar+=calc_cell_density(i,j,k-1,grid_vertex);
    scalar+=calc_cell_density(i-1,j,k-1,grid_vertex);
    scalar+=calc_cell_density(i,j-1,k-1,grid_vertex);
    scalar+=calc_cell_density(i-1,j-1,k-1,grid_vertex);

    return scalar;
}

void ParticleSystem::get_lambda() {
    for (auto i : particles) {
        //if (i->boundary) {continue;}
        double density_i = 0.0f;
        for (auto j : i->neighbors) {
            density_i+= poly6(i->x_next - j->x_next);
        }
        i->density = density_i;

        double constraint_i = density_i/rest_density - 1.0f;
        double ci_gradient = 0.0f;
        for (auto j : i->neighbors) {
            /*
            if (glm::l2Norm(i->x_next,j->x_next)>h) {
                std::cout << glm::l2Norm(i->x_next,j->x_next) << std::endl;
            }*/
            ci_gradient+=glm::length2(-1.0f/rest_density* spiky_prime(i->x_next - j->x_next));
        }
        glm::dvec3 accum = glm::dvec3(0.0f);
        for (auto j : i->neighbors) {
            accum+= spiky_prime(i->x_next - j->x_next);
            //std::cout <<glm::to_string(spiky_prime(i->x_next - j->x_next))<<","<<glm::to_string(i->x_next)<<","<< glm::to_string(j->x_next)<< std::endl;
        }
        ci_gradient+=glm::length2((1.0f/rest_density)*accum);
        ci_gradient+=epsilon;
        i->lambda=-1.0f * (constraint_i/ci_gradient);
        //std::cout << i->lambda << std::endl;
    }
}

glm::dvec3 ParticleSystem::get_delta_pos(Particle *i) {
    double w_dq = poly6(delta_q*glm::dvec3(1.0f));
    //std::cout << w_dq << std::endl;
    glm::dvec3 delta_pos(0.0f);
    for (auto j : i->neighbors) {
        double kernel_ratio = poly6(i->x_next-j->x_next)/w_dq;
        if (w_dq<glm::epsilon<double>()) {kernel_ratio=0.0f;}
        double scorr = -k*(kernel_ratio*kernel_ratio*kernel_ratio*kernel_ratio*kernel_ratio*kernel_ratio);
        //std::cout << kernel_ratio<< std::endl;
        delta_pos+=(i->lambda+j->lambda+scorr)*spiky_prime(i->x_next-j->x_next);
        //std::cout << j->lambda << std::endl;
    }

    return (1.0f/rest_density)*delta_pos;
}

void ParticleSystem::collision_check(Particle *i) {
    if (i->x_next.x<bounds_min.x) {
        i->x_next.x = bounds_min.x+dist_from_bound;
        i->boundary=true;
        i->v.x=-0.0001*i->v.x;
    }
    if (i->x_next.x>bounds_max.x) {
        i->x_next.x = bounds_max.x-dist_from_bound;
        i->boundary=true;
        i->v.x=-0.0001*i->v.x;
    }
    if (i->x_next.y<bounds_min.y) {
        i->x_next.y = bounds_min.y+dist_from_bound;
        i->boundary=true;
        i->v.y=-0.0001*i->v.y;
    }
    if (i->x_next.y>bounds_max.y) {
        i->x_next.y = bounds_max.y-dist_from_bound;
        i->boundary=true;
        i->v.y=-0.0001*i->v.y;
    }
    if (i->x_next.z<bounds_min.z) {
        i->x_next.z = bounds_min.z+dist_from_bound;
        i->boundary=true;
        i->v.z=-0.0001*i->v.z;
    }
    if (i->x_next.z>bounds_max.z) {
        i->x_next.z = bounds_max.z-dist_from_bound;
        i->boundary=true;
        i->v.z=-0.0001*i->v.z;
    }

}

void ParticleSystem::apply_pressure() {
    for (auto i : particles) {
        //if (i->boundary) {continue;}
        glm::dvec3 dp = get_delta_pos(i);
        i->x_next+=dp;
        collision_check(i);
    }
}

glm::dvec3 ParticleSystem::get_viscosity(Particle *i) {
    glm::dvec3 visc = glm::dvec3(0.0f);
    for (auto j : i->neighbors) {
        visc+=(i->v-j->v)*poly6(i->x-j->x);
    }
    return c*visc;
}

void ParticleSystem::step() {
    apply_forces();

    find_neighbors();
    //get_scalar();
    for (int iter = 0;iter<iterations;iter++) {
        get_lambda();
        apply_pressure();
    }
    for (auto i : particles) {
        i->v = (1.0f/dt)*(i->x_next-i->x);
        i->v+=get_viscosity(i);
        i->x = i->x_next;
    }
}


