#pragma once
#include <vector>
#include "datatypes.h"
#include "particles.h"

//#pragma message("Compiling the required .cpp file")
// Presets for particles
namespace physics {
    particles_t newParticles() {
        particles_t newParticles;
        newParticles.quant = 0;
        newParticles.alocated = ALOC_STEP;
        newParticles.particle.reserve(ALOC_STEP); // Using reserve instead of malloc
        return newParticles;
    }

    // Add new particle
    void addParticle(particles_t& parts, double x, double y, double vx, double vy, double m, double rho_0, double h_size) {
        particle_t newParticle;

        newParticle.r[0] = x;
        newParticle.r[1] = y;
        newParticle.v[0] = vx;
        newParticle.v[1] = vy;
        newParticle.m = m;
        newParticle.rho = rho_0;
        newParticle.rho_0 = rho_0;
        newParticle.e = 0;
        newParticle.dedt = 0;
        newParticle.h = h_size;
        newParticle.virt = 0;

        for (int alfa = 0; alfa < DIM; alfa++) {
            for (int beta = 0; beta < DIM; beta++) {
                newParticle.tau[alfa][beta] = 0;
            }
            newParticle.av[alfa] = 0.0;
        }

        parts.quant++;
        if (parts.quant > parts.alocated) {
            alocMoreParticles(parts); // Allocate more if needed
        }
        parts.particle.push_back(newParticle); // Using push_back instead of manual assignment
    }

    // Allocates more memory for particles
    void alocMoreParticles(particles_t& parts) {
        parts.alocated += ALOC_STEP;
        parts.particle.reserve(parts.alocated); // Using reserve instead of realloc
    }

    // Presets for interaction pairs
    int_pairs_t newIntPairs() {
        int_pairs_t newPairs;
        newPairs.quant = 0;
        newPairs.alocated = ALOC_STEP;
        newPairs.int_pair.reserve(ALOC_STEP);
        return newPairs;
    }

    // Allocates more memory for interaction pairs
    void alocMorePairs(int_pairs_t& pairs) {
        pairs.alocated += ALOC_STEP;
        pairs.int_pair.reserve(pairs.alocated);
    }
}