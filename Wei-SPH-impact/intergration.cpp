#include "datatypes.h"
#include "integration.h"
#include "vector.h"
#include "particles.h"
#include "kernel.h"
#include "density.h"
#include "stress.h"
#include "pressureandsound.h"
#include "force.h"
#include "viscosity.h"
#include "temperature.h"
#include <stdio.h>
#include <cmath>
#include <numeric>


void singleStep(particles_t& particles, int_pairs_t& int_pairs) {

    directFind(particles, int_pairs);
    conDensity(particles, int_pairs);
    solidsPressureAndSoundSpeed(particles, 1.81, 3630.0, 1.8, 8.0e10);
    totalStressTensor(particles, int_pairs, 8.0e10);
    resetForces(particles);
    artificialViscosity(particles, int_pairs);
    internalForce(particles, int_pairs);
    hUpgrade(particles);
    averageVelocity(particles, int_pairs);

}

void directFind(particles_t& particles, int_pairs_t& int_pairs) {
    int_pairs.int_pair.clear(); // 清空交互对列表

    for (size_t i = 0; i < particles.particle.size(); ++i) {
        for (size_t j = i + 1; j < particles.particle.size(); ++j) {
            double dist[DIM]; // 定义兼容的dist数组
            // 计算两个粒子之间的距离
            for (int d = 0; d < DIM; ++d) {
                dist[d] = particles.particle[i].r[d] - particles.particle[j].r[d];
            }
            double m_dist = std::sqrt(std::inner_product(dist, dist + DIM, dist, 0.0));
            double mh = (particles.particle[i].h + particles.particle[j].h) / 2.0;

            if (m_dist < 2 * mh) { // 如果距离小于 SPH 半径的两倍
                // 添加新的交互对
                int_pairs.int_pair.emplace_back(); // 使用emplace_back而非push_back来避免额外的复制操作
                // 初始化int_pair的成员
                int_pairs.int_pair.back().i = i;
                int_pairs.int_pair.back().j = j;
                // 计算 kernel
                kernel(m_dist, dist, int_pairs.int_pair.back(), mh);
            }
        }
    }
}

// Calculates average speed to fix speed and prevent penetration (Monaghan, 1992)
void averageVelocity(particles_t& parts, int_pairs_t& pairs) {
    double epsilon = 0.1; // A small constant chosen by experience

    // Reset average velocity for all particles
    for (auto& part : parts.particle) {
        std::fill(std::begin(part.av), std::end(part.av), 0.0);
    }

    // Calculate average velocity based on interactions
    for (auto& pair : pairs.int_pair) {
        int i = pair.i;
        int j = pair.j;

        double aux = pair.w / (parts.particle[i].rho + parts.particle[j].rho);

        for (int alfa = 0; alfa < DIM; alfa++) {
            double dv = parts.particle[i].v[alfa] - parts.particle[j].v[alfa];
            parts.particle[i].av[alfa] -= parts.particle[j].m * dv * aux;
            parts.particle[j].av[alfa] += parts.particle[i].m * dv * aux;
        }
    }

    // Apply the small constant epsilon to the average velocity
    epsilon *= 2.0;
    for (auto& part : parts.particle) {
        for (int alfa = 0; alfa < DIM; alfa++) {
            part.av[alfa] *= epsilon;
        }
    }
}

void hUpgrade(particles_t& particles) {
    for (auto& p : particles.particle) {
        if (!p.virt) {
            double aux = -(p.h / (DIM * p.rho)) * p.drhodt;
            p.h += aux * dt;
            if (p.h <= 0.0) {
                p.h -= aux * dt;  // 保持 h 为正值
            }
        }
    }
}
