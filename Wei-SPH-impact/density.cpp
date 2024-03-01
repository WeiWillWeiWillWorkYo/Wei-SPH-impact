#include <cmath>
#include "datatypes.h"
#include "density.h"
#include "kernel.h"
#include "vector.h"

// 计算每个粒子的归一化密度
void normSumDensity(particles_t& parts, int_pairs_t& pairs) {
    double* wi = new double[parts.quant]; // 用于归一化的数组

    // 首先计算空间上权重函数的积分
    for (int k = 0; k < parts.quant; ++k) {
        wi[k] = weight(0.0) * parts.particle[k].m / parts.particle[k].rho;
    }
    for (int k = 0; k < pairs.quant; ++k) {
        int i = pairs.int_pair[k].i;
        int j = pairs.int_pair[k].j;
        wi[i] += pairs.int_pair[k].w * parts.particle[i].m / parts.particle[i].rho;
        wi[j] += pairs.int_pair[k].w * parts.particle[j].m / parts.particle[j].rho;
    }

    // 其次，计算空间上的密度积分
    sumDensity(parts, pairs);

    // 最后，计算归一化密度
    for (int k = 0; k < parts.quant; ++k) {
        parts.particle[k].rho /= wi[k];
    }

    delete[] wi; // 释放分配的内存
}

// 通过简单求和计算密度
void sumDensity(particles_t& parts, int_pairs_t& pairs) {
    for (int k = 0; k < parts.quant; ++k) {
        parts.particle[k].rho = weight(0.0) * parts.particle[k].m;
    }
    for (int k = 0; k < pairs.quant; ++k) {
        int i = pairs.int_pair[k].i;
        int j = pairs.int_pair[k].j;
        parts.particle[i].rho += parts.particle[i].m * pairs.int_pair[k].w;
        parts.particle[j].rho += parts.particle[j].m * pairs.int_pair[k].w;
    }
}

// 计算密度的变化率
void conDensity(particles_t& parts, int_pairs_t& pairs) {
    for (int k = 0; k < parts.quant; ++k) {
        parts.particle[k].drhodt = 0.0;
    }
    for (int k = 0; k < pairs.quant; ++k) {
        int i = pairs.int_pair[k].i;
        int j = pairs.int_pair[k].j;
        double vcc = 0.0;
        for (int beta = 0; beta < DIM; ++beta) {
            double dv = parts.particle[i].v[beta] - parts.particle[j].v[beta];
            vcc += dv * pairs.int_pair[k].dwdx[beta];
        }
        parts.particle[i].drhodt += parts.particle[j].m * vcc;
        parts.particle[j].drhodt += parts.particle[i].m * vcc;
    }
}

// 将当前密度复制到初始密度
void copy2InitDensity(particles_t& parts) {
    for (int k = 0; k < parts.quant; ++k) {
        parts.particle[k].rho_0 = parts.particle[k].rho;
    }
}