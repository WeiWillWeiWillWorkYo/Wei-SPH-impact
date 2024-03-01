#include "datatypes.h"
#include "force.h"
#include "vector.h"

// 重置所有粒子的加速度
void resetForces(particles_t& parts) {
    for (int k = 0; k < parts.quant; ++k) {
        for (int alfa = 0; alfa < DIM; ++alfa) {
            parts.particle[k].a[alfa] = 0.0; // 初始化加速度为0
        }
        parts.particle[k].dedt = 0.0; // 初始化能量变化率为0
    }
}

// 计算由内部材料力产生的加速度
void internalForce(particles_t& parts, int_pairs_t& pairs) {
    for (int k = 0; k < pairs.quant; ++k) {
        int i = pairs.int_pair[k].i;
        int j = pairs.int_pair[k].j;
        double rhoirhoj = parts.particle[i].rho * parts.particle[j].rho;
        double dv[DIM];
        subVector(parts.particle[j].v, parts.particle[i].v, dv); // 计算速度差

        double aux_e = 0.0;
        for (int alfa = 0; alfa < DIM; ++alfa) {
            aux_e += (parts.particle[i].p + parts.particle[j].p) * dv[alfa] * pairs.int_pair[k].dwdx[alfa] / rhoirhoj; // 计算能量变化率的辅助变量

            double aux_a = 0.0;
            for (int beta = 0; beta < DIM; ++beta) {
                aux_a += (parts.particle[i].sigma[alfa][beta] + parts.particle[j].sigma[alfa][beta]) * pairs.int_pair[k].dwdx[beta] / rhoirhoj; // 计算加速度的辅助变量
            }

            parts.particle[i].a[alfa] += parts.particle[j].m * aux_a;
            parts.particle[j].a[alfa] -= parts.particle[i].m * aux_a;
        }

        parts.particle[i].dedt += 0.5 * parts.particle[j].m * aux_e;
        parts.particle[j].dedt += 0.5 * parts.particle[i].m * aux_e;
    }

    // 更新能量变化率
    for (int k = 0; k < parts.quant; ++k) {
        for (int alfa = 0; alfa < DIM; ++alfa) {
            for (int beta = 0; beta < DIM; ++beta) {
                parts.particle[k].dedt += parts.particle[k].tau[alfa][beta] * parts.particle[k].epsilon[alfa][beta] / parts.particle[k].rho;
            }
        }
    }
}

// 计算由外部力产生的加速度
void externalForce(particles_t& parts) {
    // 这里添加计算外部力的代码，如果您的模拟中使用了外部力
}
