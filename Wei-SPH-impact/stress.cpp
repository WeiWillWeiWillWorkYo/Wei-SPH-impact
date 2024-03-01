#include <cmath>
#include "datatypes.h"
#include "vector.h"
#include "stress.h"
#include <stdlib.h>

// 计算总张力张量（sigma）
void totalStressTensor(particles_t& parts, int_pairs_t& pairs, double mi) {
    for (int k = 0; k < parts.quant; ++k) {
        for (int alfa = 0; alfa < DIM; ++alfa) {
            for (int beta = 0; beta < DIM; ++beta) {
                parts.particle[k].epsilon[alfa][beta] = 0.0;
                parts.particle[k].R[alfa][beta] = 0.0;
            }
        }
    }

    for (int k = 0; k < pairs.quant; ++k) {
        int i = pairs.int_pair[k].i;
        int j = pairs.int_pair[k].j;

        double mprhoi = 0.5 * (parts.particle[i].m / parts.particle[i].rho);
        double mprhoj = 0.5 * (parts.particle[j].m / parts.particle[j].rho);

        double dv[DIM];
        subVector(parts.particle[j].v, parts.particle[i].v, dv);

        for (int alfa = 0; alfa < DIM; ++alfa) {
            for (int beta = 0; beta < DIM; ++beta) {
                double aux_epsilon = dv[alfa] * pairs.int_pair[k].dwdx[beta] + dv[beta] * pairs.int_pair[k].dwdx[alfa];
                double aux_R = dv[alfa] * pairs.int_pair[k].dwdx[beta] - dv[beta] * pairs.int_pair[k].dwdx[alfa];

                parts.particle[i].epsilon[alfa][beta] += mprhoj * aux_epsilon;
                parts.particle[j].epsilon[alfa][beta] += mprhoi * aux_epsilon;
                parts.particle[i].R[alfa][beta] += mprhoj * aux_R;
                parts.particle[j].R[alfa][beta] += mprhoi * aux_R;
            }
        }
    }

    // 根据张量更新应力张量和Sigma
    for (int k = 0; k < parts.quant; ++k) {
        double avarage = 0.0;
        for (int alfa = 0; alfa < DIM; ++alfa) {
            avarage += parts.particle[k].epsilon[alfa][alfa];
        }
        avarage /= static_cast<double>(DIM);

        for (int alfa = 0; alfa < DIM; ++alfa) {
            for (int beta = 0; beta < DIM; ++beta) {
                double tau_dot = mi * parts.particle[k].epsilon[alfa][beta];
                if (alfa == beta) {
                    tau_dot -= mi * avarage;
                }
                parts.particle[k].tau_dot[alfa][beta] = tau_dot;

                parts.particle[k].sigma[alfa][beta] = parts.particle[k].tau[alfa][beta];
                if (alfa == beta) {
                    parts.particle[k].sigma[alfa][beta] -= parts.particle[k].p;
                }
            }
        }
    }
}

// 检查牵引力是否超过限制并缩放以避免这种情况
void plasticYieldModel(particles_t& parts, double J_0) {
    for (int k = 0; k < parts.quant; ++k) {
        double J = 0.0;
        for (int alfa = 0; alfa < DIM; ++alfa) {
            for (int beta = 0; beta < DIM; ++beta) {
                J += parts.particle[k].tau[alfa][beta] * parts.particle[k].tau[alfa][beta];
            }
        }
        J = std::sqrt(J);

        if (J > J_0) {
            for (int alfa = 0; alfa < DIM; ++alfa) {
                for (int beta = 0; beta < DIM; ++beta) {
                    parts.particle[k].tau[alfa][beta] *= std::sqrt(J_0 / (3.0 * J * J));
                }
            }
        }
    }
}
