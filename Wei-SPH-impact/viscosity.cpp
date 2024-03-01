#include <cmath>
#include "datatypes.h"
#include "vector.h"
#include "pressureandsound.h"

// 计算由于人工粘性导致的加速度
void artificialViscosity(particles_t& parts, int_pairs_t& pairs) {
    double alfa = 2.5;
    double beta = 2.5;

    for (int k = 0; k < pairs.quant; ++k) {
        int i = pairs.int_pair[k].i;
        int j = pairs.int_pair[k].j;

        double vintr = 0.0;  // 内积 v 和 r
        double modr2 = 0.0;  // r 的模平方
        double dv[DIM];      // 速度差
        double dr;           // 位置差

        for (int d = 0; d < DIM; ++d) {
            dv[d] = parts.particle[i].v[d] - parts.particle[j].v[d];
            dr = parts.particle[i].r[d] - parts.particle[j].r[d];
            vintr += dv[d] * dr;
            modr2 += dr * dr;
        }

        if (vintr < 0.0) {
            double h_size = (parts.particle[i].h + parts.particle[j].h) / 2.0;
            double phi = h_size * vintr / (modr2 + 0.01 * h_size * h_size);
            double cmed = (parts.particle[i].c + parts.particle[j].c) / 2.0;
            double rhomed = (parts.particle[i].rho + parts.particle[j].rho) / 2.0;
            double pi = ((-alfa * cmed * phi) + (beta * phi * phi)) / rhomed;

            for (int d = 0; d < DIM; ++d) {
                double aux = -pi * pairs.int_pair[k].dwdx[d];
                parts.particle[i].a[d] += aux * parts.particle[j].m;
                parts.particle[j].a[d] -= aux * parts.particle[i].m;

                // 更新能量的变化率
                parts.particle[i].dedt -= 0.5 * aux * dv[d] * parts.particle[j].m;
                parts.particle[j].dedt -= 0.5 * aux * dv[d] * parts.particle[i].m;
            }
        }
    }
}
