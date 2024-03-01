#include <cmath>
#include <vector>
#include "datatypes.h"

// ?Z—R˜°lH?—Ê?’v“I”\—Ê–ÕU
void artificialHeat(particles_t& parts, int_pairs_t& pairs) {
    double g1 = 0.5;
    double g2 = 1.0;

    std::vector<double> divv(parts.quant, 0.0); // g—pvector‰n‰»divv›ó“U[0

    for (int k = 0; k < pairs.quant; ++k) {
        int i = pairs.int_pair[k].i;
        int j = pairs.int_pair[k].j;

        double aux_divv = 0.0;
        for (int alfa = 0; alfa < DIM; ++alfa) {
            double dv = parts.particle[j].v[alfa] - parts.particle[i].v[alfa];
            aux_divv += dv * pairs.int_pair[k].dwdx[alfa];
        }

        divv[i] += parts.particle[j].m * aux_divv / parts.particle[j].rho;
        divv[j] += parts.particle[i].m * aux_divv / parts.particle[i].rho;
    }

    for (int k = 0; k < pairs.quant; ++k) {
        int i = pairs.int_pair[k].i;
        int j = pairs.int_pair[k].j;

        double mrho = (parts.particle[i].rho + parts.particle[j].rho) / 2.0;
        double h_size = (parts.particle[i].h + parts.particle[j].h) / 2.0;

        double modr2 = 0.0;
        double rdwdx = 0.0;
        for (int alfa = 0; alfa < DIM; ++alfa) {
            double dr = parts.particle[i].r[alfa] - parts.particle[j].r[alfa];
            modr2 += dr * dr;
            rdwdx += dr * pairs.int_pair[k].dwdx[alfa];
        }

        double mui = g1 * h_size * parts.particle[i].c + g2 * h_size * h_size * (std::abs(divv[i]) - divv[i]);
        double muj = g1 * h_size * parts.particle[j].c + g2 * h_size * h_size * (std::abs(divv[j]) - divv[j]);
        double muij = (mui + muj);

        double aux = muij * rdwdx / (mrho * (modr2 + 0.01 * h_size * h_size));

        parts.particle[i].dedt += parts.particle[j].m * aux * (parts.particle[i].e - parts.particle[j].e);
        parts.particle[j].dedt += parts.particle[i].m * aux * (parts.particle[j].e - parts.particle[i].e);
    }
}
