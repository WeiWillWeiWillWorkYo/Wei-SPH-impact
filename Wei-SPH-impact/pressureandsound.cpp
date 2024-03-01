#include <cmath>
#include "datatypes.h"
#include "pressureandsound.h"

void solidsPressureAndSoundSpeed(particles_t& parts, double gamma, double sound_speed, double slope, double mi) {
    for (int k = 0; k < parts.quant; ++k) {
        double a_0 = parts.particle[k].rho_0 * sound_speed * sound_speed;
        double b_0 = a_0 * (1.0 + 2.0 * (slope - 1.0));
        double c_0 = a_0 * (2.0 * (slope - 1.0) + 3.0 * (slope - 1.0) * (slope - 1.0));

        double eta = (parts.particle[k].rho / parts.particle[k].rho_0) - 1.0;

        double p_H;
        if (eta > 0.0) {
            p_H = a_0 * eta + b_0 * eta * eta + c_0 * eta * eta * eta;
        }
        else {
            p_H = a_0 * eta;
        }

        parts.particle[k].p = (1.0 - gamma * eta / 2.0) * p_H + gamma * parts.particle[k].rho * parts.particle[k].e;

        parts.particle[k].c = std::sqrt(4.0 * mi / (3.0 * parts.particle[k].rho_0));
    }
}
