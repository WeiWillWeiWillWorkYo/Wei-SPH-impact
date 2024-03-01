#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "datatypes.h"
#include "particles.h"
#include "integration.h"
#include "stress.h"

using namespace std;
using namespace physics;

double norm; // Normalization of the weight function
double dt; // Integration interval
constexpr double M_PI = 3.14159265358979323846;


int main() {
    std::ofstream outputFile("data.dat", std::ios::app);
    if (!outputFile) {
        std::cerr << "Failed to open the output file." << std::endl;
        return EXIT_FAILURE;
    }

    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    particles_t particles = newParticles();
    int_pairs_t int_pairs = newIntPairs();

    double h_size = 0.00038 * 1.2;
    norm = (DIM == 2) ? 15.0 / (7.0 * M_PI * h_size * h_size) : 3.0 / (2.0 * M_PI * h_size * h_size * h_size);

    for (int j = 0; j < 67; ++j) {
        for (int i = 0; i < 20; ++i) {
            addParticle(particles, i * 0.00038 - 0.00361, j * 0.00038 + 0.00019, 0.0, -221.0, 1.13354e-3, 7850, h_size);
        }
    }

    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 76; ++i) {
            addParticle(particles, i * 0.00038 - 0.01425, j * 0.00038 - 0.00171, 0.0, 0.0, 1.13354e-3, 7850, h_size);
            particles.particle.back().virt = 1;
        }
    }

    dt = 0.00000001;
    double t = 0.0;

    singleStep(particles, int_pairs);

    for (auto& p : particles.particle) {
        if (!p.virt) {
            p.rho += dt * p.drhodt / 2.0;
            p.e += dt * p.dedt / 2.0;
            for (int alfa = 0; alfa < DIM; ++alfa) {
                p.v[alfa] += dt * p.a[alfa] / 2.0 + p.av[alfa];
                p.r[alfa] += dt * p.v[alfa];
                for (int beta = 0; beta < DIM; ++beta) {
                    p.tau[alfa][beta] += dt * p.tau_dot[alfa][beta] / 2.0;
                }
            }
        }
    }

    int c = 0;
    for (t = 0.0 + dt; t < 0.0001; t += dt) {
        for (auto& p : particles.particle) {
            if (!p.virt) {
                p.rho_prev = p.rho;
                p.e_prev = p.e;
                for (int alfa = 0; alfa < DIM; ++alfa) {
                    p.v_prev[alfa] = p.v[alfa];
                    p.v[alfa] += dt * p.a[alfa] / 2.0;
                    p.r[alfa] += dt * p.v[alfa];
                    for (int beta = 0; beta < DIM; ++beta) {
                        p.tau_prev[alfa][beta] = p.tau[alfa][beta];
                        p.tau[alfa][beta] += dt * p.tau_dot[alfa][beta] / 2.0;
                    }
                }
            }
        }

        singleStep(particles, int_pairs);

        for (auto& p : particles.particle) {
            if (!p.virt) {
                p.rho = p.rho_prev + dt * p.drhodt;
                p.e = p.e_prev + dt * p.dedt;
                for (int alfa = 0; alfa < DIM; ++alfa) {
                    p.v[alfa] = p.v_prev[alfa] + dt * p.a[alfa] + p.av[alfa];
                    for (int beta = 0; beta < DIM; ++beta) {
                        p.tau[alfa][beta] = p.tau_prev[alfa][beta] + dt * p.tau_dot[alfa][beta];
                    }
                }
            }
        }

        plasticYieldModel(particles, 5.0e8);

        if (c % 100 == 0) {
            for (const auto& p : particles.particle) {
                if (!p.virt) {
                    outputFile << p.r[0] << " " << p.r[1] << "\n";
                }
            }
        }
        ++c;
    }

    outputFile.close();
    return 0;
}
