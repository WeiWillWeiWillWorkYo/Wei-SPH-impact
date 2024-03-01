#ifndef DATATYPES_H
#define DATATYPES_H

#include <vector>

constexpr int DIM = 2; // Number of dimensions of the problem

// Global variables can be replaced with a configuration class or structure if needed
extern double norm; // Weight function normalization
extern double dt; // Integration interval

// Particle structure
struct particle_t {
    double a[DIM];                // Acceleration (m/s^2)
    double v[DIM];                // Speed (m/s)
    double r[DIM];                // Position (m)

    double v_prev[DIM];           // Previous speed (for LeapFrog integration)
    double av[DIM];               // Average speed (for correction)

    double m;                     // Mass (kg)
    double p;                     // Pressure (N/m^3)
    double rho;                   // Density (kg/m^3)
    double rho_0;                 // Initial density (kg/m^3)
    double drhodt;                // Density variation (kg/s*m^3)
    double e;                     // Energy (J)
    double dedt;                  // Energy variation (J/s)
    double c;                     // Sound speed (m/s)

    double rho_prev;              // Previous density
    double e_prev;                // Previous energy

    double h;                     // Smoothing size

    double epsilon[DIM][DIM];     // Strain rate tensor
    double tau[DIM][DIM];         // Shear stress tensor
    double tau_dot[DIM][DIM];     // Shear stress rate tensor
    double sigma[DIM][DIM];       // Total stress tensor

    double tau_prev[DIM][DIM];    // Previous shear stress tensor
    double R[DIM][DIM];

    short virt;                   // Virtual particle flag
};

// Particles structure
struct particles_t {
    int quant;                        // Particles quantity
    std::vector<particle_t> particle; // Using vector for automatic memory management
    int alocated;                     // Allocated memory quantity
};

// Interaction pair structure
struct int_pair_t {
    double w;                     // Weight function (0.0 to 1.0)
    double dwdx[DIM];             // Derivative weight function (-1.0 to 1.0)
    int i;                        // First particle of the pair
    int j;                        // Second particle of the pair
};

// Interaction pairs structure
struct int_pairs_t {
    int quant;                        // Pairs quantity
    std::vector<int_pair_t> int_pair; // Using vector for automatic memory management
    int alocated;                     // Allocated memory quantity
};

#endif // DATATYPES_H
