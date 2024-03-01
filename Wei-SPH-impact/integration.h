#pragma once
#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "datatypes.h"

// Calculate the values for a time step
void singleStep(particles_t& particles, int_pairs_t& int_pairs);

// Find the particles that will make the interaction (that are within the SPH radius) and calculate the weight function for each pair
void directFind(particles_t& particles, int_pairs_t& int_pairs);

// Calculates average speed to fix speed and prevent penetration (Monaghan, 1992)
void averageVelocity(particles_t& parts, int_pairs_t& pairs);

// Update the h size of smooth according to the article of Liu 2005 Eq. 11
void hUpgrade(particles_t& parts);

#endif // INTEGRATION_H




