#pragma once
#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "datatypes.h"

// ?算由人工?量?致的能量耗散（Libersky等人提出的高??拉格朗日流体?力学中的Eq. 28）
void artificialHeat(particles_t& particles, int_pairs_t& int_pairs);

#endif // TEMPERATURE_H
