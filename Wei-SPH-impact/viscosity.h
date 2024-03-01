#pragma once
#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "datatypes.h"

// 计算由于人工粘性导致的加速度（根据Liu的Eq. 4.66）
void artificialViscosity(particles_t& parts, int_pairs_t& pairs);

#endif // VISCOSITY_H
