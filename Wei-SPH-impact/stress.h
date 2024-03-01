#pragma once
#ifndef STRESS_H
#define STRESS_H

#include "datatypes.h"

// 计算总张力张量（sigma）
void totalStressTensor(particles_t& parts, int_pairs_t& pairs, double mi);

// 检查牵引力是否超过限制，并缩放以避免这种情况（von Mieses, Liu 2003 eq. 8.8）
void plasticYieldModel(particles_t& parts, double J_0);

#endif // STRESS_H

