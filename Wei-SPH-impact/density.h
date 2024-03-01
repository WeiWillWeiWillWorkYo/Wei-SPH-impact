#pragma once
#ifndef DENSITY_H
#define DENSITY_H

#include "datatypes.h"

// 计算每个粒子的归一化密度（根据 Liu 2003 - Randies 和 Libersky, 1996; Chen et al., 1999a; 2000 的公式4.35）
void normSumDensity(particles_t& parts, int_pairs_t& pairs);

// 通过简单求和计算密度（根据 Liu 2003 的公式4.26）
void sumDensity(particles_t& parts, int_pairs_t& pairs);

// 通过密度的连续性计算密度的变化率（根据 Liu 2003 的公式4.31）
void conDensity(particles_t& parts, int_pairs_t& pairs);

// 将当前密度复制到初始密度
void copy2InitDensity(particles_t& parts);

#endif // DENSITY_H

