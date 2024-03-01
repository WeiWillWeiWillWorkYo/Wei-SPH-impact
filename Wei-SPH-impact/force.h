#pragma once
#ifndef FORCE_H
#define FORCE_H

#include "datatypes.h"

// 重置所有粒子的加速度
void resetForces(particles_t& parts);

// 计算由内部材料力产生的加速度
void internalForce(particles_t& parts, int_pairs_t& pairs);

// 计算由外部力产生的加速度
void externalForce(particles_t& parts);

#endif // FORCE_H
