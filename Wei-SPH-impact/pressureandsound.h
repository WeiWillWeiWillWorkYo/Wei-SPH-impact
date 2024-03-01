#pragma once
#ifndef PRESSUREANDSOUND_H
#define PRESSUREANDSOUND_H

#include "datatypes.h"

// 计算粒子的压力和声速（使用 Liu 的Mie-Gruneisen状态方程为固体计算的公式8.9和8.17）
void solidsPressureAndSoundSpeed(particles_t& parts, double gamma, double sound_speed, double slope, double mi);

#endif // PRESSUREANDSOUND_H
