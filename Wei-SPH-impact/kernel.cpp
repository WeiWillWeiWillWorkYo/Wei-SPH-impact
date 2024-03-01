#include <cmath>
#include "datatypes.h"
#include "kernel.h"

// 在C++中，通常更倾向于使用引用而不是指针
void kernel(double m_dist, const double(&dist)[DIM], int_pair_t& int_pair, double h_size) {
    double q = m_dist / h_size;
    int_pair.w = weightAndGrad(q, dist, m_dist, int_pair.dwdx, h_size); // Weight function and gradient
}

// C++中通常不需要声明d为int，for循环会自动推断
double weightAndGrad(double q, const double(&dist)[DIM], double m_dist, double(&grad)[DIM], double h_size) {
    if (q >= 0.0 && q <= 1.0) {
        for (int d = 0; d < DIM; ++d) {
            grad[d] = norm * (-2.0 + 3.0 / 2.0 * q) / (h_size * h_size) * dist[d];
        }
        return norm * (2.0 / 3.0 - q * q + 0.5 * q * q * q);
    }
    else if (q > 1.0 && q <= 2.0) {
        for (int d = 0; d < DIM; ++d) {
            grad[d] = -norm * 1.0 / 6.0 * 3.0 * std::pow(2.0 - q, 2) / h_size * (dist[d] / m_dist);
        }
        return norm * 1.0 / 6.0 * std::pow(2.0 - q, 3);
    }
    else {
        for (int d = 0; d < DIM; ++d) {
            grad[d] = 0.0;
        }
        return 0.0;
    }
}

// 使用C++的重载函数特性，我们可以定义两个名称相同但参数不同的函数
double weight(double q) {
    if (q >= 0.0 && q <= 1.0) {
        return norm * (2.0 / 3.0 - q * q + 0.5 * q * q * q);
    }
    else if (q > 1.0 && q <= 2.0) {
        return norm * 1.0 / 6.0 * std::pow(2.0 - q, 3);
    }
    else {
        return 0.0;
    }
}