#include "vector.h"
#include <cmath> // 引入cmath库用于sqrt函数
#include "datatypes.h"

// 减去两个向量
void subVector(const double a[DIM], const double b[DIM], double result[DIM]) {
    for (int i = 0; i < DIM; ++i) {
        result[i] = a[i] - b[i];
    }
}

// 计算向量的模
double modVector(const double v[DIM]) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1]); // 假设是2维向量
    //return sqrt(v[0]*v[0] + v[1]* v[1];
}
