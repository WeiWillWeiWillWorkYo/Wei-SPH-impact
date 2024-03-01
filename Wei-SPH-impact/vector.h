#ifndef VECTOR_H
#define VECTOR_H

#include "datatypes.h"

// 减去两个向量
void subVector(const double a[DIM], const double b[DIM], double result[DIM]);

// 计算向量的模
double modVector(const double v[DIM]);

#endif // VECTOR_H
