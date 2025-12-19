#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix, matrix);
matrix ff1R(matrix, matrix, matrix);
matrix df1(double, matrix, matrix, matrix);
matrix dff1R(double, matrix, matrix, matrix);
matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);
matrix df3(double t, matrix Y, matrix ud1, matrix omega_m);
matrix ff3T_inside(matrix x, matrix a, matrix ud2);
matrix ff3T_outside(matrix x, matrix a, matrix ud2);
matrix ff3R(matrix x, matrix ud1, matrix /*unused*/);

matrix ff4T(matrix, matrix = NAN, matrix = NAN);
matrix gf4T(matrix, matrix = NAN, matrix = NAN);
matrix hf4T(matrix, matrix = NAN, matrix = NAN);

matrix ff4R(matrix, matrix = NAN, matrix = NAN);
matrix gf4R(matrix, matrix = NAN, matrix = NAN);