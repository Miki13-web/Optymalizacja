//Ten plik nie powinien byæ edytowany

#pragma once

#include"matrix.h"
#include"user_funs.h"

matrix* solve_ode(matrix(*)(double, matrix, matrix, matrix), double, double, double, matrix, matrix = 0, matrix = 0); // throw (string);
