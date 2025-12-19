cmake_minimum_required(VERSION 3.15)
project(Pliki_rodowe)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Pliki_rodowe
        main.cpp
        matrix.cpp
        matrix.h
        ode_solver.cpp
        ode_solver.h
        opt_alg.cpp
        opt_alg.h
        solution.cpp
        solution.h
        user_funs.cpp
        user_funs.h)
