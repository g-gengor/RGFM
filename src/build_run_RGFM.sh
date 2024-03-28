#!/bin/bash
#g++ -std=c++14 -c Greens_PWexpansion.cpp disloc_h_solver_LP.cpp -I /home/gorkem/anaconda3/envs/GFM_c/include/eigen3 -I/home/gorkem/anaconda3/envs/GFM_c/include -O2
#g++ -std=c++14 -o LP_h Greens_PWexpansion.o disloc_h_solver_LP.o -lpthread -O2
g++ -std=c++14 -c RGFM_parser.cpp run_RGFM.cpp Greens_PWexpansion.cpp -I /home/gorkem/anaconda3/envs/GFM_c/include/eigen3 -I/home/gorkem/anaconda3/envs/GFM_c/include -O2
g++ -std=c++14 -o run_RGFM RGFM_parser.o run_RGFM.o Greens_PWexpansion.o -lpthread -O2
mv run_RGFM ../bin/run_RGFM
