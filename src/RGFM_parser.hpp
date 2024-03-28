#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

#include "Greens_PWexpansion.hpp"

using namespace std;
using namespace Eigen;

class RGFM_parser
{
private:
    int periodicity_dim;
    Eigen::Vector3i max_images;
    Eigen::Matrix3d periodic_vecs;//every row is a periodic vector
    int gpw_l_max;
    std::vector<std::string> gpw_species;
    //species must be in the same order they are given to the cosntructor for all instances
    int thread_num;
    int mpi_partition_num;//this is for naively parellizing using mpi. this is th grid num
    std::string mode;
public:
    RGFM_parser(int grid_num);
    void disloc_solve_h();
    void disloc_solve_grid();
};