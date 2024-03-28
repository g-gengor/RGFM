#define _TURN_OFF_PLATFORM_STRING
#define _USE_MATH_DEFINES
#include <string>
#include <vector>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include <fstream>
#include <iterator>
#include <iostream>
#include <Eigen/Dense>
#include "boost/math/special_functions/spherical_harmonic.hpp"
#include "boost/math/interpolators/barycentric_rational.hpp"
#include <chrono>

#include <thread>
#include <vector>
#include <functional>
#include <future>
using namespace std;
using namespace Eigen;

class Greens_PWexpansion
{
private:
    int gpw_l_max;
    std::vector<std::string> gpw_species;
    int gpw_total_pw_loop_size;
    Eigen::Matrix<Eigen::dcomplex,Eigen::Dynamic,Eigen::Dynamic> A_matrices;
    std::vector<Eigen::MatrixXd> F_l;
    Eigen::MatrixXd disps; //displacements: species (starting from 1), coords xyz, disps xyz
    //species must be in the same order they are given to the cosntructor for all instances
    int total_atom_num;
    Eigen::VectorXd h;
    Eigen::Vector3d disloc_per_vec; //this is the periodicity vector along the periodic direction of the dislocation
    bool thread_num_set;
    int thread_num;
public:
    Greens_PWexpansion(int l_max, std::vector<std::string> species);
    Greens_PWexpansion(int l_max, std::vector<std::string> species, bool system_solved);
    ~Greens_PWexpansion();
    MatrixXd readMatrix(const char *filename);
    Matrix3d getGreens(int atom_index, double px, double py, double pz);
    VectorXd solve_h();
    MatrixX3d solve4grid(std::string grid_path);
    void construct_eigenstrain_mat();
    void set_dislocation_per_vec(double vx, double v2, double v3);
    Matrix3d getDislocGreens(int atom_index, double px, double py, double pz, int max_image_num);
    VectorXd solve_h_disloc(int max_image_num);
    MatrixX3d solve4grid_disloc(std::string grid_path, int max_image_num);
    MatrixX3d solve4grid_disloc_threaded(std::string grid_path, int max_image_num);
    void G_row_par(int s, double x, double y, double z, int max_image_num, Eigen::Ref<Eigen::MatrixXd> sys_mat); 
    void sys_mat_row_par(int s, int t, double x, double y, double z, int max_image_num, Eigen::Ref<Eigen::MatrixXd> sys_mat); 
    void async_par_for(unsigned start, unsigned end, std::function<void(unsigned i)> fn);
    void thread_par_for(unsigned start, unsigned end, std::function<void(unsigned i)> fn);
    void set_thread_num(int num);
};



