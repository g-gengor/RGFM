Central compilation of RGFM code.

src/ folder contains the source code and building bash files. 

RGFM calculations are configured via RGFM_config file. There are several tags that needs to be given in this file

PERIODIC_DIM: dimension of periodicity (Default 0). For dislocations, it is 1. For point defects, it is 0.

TO DO: imlement point defect codes

L_MAX: the maximum degree of spherical harmonics in spherical harmonic expansion (Default 8). For 4H-SiC, 8 is good enough.

MAX_IMAGES: maximum number of images (both positive and negative) along the periodicity vector (Default 0). This is relevant for dislcoations, and irrelevant for point defects. 

SPECIES: the number of species and species names (Mandatory). It should be given in the file in the following fashion: SPECIES 3 C Si N

PERIODIC_VECS: these are the vectors along which the system is periodic. For a dislocation, for example, PERIODIC_VECS 3.062 0.0 0.0

TO DO: Implement higher periodicities in RGFM code. Extra periodic vectors will be added to the above example string.

THREAD_NUM: number of threads to run parallel fors on in RGFM code (Default 1).

PROBLEM: type of problem to be solved (Mandatory). PROBLEM H solves the H scaling vectors for given configuration. PROBLEM GRID solves the grid problem.

TO DO: Implement point defect solutions wrapper here. Currently these only work for dislocation cases.

Apart from the RGFM_config file, you can pass a grid number through CLI as follows.
./run_RGFM 127
Above command will look for GRID_127 file in GRIDS folder and write the dislplacement field into US/U_127 when configured for a grid solution. For H calculations, this parameter does not do anything.

Apart form the RGFM_config file, for RGFM, you have to provide DISP file that gives the atomic displacment data, real and imaginary parts of anisotropy tensors in A/ folder, and force profile data in F/F_species0, F/F_species1, ...

DISP file is organized as a matrix of Nx8 where N is the number of atoms. This file is organized as follows.

Column 0: atom id numbers. Currently, RGFM does not use this information.

Column 1: species numbers. This enumerations starts from 1 and follows the order given by SPECIES tag in RGFM_config file. For example, for SPECIES 3 C Si N, C is 1, Si is 2, and N is 3.

Column 2, 3, 4: displaced positions of the atoms. These positions must be given in the same coordinate system that anisotropy tensors are calculated. 

Column 5, 6, 7: displacement vectors. These displacements must be given in the same coordinate system that anisotropy tensors are calculated. 
