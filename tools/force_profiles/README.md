Calculation here follows the formulation given in Section 2.2 of [RGFM](https://doi.org/10.1016/j.jmps.2024.105653)

The force density has 3 main contribution to it:

1. Kinetic component. Calculated using WAVECAR file. See eq 12 in [RGFM](https://doi.org/10.1016/j.jmps.2024.105653)
2. Electrostatic component. Calculated using LOCPOT file
3. Exchange-Correlation component. Calculated uisng AECCAR file

## Steps for Calculation

1. Extract the occupation of bands. VASP decides how many electrons are in each band writes this information in its .h5 output file. This file can be read by using py4vasp package developed by VASP. band_extractor.py script does this and writes the results in a txt file.
2. Kinetic component of the force density calculated as divergence of kinetic component of quantum mechanical stress density. This is calculated using quantum_kinetic_stress.py script. This script uses VASPBandUnfolding package in github (you should be able to find it by googling and install it by following their instructions). This package let's you calculate the electron wavefunctions as a filed in real space based on band and kpoint numbers. The code iterates through k points and band numbers and sums the results. Parallel computation is possible by giving different sets of k points to each job. This is done using the files in k_parts folder. An example slurm script to run these parallel jobs is also given here. 
3. Calculate the force density field and profiles using the force_density notebook. You cna follow along this notebook to see how the calculation works. This notebook first calculates the force density field, then calculates the average profiles. Calculation of average profiles requires DISP file, i.e., atomic positions. Explanation of DISP file is given in the page for main RGFM code. 

In addition to aforementioned scripts, I have added two utility classes: poscar_ops.py and vasp_tools.py. These classes implements functions for reading/writing poscar files and plotting etc. 