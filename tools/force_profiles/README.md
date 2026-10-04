The calculation here follows the formulation in Section 2.2 of the
[RGFM paper](https://doi.org/10.1016/j.jmps.2024.105653).

The force density has three main contributions:

1. Kinetic component, calculated from the WAVECAR file (see Eq. 12 of the
   [RGFM paper](https://doi.org/10.1016/j.jmps.2024.105653)).
2. Electrostatic component, calculated from the LOCPOT file.
3. Exchange-correlation component, calculated from the AECCAR files.

## Dependencies

`pip install -r ../requirements.txt` installs the Python packages imported here
(numpy, scipy, pandas, matplotlib, py4vasp). The kinetic-stress script also
needs `vaspwfc` and `aewfc` from
[VaspBandUnfolding](https://github.com/QijingZheng/VaspBandUnfolding), which is
not on PyPI: install it by following the instructions in its repository and
make sure it is on your `PYTHONPATH`.

## Steps for the calculation

1. **Extract the band occupations.** VASP decides how many electrons are in
   each band and writes this information to its `.h5` output file, which can be
   read with VASP's `py4vasp` package. `band_extractor.py` does this and writes
   the occupations to `occupancies.txt`.
2. **Kinetic component.** The kinetic component of the force density is the
   divergence of the kinetic part of the quantum-mechanical stress density,
   calculated by `quantum_kinetic_stress.py`. The script uses VaspBandUnfolding
   (see above) to build the electron wavefunctions as real-space fields from
   their band and k-point indices, then iterates over k-points and bands and
   sums the results. The calculation can be parallelized by giving each job a
   different set of k-points through the files in `k_parts/`
   (`python quantum_kinetic_stress.py <part id>` reads `k_parts/k_part_<part id>.txt`).
   `sum_kin_forces.py` then adds the partial results. Input paths
   (for example `../vasp_NV/WAVECAR`) and the thread count are set at the top
   of the scripts; adjust them to your directory layout. An example batch
   (SLURM) script for these parallel jobs is not included in this repository.
3. **Force density field and profiles.** Use the `force_density.ipynb`
   notebook and follow along to see how the calculation works. It first
   calculates the force density field, then the averaged profiles. The averaged
   profiles need the `DISP` file (atomic positions); its format is described in
   the [main README](../../README.md#disp).

`poscar_ops.py` and `vasp_tools.py` are utility classes for reading and
writing POSCAR files, plotting, and similar tasks.
