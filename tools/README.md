# Pre-processing tools for RGFM

These tools compute the input data that RGFM needs. Install the Python
dependencies with `pip install -r tools/requirements.txt`; see
[force_profiles](force_profiles/README.md) for the one dependency that has to be
installed separately (VaspBandUnfolding).

## Anisotropy matrix calculation (`aniso_matrices/`)

Computes the A matrices in Eq. 10 of the
[RGFM paper](https://doi.org/10.1016/j.jmps.2024.105653).

## Force profile calculation (`force_profiles/`)

Computes the force profiles of the different species in the system, following
Section 2.2 of the [RGFM paper](https://doi.org/10.1016/j.jmps.2024.105653).
