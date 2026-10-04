The A matrices for the RGFM model are calculated with the notebook
`A_matrix_Calculator.ipynb`. It needs two inputs:

1. The elastic tensor of the bulk material in Voigt form. Its units do not
   matter, but using TPa helps avoid floating-point precision problems.
2. `l_max`: the truncation limit for the degree l in the spherical harmonic
   expansion. It is set to 20 in the notebook. This value caps the `L_MAX` you
   can use when running RGFM: it must be greater than or equal to the `L_MAX`
   in `RGFM_config`.

The A matrices are complex-valued tensors, so their real and imaginary parts are
saved to separate text files (`A_<l>_<m>_r.txt` and `A_<l>_<m>_i.txt`). The main
RGFM code reads and combines them, as long as you keep the file format this tool
writes.
