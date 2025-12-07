A matrices for RGFM model can be calculated by using this notebook. Two inputs are necessary:

1. Elastic tensor for bulk material in Voigt form. The units of this is not important; however, using TPa helps avoiding floting point number precision problems.
2. l_max: this is the truncation limit for parameter l in spherical harmonic expansion. It is given as 20 here. This number limits what you can choose as l_max whe running RGFM code. This number has to be greater than equal to l_max you choose while running RGFM.

The calculated A matrices are complex valued tensors. Hence, we save their imaginary and real parts separetely in text files. Main RGFM code is developed to read and combine them as long as you follow format this tool uses. 