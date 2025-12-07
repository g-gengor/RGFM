import sys
sys.dont_write_bytecode = True
import py4vasp
import numpy as np
calc = py4vasp.Calculation.from_path('./')
band_data = calc.band
band_dict = band_data.to_dict()
np.savetxt('occupancies.txt',band_dict['occupations'])
print('shape: ',band_dict['occupations'].shape)
