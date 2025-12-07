import sys
sys.dont_write_bytecode = True
import numpy as np
import gc
from poscar_ops import poscar_ops
from vasp_tools import vasp_tools

#utility objects
posc = poscar_ops()
tools = vasp_tools()

#read the poscar file so we know how to differentiate
posc.poscar_reader('POSCAR')

voigt2tensor = [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]

#read and add the data
print('reading the calculated data')
sigma_voigt = np.load('phi_stress_0.00_0.npz')['phi_stress']
sigma_voigt = sigma_voigt + np.load('phi_stress_0.00_1.npz')['phi_stress']
sigma_voigt = sigma_voigt + np.load('phi_stress_0.00_2.npz')['phi_stress']
np.savez_compressed('phi_stress_0.00.npz',phi_stress=sigma_voigt)
'''
#convert to tensor
print('converting to tensorial form')
ilen, jlen, klen, alpha_len = sigma_voigt.shape
sigma_tensor = np.zeros((ilen,jlen,klen,3,3))
for voigt, voigt_pair in enumerate(voigt2tensor):
    ind1, ind2 = voigt_pair
    sigma_tensor[:,:,:,ind1,ind2] = sigma_voigt[:,:,:,voigt]
del sigma_voigt

#take the derivative
print('taking gradient')
nabla_sigma = np.zeros((ilen,jlen,klen,3,3,3))
for ind1 in range(3):
    for ind2 in range(3):
        nabla_sigma[:,:,:,ind1,ind2,:] = tools.get_grad(sigma_tensor[:,:,:,ind1,ind2],posc.box)
del sigma_tensor

#take the divergence
print('taking divergence')
Q = np.einsum('xyzijj->xyzi',nabla_sigma)

np.savez_compressed('Q.npz', Q = Q)'''
print('All done')
