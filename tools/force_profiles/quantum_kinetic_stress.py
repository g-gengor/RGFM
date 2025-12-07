from vaspwfc import vaspwfc
from aewfc import vasp_ae_wfc
import numpy as np
import gc
from poscar_ops import poscar_ops
import sys

#this code follows the formula for stress given in Rappe's paper

posc = poscar_ops()
posc.poscar_reader('POSCAR')
box_vectors = posc.box

thread_num = 64
aecut_num = -4

###DONT CHANGE BELOW
k_id = sys.argv[1]
print('k id: ',k_id)
k_parts = np.loadtxt('k_parts/k_part_{}.txt'.format(int(k_id)))
if len(k_parts.shape) == 0:
    k_parts = [int(k_parts)]
print('k parts: \n', k_parts)

print('gettin ready for partial wfc')
ps_wfc = vaspwfc('../vasp_NV/WAVECAR',omp_num_threads=thread_num)

ae_wfc = vasp_ae_wfc(ps_wfc, aecut=aecut_num,ikpt=1)
print('AE grid shape: ',ae_wfc._aegrid)

k_num = ps_wfc._nkpts
band_num = ps_wfc._nbands

#load the occupancies
#occupancies = np.zeros((k_num,band_num))
occupancies = np.loadtxt('occupancies.txt')
occupancies = occupancies * 2
#for i in range(k_num):
#    data = np.loadtxt('k{}_band_occ'.format(int(i+1)))
#    occupancies[i,:] = data[:,2]
#del data
#gc.collect()

ilen, jlen, klen = ae_wfc._aegrid
#we will write the stress in this vector
sigma = np.zeros((ilen,jlen,klen,6))
#this vector is for density
n = np.zeros((ilen,jlen,klen))


step_x1, stepsize_x1 = np.linspace(0,1.0,endpoint=False,num=ilen,retstep=True)
step_x2, stepsize_x2 = np.linspace(0,1.0,endpoint=False,num=jlen,retstep=True)
step_x3, stepsize_x3 = np.linspace(0,1.0,endpoint=False,num=klen,retstep=True)

h1 = np.linalg.norm(box_vectors[0,:] * stepsize_x1)
h2 = np.linalg.norm(box_vectors[1,:] * stepsize_x2)
h3 = np.linalg.norm(box_vectors[2,:] * stepsize_x3)
n1 = box_vectors[0,:] / np.linalg.norm(box_vectors[0,:])
n2 = box_vectors[1,:] / np.linalg.norm(box_vectors[1,:])
n3 = box_vectors[2,:] / np.linalg.norm(box_vectors[2,:])
G = np.array([[n1[0],n1[1],n1[2]],
                [n2[0],n2[1],n2[2]],
                [n3[0],n3[1],n3[2]]])
G_inv = np.linalg.inv(G)

def nabla_finder(data):
    #nabla_1
    nabla1_psi = (np.roll(data,-1,axis=0) - np.roll(data,1,axis=0)) / (2 * h1)
    #nabla_3
    nabla3_psi = (np.roll(data,-1,axis=2) - np.roll(data,1,axis=2)) / (2 * h3)
    #nabla\psi dot v2
    nabla2_psi = (np.roll(data,-1,axis=1) - np.roll(data,1,axis=1)) / (2 * h2)
    nabla_psi_tens = np.zeros((nabla1_psi.shape[0],nabla1_psi.shape[1],nabla1_psi.shape[2],3))
    nabla_psi_tens[:,:,:,0] = nabla1_psi
    nabla_psi_tens[:,:,:,1] = nabla2_psi
    nabla_psi_tens[:,:,:,2] = nabla3_psi
    nabla_cartesian = np.einsum('ij,xyzj->xyzi',G_inv,nabla_psi_tens)
    return [nabla_cartesian[:,:,:,0],nabla_cartesian[:,:,:,1],nabla_cartesian[:,:,:,2]], [np.conjugate(nabla_cartesian[:,:,:,0]),np.conjugate(nabla_cartesian[:,:,:,1]),np.conjugate(nabla_cartesian[:,:,:,2])]

def voigt2tensor_index(ind):
    if ind == 0:
        return 0, 0
    elif ind == 1:
        return 1, 1
    elif ind == 2:
        return 2, 2
    elif ind == 3:
        return 1, 2
    elif ind == 4:
        return 0, 2
    elif ind == 5:
        return 0, 1

def sigma_kin(beta,k_num=k_num,band_num=band_num):
    phi_stress = np.zeros((ilen,jlen,klen,6))
    n_tot = np.zeros((ilen,jlen,klen))
    for kcounter_dum in k_parts:
        kcounter = int(kcounter_dum)
        print('sigma  loop kcounter: {}'.format(int(kcounter)))
        ae_wfc = vasp_ae_wfc(ps_wfc, aecut=aecut_num,ikpt=kcounter+1)
        for bandcounter in range(band_num):
            if occupancies[kcounter,bandcounter] == 0.0:
                continue
            print('sigma loop bandcounter: {}'.format(int(bandcounter)))
            phi = ae_wfc.get_ae_wfc(iband=bandcounter+1)

            nabla_phi, trash = nabla_finder(phi)
            nabla_phi_conj, trash = nabla_finder(np.conjugate(phi))

            n = np.multiply(np.conjugate(phi),phi)
            n_tot = n_tot + n * occupancies[kcounter,bandcounter]

            nabla_n, trash = nabla_finder(n)
            #now to nabla nabla
            nabla_nabla1_n, trash = nabla_finder(nabla_n[0])
            nabla_nabla2_n, trash = nabla_finder(nabla_n[1])
            nabla_nabla3_n, trash = nabla_finder(nabla_n[2])

            nabla_nabla_n = [nabla_nabla1_n,nabla_nabla2_n,nabla_nabla3_n] #this has the inside index first

            #terms that are dependent on phi only
            for alpha in range(6):
                i, j = voigt2tensor_index(alpha)
                phi_stress[:,:,:,alpha] = phi_stress[:,:,:,alpha] - np.multiply(nabla_phi_conj[i],nabla_phi[j]) * occupancies[kcounter,bandcounter]
                if i == j:
                    phi_stress[:,:,:,alpha] = phi_stress[:,:,:,alpha] + 0.25 * (nabla_nabla_n[0][0] + nabla_nabla_n[1][1] + nabla_nabla_n[2][2]) * occupancies[kcounter,bandcounter]
                #now beta terms
                if beta != 0:
                    phi_stress[:,:,:,alpha] = phi_stress[:,:,:,alpha] + beta * nabla_nabla_n[i][j] * occupancies[kcounter,bandcounter]
                if i == j and beta != 0:
                    phi_stress[:,:,:,alpha] = phi_stress[:,:,:,alpha] - beta * ((nabla_nabla_n[0][0] + nabla_nabla_n[1][1] + nabla_nabla_n[2][2])) * occupancies[kcounter,bandcounter]
    return phi_stress, n_tot
            

beta = 0
phi_stress, n_tot = sigma_kin(beta)
np.savez_compressed('phi_stress_{:.2f}_{}.npz'.format(beta,int(k_id)), phi_stress = phi_stress)
np.savez_compressed('n.npz', n = n_tot)
print('All done!')


