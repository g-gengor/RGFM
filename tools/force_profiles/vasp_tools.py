import sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import scipy
from poscar_ops import poscar_ops
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import cm
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

class vasp_tools(object):
    
    def __init__(self) -> None:
        self.Eh2eV = 27.2114
        self.a02A = 0.529177
        self.Ha_a0_2TPa = 2.9421015697 * 10
    def read_density_file(self,path):
        #read the locpot
        #MAKE SURE YOU HAVE DELETED THE SIMBOX INFO AT THE BEGINNING OF LOCPOT
        #Make a table, elements separated by whitespace,
        #Always five columns
        df = pd.read_table(path, sep='\s+', names=range(5))
        #print(df)
        #check if the header is removed
        print('comment line check: ', df[0].str.isnumeric()[0]==False)#df.isnull().values[0][1])
        if df[0].str.isnumeric()[0]==False:#df.isnull().values[0][1] == True:
            #get the number of atoms
            atom_nums = df.iloc[6].values
            total_num = 0
            for num in atom_nums:
                if np.isnan(float(num)) == False:
                    total_num = total_num + float(num)
            print('total atom number: ', total_num)
            drop_row_lim = 7 + total_num
            drop_list = []
            for i in range(int(drop_row_lim+1)):
                #print('dropping row', i)
                drop_list.append(i)
            df.drop(drop_list,inplace=True)
        #df = update_df
        #print(df)
        #Grab the dimensions, first three elements of the first row
        dims=tuple(int(d) for d in df.iloc[0][:3])
        print('tesselation: ', dims)

        #Find number of missing values, displayed as nan
        #remove two for "missing values" in dimension row
        nan_count=df.isna().sum().sum()-2
        print('nan count: ', nan_count)
        #Convert to 1D numpy array, removing excess nans/dimensions
        if nan_count==0:
            #D=df.to_numpy().flatten()[5:]
            D=df.values.flatten()[5:]
        else:
            D=df.to_numpy().flatten()[5:-nan_count]
        D = D.astype(float)
        #Reshape to 3D array
        V=D.reshape(dims,order='F')#eV
        #print(D)
        del df
        del D
        
        return V

    def get_grad(self,data,box):
        total_steps = data.shape
        total_step_x1, total_step_x2, total_step_x3 = total_steps
        #print(total_steps)
        step_x1, stepsize_x1 = np.linspace(0,1.0,endpoint=False,num=total_step_x1,retstep=True)
        step_x2, stepsize_x2 = np.linspace(0,1.0,endpoint=False,num=total_step_x2,retstep=True)
        step_x3, stepsize_x3 = np.linspace(0,1.0,endpoint=False,num=total_step_x3,retstep=True)
        
        h1 = np.linalg.norm(box[0,:] * stepsize_x1)
        h2 = np.linalg.norm(box[1,:] * stepsize_x2)
        h3 = np.linalg.norm(box[2,:] * stepsize_x3)
        n1 = box[0,:] / np.linalg.norm(box[0,:])
        n2 = box[1,:] / np.linalg.norm(box[1,:])
        n3 = box[2,:] / np.linalg.norm(box[2,:])
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
            return nabla_cartesian
        return nabla_finder(data)
    
    def calcualte_F(self,path,box):
        V = self.read_density_file(path)#eV
        V = V / self.Eh2eV #Ha
        F = self.get_grad(-V,box)
        return F #Ha/e*a0
    def pbc_enforcer(self,r,box_vectors):
        #r = A * f, f is fractional coords 
        A = np.transpose(box_vectors)
        #f = B * r, B is A^-1
        B = np.linalg.inv(A)
        resx = 0.0
        resy = 0.0
        resz = 0.0
        f = np.dot(B,r)
        g = f - np.floor(f)
        t = np.dot(A,g)
        
        return t

    def shepards_pointwise(self,V,box_vectors,rs):
        #total_step_x1,total_step_x2,total_step_x3 = target_shape
        ilen,jlen,klen = V.shape
        V_coordinates_x,V_coordinates_y,V_coordinates_z = np.meshgrid(np.arange(ilen),np.arange(jlen),np.arange(klen),indexing='ij')
        V_coordinates_x = np.ravel(V_coordinates_x,order='C').reshape(-1,1)
        V_coordinates_y = np.ravel(V_coordinates_y,order='C').reshape(-1,1)
        V_coordinates_z = np.ravel(V_coordinates_z,order='C').reshape(-1,1)
        V_coordinates_proto = np.multiply(V_coordinates_x, box_vectors[0,:]) / ilen +\
                        np.multiply(V_coordinates_y, box_vectors[1,:]) / jlen +\
                        np.multiply(V_coordinates_z, box_vectors[2,:]) / klen
        V_coordinates = np.copy(V_coordinates_proto)
        V_flattend = np.ravel(V,order='C')
        V_flat = np.copy(V_flattend)
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    V_flat = np.append(V_flat,V_flattend,axis=0)
                    V_coordinates = np.append(V_coordinates,V_coordinates_proto,axis=0)
        #form the KD tree for neighbor lookup
        tree = scipy.spatial.KDTree(V_coordinates)
        
        #now we are going to apply modified inverse distance weghing interpolation
        #https://en.wikipedia.org/wiki/Inverse_distance_weighting
        def weight_function(query_point,tree,radius):
            neigh = tree.query_ball_point(query_point,radius)
            if len(neigh) == 0:
                print('error no neighbours')
                print(query_point)
                raise SystemExit("Stop right there!")
            w = np.zeros((len(neigh),),dtype='float')
            dists = np.linalg.norm(tree.data[neigh] - query_point,axis=1)
            if np.any(dists==0.0):
                w[np.where(dists==0.0)] = 1.0
                return w, neigh
            w = np.square((np.maximum(0.0,radius - dists) / (radius * dists)))
            return w, neigh
        
        
        psi_coordinates = rs
        
        V_interpolated = np.zeros((len(psi_coordinates[:,0]),))
        radius = 0.12
        stopcount = 0
        for i in range(len(psi_coordinates[:,0])):
            psi_reduced = self.pbc_enforcer(psi_coordinates[i,:],box_vectors)
            weights, neighs = weight_function(psi_reduced,tree,radius)
            if len(weights) == 0:
                print('No neighbor situation: ', i)
            V_interpolated[i] = np.dot(V_flat[neighs], weights / sum(weights))
        return V_interpolated
    
    def plot_plane(self,data,v1,v2,plane_normal,plane_offset,atom_offset,poscar,tesselation_step,box_vectors,atomic_radii,atomic_colors,colormap_type='seismic',colormap_lb=0,colormap_ub=1.0,norm=None):
        #plane_offset moves plane for point generation
        #atom offset choses the atoms on a specific plane so that plane normal * atom pos gives this number
        print('cmap: ', colormap_type)
        #read the DISP file
        df = pd.read_table('DISP', sep='\s+', names=range(8))#id type x y z DisplacementX DisplacementY DisplacementZ
        disp_data = df.to_numpy()
        #plot plane parameters
        e3 = np.cross(v1,v2)
        e3 = e3 / np.linalg.norm(e3)
        
        e1 = v1 / np.linalg.norm(v1)
        e2 = v2 / np.linalg.norm(v2)
        
        R = np.einsum('i,j->ij',e1,np.array([1,0,0])) + np.einsum('i,j->ij',e2,np.array([0,1,0])) + np.einsum('i,j->ij',e3,np.array([0,0,1]))
        #create the interpolation grid
        step_number = tesselation_step
        steps = np.linspace(0,1.0,step_number,endpoint=False)
        points = []
        offset = e3 * plane_offset#np.zeros((3))#np.array([a/2,0,a/2])
        periodic_points = []
        for steps1 in steps:
            for steps2 in steps:
                point = steps1 * v1 + steps2 * v2
                point = point + offset
                points.append(point)
                point = self.pbc_enforcer(point,box_vectors)
                periodic_points.append(point)
        points = np.array(points)
        periodic_points = np.array(periodic_points)
        #contour points to be used in plotting
        ccc = np.einsum('ij,xj->xi',R.T,points)
        #interpolate
        data_interp = self.shepards_pointwise(data,poscar.box,points)
        cmap = plt.get_cmap(colormap_type)

        #plotting
        qfnorm2plot = data_interp
        print('range max: ',np.max(qfnorm2plot))
        print('range min: ',np.min(qfnorm2plot))
        plt.close('all')
        fig1, ax1 = plt.subplots()
        fig1.set_dpi(100)
        
        #norm= cm.colors.Normalize(vmin=np.min(qfnorm2plot), vmax=np.max(qfnorm2plot))
        
        #norm= cm.colors.CenteredNorm(halfrange=3.4)
        if norm == None:
            norm = cm.colors.CenteredNorm()
        new_cmap = self.truncate_colormap(cmap, colormap_lb, colormap_ub)
        plt.tricontourf(ccc[:,0],ccc[:,1],qfnorm2plot,100,cmap=new_cmap,norm=norm)
        
        for atomic_disp in disp_data:
            pos = atomic_disp[2:5]
            if np.abs(np.dot(pos,plane_normal) - atom_offset) > 0.1:
                continue
            print('pos: ', atomic_disp[2:5])
            species_code = atomic_disp[1]
            r = atomic_radii[int(species_code)-1]
            color = atomic_colors[int(species_code)-1]
                
            coords = np.zeros((2))
            coords[0] = np.einsum('i,i->',atomic_disp[2:5],e1)
            coords[1] = np.einsum('i,i->',atomic_disp[2:5],e2)
            cc = plt.Circle((coords[0],coords[1]),r,color=color,fill=False)#,ls='--')
            
            ax1.add_patch(cc)
        scalebar = AnchoredSizeBar(ax1.transData,
                                1, r'1 $\AA$', 'lower center', 
                                pad=0.1,
                                color='black',
                                frameon=False,
                                size_vertical=0.1,)
        
        ax1.add_artist(scalebar)
        
        sm = plt.cm.ScalarMappable(norm=norm, cmap = new_cmap)
        cbar = fig1.colorbar(sm)
        cbar.formatter.set_useMathText(True)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.ax.yaxis.get_offset_text().set(size=18)
        cbar.ax.tick_params(labelsize=18)
        #plt.colorbar()
        #ax1.axis('equal')
        #plt.colorbar()
        plt.axis('off')
        return fig1, ax1

    def truncate_colormap(self,cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    def find_atoms(self,plane_normal,atom_offset):
        #read the DISP file
        df = pd.read_table('DISP', sep='\s+', names=range(8))#id type x y z DisplacementX DisplacementY DisplacementZ
        disp_data = df.to_numpy()

        res_data = []
        for atomic_disp in disp_data:
            pos = atomic_disp[2:5]
            if np.abs(np.dot(pos,plane_normal) - atom_offset) > 0.1:
                continue
            #print('pos: ', atomic_disp[2:5])
            res_data.append(atomic_disp)

        res_data = np.array(res_data)
        return res_data

    def voronoi_integrate(self,data,box_vectors):
        #box is taken in A^3 but the result has a0^3
        #read the DISP file
        df = pd.read_table('DISP', sep='\s+', names=range(8))#id type x y z DisplacementX DisplacementY DisplacementZ
        disp_data = df.to_numpy()
        #We are adding some images so that we can get full tesselation in the primary simulation cell
        adding_array = [-1.0,0.0,1.0]
        atomic_pos_pb = disp_data[:,2:5]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    atomic_pos_temp = disp_data[:,2:5] + adding_array[i] * box_vectors[0,:] + adding_array[j] * box_vectors[1,:] + adding_array[k] * box_vectors[2,:]
                    np.append(atomic_pos_pb,atomic_pos_temp,axis=0)
        #use the KD tree algorithm to decide which grid point belongs to which atom
        #form the KD tree for neighbor lookup
        tree = scipy.spatial.KDTree(atomic_pos_pb)

        #volume per voxel
        dV = np.linalg.det(box_vectors/self.a02A) / (data.shape[0] * data.shape[1] * data.shape[2]) 
        #create the vector to hold the strains
        vor_integs = np.zeros((len(disp_data[:,0]),2))#this will hold the sum in first index and total number of voxels in the second
        #create the grid points
        ilen,jlen,klen = data.shape
        V_coordinates_x,V_coordinates_y,V_coordinates_z = np.meshgrid(np.arange(ilen),np.arange(jlen),np.arange(klen),indexing='ij')
        V_coordinates_x = np.ravel(V_coordinates_x,order='C').reshape(-1,1)
        V_coordinates_y = np.ravel(V_coordinates_y,order='C').reshape(-1,1)
        V_coordinates_z = np.ravel(V_coordinates_z,order='C').reshape(-1,1)
        points = np.multiply(V_coordinates_x, box_vectors[0,:]) / ilen +\
                        np.multiply(V_coordinates_y, box_vectors[1,:]) / jlen +\
                        np.multiply(V_coordinates_z, box_vectors[2,:]) / klen
        #associate every grid point with an atom
        data_flat = np.ravel(data,order='C')
        for point_index in range(len(points[:,0])):
            point = points[point_index,:]
            dd, ii = tree.query(point,k=1)
            #print(dd)
            #print(ii)
            #if point_index == 3:
            #    break
            index = ii % len(disp_data[:,0])
            vor_integs[index,0] = vor_integs[index,0] + data_flat[point_index] #/ (data.shape[0] * data.shape[1] * data.shape[2]) #integrate the data
            vor_integs[index,1] = vor_integs[index,1] + 1 #keep track of total volume of every voxel
        return vor_integs       