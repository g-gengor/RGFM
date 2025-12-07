import numpy as np

class poscar_ops(object):

    def __init__(self) -> None:
        pass

    def poscar_reader(self,path):
        f = open(path, 'r')
        self.box = np.zeros((3,3))
        atom_count = 0
        for i, line in enumerate(f):
            if i == 0:
                #first line is comment
                comment = line
                continue
            elif i == 1:
                #multiplying factor
                self.multiplier = float(line.split()[0])
                continue
            elif i == 2:
                #first line of box vectors
                vec_string = line.split()
                self.box[0,:] = np.array([float(vec_string[0]), float(vec_string[1]), float(vec_string[2])])
                continue
            elif i == 3:
                #second line of box vectors
                vec_string = line.split()
                self.box[1,:] = np.array([float(vec_string[0]), float(vec_string[1]), float(vec_string[2])])
                continue
            elif i == 4:
                #third line of box vectors
                vec_string = line.split()
                self.box[2,:] = np.array([float(vec_string[0]), float(vec_string[1]), float(vec_string[2])])
                continue
            elif i == 5:
                self.species = line.split()
                continue
            elif i == 6:
                self.species_nums = []
                self.total_atom_num = 0
                for num in line.split():
                    self.species_nums.append(float(num))
                    self.total_atom_num = self.total_atom_num + float(num)
                #self.species_nums = [float(line.split()[0]),float(line.split()[1])]
                self.coords = []
                self.species_counters = []
                for spec in self.species:
                    self.coords.append(np.zeros((int(self.species_nums[0]),3)))
                    self.species_counters.append(0)
                continue
            elif i == 7:
                self.dynamics_flag = line
                continue
            elif i == 8:
                self.coord_mode = line
                continue
            #start reading the coordinates
            line_data = line.split()
            if atom_count < self.total_atom_num:
                coords = np.array([float(line_data[0]), float(line_data[1]), float(line_data[2])])
            else:
                break
            atom_count = atom_count + 1
            for spec_indicator in range(len(self.species)):
                if self.species_counters[spec_indicator] < self.species_nums[spec_indicator]:
                    self.coords[spec_indicator][self.species_counters[spec_indicator],:] = coords
                    self.species_counters[spec_indicator] = self.species_counters[spec_indicator] + 1
                    break
        f.close()
        print('POSCAR read')
    
    def set_structure(self, box, species, coords,coord_mode):
        self.box = box
        self.species = species
        self.coords = coords
        self.coord_mode = coord_mode

    def write_poscar(self,path):
        #currently, only writes selective dynamics, cartesian
        f = open(path,'w')
        f.write('POSCAR file written by GG\n')
        f.write('1\n')
        f.write('{} {} {}\n'.format(self.box[0,0],self.box[0,1],self.box[0,2]))
        f.write('{} {} {}\n'.format(self.box[1,0],self.box[1,1],self.box[1,2]))
        f.write('{} {} {}\n'.format(self.box[2,0],self.box[2,1],self.box[2,2]))
        
        line = ''
        for spec in self.species:
            line = line + spec + ' '
        line = line + '\n'
        f.write(line)
        #write species numbers
        line = ''
        for i, spec in enumerate(self.species):
            line = line + '{} '.format(int(len(self.coords[i][:,0])))
        line = line + '\n'
        f.write(line)

        f.write('Selective Dynamics\n')
        f.write(self.coord_mode+'\n')

        #start writing the coordinates
        for coord_set in self.coords:
            for coord in coord_set:
                line = '{} {} {} T T T\n'.format(coord[0],coord[1],coord[2])
                f.write(line)
        f.close()
        print('POSCAR written')
