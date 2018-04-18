#! /usr/bin/env python

import math
import os
import sys
import numpy as np
import colours
from BigDFT import Logfiles as lf 

# build graphene posinp in the xy plane
def build_graphene(nx,nz):
    xlength=4.61388
    zlength=7.99146
    #C  0.00000  0.00000  0.00000
    #C  0.00000  0.00000  2.66382
    #C  2.30694  0.00000  3.99573
    #C  2.30694  0.00000  6.65955
    posinp=[]
    posinpx=[]
    for i in range(nx):
       posinpx.append([i*xlength+0.0,0.0,0.0])
       posinpx.append([i*xlength+0.0,0.0,2.66382])
       posinpx.append([i*xlength+2.30694,0.0,3.99573])
       posinpx.append([i*xlength+2.30694,0.0,6.65955])

    for k in range(nz):
       for pos in posinpx:
          posinp.append([pos[0],0.0,k*zlength+pos[2]])

    return [xlength*nx,0.0,zlength*nz],posinp




def get_defect_centre(log):
    #this should be generalized somehow - only works for Si substitutional or di-vacancy
    if 'Si' in log.log['posinp']['positions'][0]:
        #Si substitutional
        centre = log.log['posinp']['positions'][0]['Si']
    else:
        #di-vacancy: find centre as the average of first 4 atoms
        #X 26.5317271  9.5000000  27.30427002
        centre = [0.0, 0.0, 0.0]
        for i in range(4):
            for xyz in range(3):
                centre[xyz] = centre[xyz] + log.log['posinp']['positions'][i]['C'][xyz]
        for xyz in range(3):
            centre[xyz] = 0.25*centre[xyz]
    
    #print(centre)
    return centre


# simple criterion - if the closest vertically aligned atom is below, atom is A, otherwise B
# have to look at next closest for defect and edge of cell
# considered vertically aligned if x dist is <0.35
# return true for type a, false for b
def atom_is_a(index,atoms):
    #atoms=list()
    #for position in positions:
    #    atoms.append(position.values()[0])          
    
    #print(atoms)
    xt=0.5

    atom_test=atoms[index]
    #print(atom_test)
       
    lt3=[]
    lt6=[]
    for i, atom in enumerate(atoms):
        #print(atom,atom_test)
        x_dist=atom[0]-atom_test[0]
        if abs(x_dist) < xt:
            z_dist=atom[2]-atom_test[2]
            # first check for immediate above or below (won't exist for atoms at edge of cell)          
            if abs(z_dist)<3.1 and i!=index:
                #print(i,atom,x_dist,z_dist)
                lt3.append(i)
                #print index,i,x_dist,z_dist,xt,3.1
            # for edge of cell cases
            elif abs(z_dist)<6.2 and i!=index:
                lt6.append(i)
    
    #print(lt3)
    #print(lt6)
    #print(len(lt3),len(lt6))
    
    if len(lt3)==1:
        if atoms[lt3[0]][2]-atom_test[2]>0:
            return True
        else:
            return False  
    elif len(lt3)==0:
        # look at lt6 instead
        if len(lt6)==1:
            if atoms[lt6[0]][2]-atom_test[2]<0:
                return True
            else:
                return False
        else:
            print("Error: no lt3 and lt6 of size",len(lt6))        
    else:
        print("Error: lt3 of size and lt6 of size",len(lt3),len(lt6))

# converts from frag list with e.g. 1, ..., 3 to 1, 2, 3
def remove_dots(atoms_dots):
    atoms=list()
    for i,atom_dot in enumerate(atoms_dots):
        if atom_dot != '...':
            atoms.append(atom_dot)
        else:
            first=atoms_dots[i-1]+1
            last=atoms_dots[i+1]
            for j in range(first,last):
                atoms.append(j)            
    #print(atoms_dots)    
    #print(atoms)
    #print("")
    
    return atoms


# This routine returns the distances maximum distance between a given fragment (indicated by n in Cn/Cna/Cnb...)
# and the defect centre, assuming either di-vacancy or Si substitutional
def get_max_distance_from_log(log,index):
    
    #if the fragment index is zero, it must be fragment Si0 and thus the distaance is zero
    if (index==0):
        return 0.0
    
    b2a = 0.529177249
    frags = log.log['frag']
     
    last_atoms = list()
        
    #have to check for a, b, c, d  
    strings = ['C'+str(index)+'f', 'C'+str(index)+'e', 'C'+str(index)+'d', 'C'+str(index)+'c',
	       'C'+str(index)+'b', 'C'+str(index)+'a', 'C'+str(index)]
    for string in strings:
        if string in frags:
            t=frags[string]  
            #last=string
            #print(t)
            #last_atoms.append(t[-1])
            last_atom=t[-1]
	    break
            
    #print(last_atom)
    position = log.log['posinp']['positions'][last_atom-1]['C'] #for last_atom in last_atoms]
    #print(position)
    
    centre = get_defect_centre(log)
    
    #I'm sure there is a function that can simplify life here...
    distance = 0.0
    for xyz in range(3):
        distance = distance + (position[xyz] - centre[xyz])**2
    distance = b2a*np.sqrt(distance)
    #print(distance)  
    
    #print(last_atom,position,distance)

    return distance


# this could be tidied to avoid so much code repetition
def get_all_distances_from_log(log,max_index):
    
    b2a = 0.529177249
    centre = get_defect_centre(log)
       
    frags = log.log['frag']
    #print(len(frags))
    #print(log.log['frag'])

    distances={}
    distances['frag_name']=list()
    distances['frag_num']=list()
    distances['distance']=list()
    distances['index']=list()
    distances['atom_num']=list()
    distances['a_or_b']=list()
    
    at_num=0

    pos=log.log['posinp']['positions']
    pos_list=list()
    for p in pos:
        pos_list.append(p.values()[0])
    
    #first check for Si
    if 'Si0' in frags:     
        #print(frags['Si0'])
        atoms=frags['Si0']
        positions = [log.log['posinp']['positions'][atom-1]['Si'] for atom in atoms]

        #print(positions)
        for a,position in enumerate(positions):
            at_num = at_num + 1 # this one is purely for plotting purposes, not the actual atom number
            distance=0.0
            for xyz in range(3):
                distance = distance + (position[xyz] - centre[xyz])**2
            distance = b2a*np.sqrt(distance)
            distances['distance'].append(distance)
            distances['frag_name'].append('Si0')
            distances['frag_num'].append(0)
            distances['index'].append(at_num)
            distances['atom_num'].append(atoms[a])
            #distances['a_or_b'].append(atom_is_a(atoms[a]-1,log.log['posinp']['positions']))
            distances['a_or_b'].append(atom_is_a(atoms[a]-1,pos_list))
    
    for index in range(1,max_index+1):
        atoms=list()
        #have to check for a, b, c, d
        strings = ['C'+str(index), 'C'+str(index)+'a', 'C'+str(index)+'b',
                   'C'+str(index)+'c','C'+str(index)+'d','C'+str(index)+'e','C'+str(index)+'f']
        for string in strings:
            if string in frags:
                atoms_dots=frags[string]
                #print(atoms_dots)
                atoms=remove_dots(atoms_dots)
                positions = [log.log['posinp']['positions'][atom-1]['C'] for atom in atoms]
                #print(positions)
                for a,position in enumerate(positions):
                    at_num = at_num + 1
                    distance=0.0
                    for xyz in range(3):
                        distance = distance + (position[xyz] - centre[xyz])**2
                    distance = b2a*np.sqrt(distance)
                    distances['distance'].append(distance)
                    distances['frag_name'].append(string)
                    distances['frag_num'].append(index)
                    distances['index'].append(at_num)
                    distances['atom_num'].append(atoms[a])
                    #distances['a_or_b'].append(atom_is_a(atoms[a]-1,log.log['posinp']['positions']))
                    distances['a_or_b'].append(atom_is_a(atoms[a]-1,pos_list))
                    
    #finish with CA and CB
    if 'CA' in frags:     
        atoms_dots=frags['CA']
        atoms=remove_dots(atoms_dots)
        positions = [log.log['posinp']['positions'][atom-1]['C'] for atom in atoms]
        #print(positions)
        for a,position in enumerate(positions):
            at_num = at_num + 1
            distance=0.0
            for xyz in range(3):
                distance = distance + (position[xyz] - centre[xyz])**2
            distance = b2a*np.sqrt(distance)
            distances['distance'].append(distance)
            distances['frag_name'].append('CA')
            distances['frag_num'].append(-1)
            distances['index'].append(at_num)
            distances['atom_num'].append(atoms[a])
            distances['a_or_b'].append(atom_is_a(atoms[a]-1,pos_list))

    #finish with CA and CB
    if 'CB' in frags:     
        atoms_dots=frags['CB']
        atoms=remove_dots(atoms_dots)
        positions = [log.log['posinp']['positions'][atom-1]['C'] for atom in atoms]
        #print(positions)
        for a,position in enumerate(positions):
            at_num = at_num + 1
            distance=0.0
            for xyz in range(3):
                distance = distance + (position[xyz] - centre[xyz])**2
            distance = b2a*np.sqrt(distance)
            distances['distance'].append(distance)
            distances['frag_name'].append('CB')
            distances['frag_num'].append(-2)
            distances['index'].append(at_num)
            distances['atom_num'].append(atoms[a])
            distances['a_or_b'].append(atom_is_a(atoms[a]-1,pos_list)) 
    
    return distances

# Fill in data for energies, errors, fragment distances, DoS etc

#seed=graphene_12x_07z
#   or graphene_12x_07z_d
#   or graphene_12x_07z_2v
def find_data(seed, label, smear):
    #internal to here for now, choose between cubic or linear reference
    #just for calculating errors though, get DoS for both
    ref='l'
    lin_label='linear'
    cub_label='cubic'
    
    Ha2eV = 27.211396132
    data = {'label': label}
    data['frag_max']=[]
    
    filenames=list()
    distances=list()
    # assuming maximum number of files here
    for l in range(100):
        filename="log-"+str(seed)+"_fr"+str(l)+".yaml"
        if (os.path.exists(filename)):
            filenames.append(filename)  
            data['frag_max'].append(l)
     
    filename_lin="log-"+str(seed)+"_fw.yaml"
    filename_cub="log-"+str(seed)+"_cubic.yaml"
    
    data['logfiles'] = [lf.Logfile(filename) for filename in filenames]
    data['index'] = [i for i, filename in enumerate(filenames)]
    data['distances'] = [get_max_distance_from_log(log,data['frag_max'][i])
                         for i, log in enumerate(data['logfiles'])]
    
    #print(paths)
    #print data
    if len(data['logfiles']) < 1:
        raise ValueError("No logfiles found, try a new regexp.")
        
    log_lin = lf.Logfile(filename_lin)
    log_cub = lf.Logfile(filename_cub)
       
    # Find the number of atoms
    data['atoms'] = [log.nat for log in data['logfiles']]
    data['atoms_lin'] = log_lin.nat
    data['atoms_cub'] = log_cub.nat

    # Find the energies/atom wrt ref and convert to eV
    data['energies'] = [Ha2eV * log.energy/log.nat for log in data['logfiles']]
    #print(data['energies'])
    data['energy_lin'] = Ha2eV * log_lin.energy/log_lin.nat
    data['energy_cub'] = Ha2eV * log_cub.energy/log_cub.nat
    data['max_wahba'] = [log.log['Input Hamiltonian']['Maximum Wahba cost function value']
                         for log in data['logfiles']]
    data['av_wahba'] = [log.log['Input Hamiltonian']['Average Wahba cost function value']
                         for log in data['logfiles']]
    
    if (ref=='l'):
        data['errors'] = [Ha2eV * abs(log.energy/log.nat - log_lin.energy/log_lin.nat) \
                      for log in data['logfiles']]
    else:
        data['errors'] = [Ha2eV * abs(log.energy/log.nat - log_cub.energy/log_cub.nat) \
                      for log in data['logfiles']]
    #print data['errors']
    
    npts=2500
    dos_logs = [log.get_dos(label=label, npts=npts) for log in data['logfiles']]
    data['dos_energies'] = [np.array(dos_log.range) for dos_log in dos_logs]
    dos_log_lin = log_lin.get_dos(label=label, npts=npts)
    data['dos_energies_lin'] = np.array(dos_log_lin.range)
    dos_log_cub = log_cub.get_dos(label=label, npts=npts)
    data['dos_energies_cub'] = np.array(dos_log_cub.range)
 
    data['dos'] = [np.array([dos['dos'].curve(dos_log.range, sigma=smear)[1] 
                            for dos in dos_log.ens]).reshape(npts,1) for dos_log in dos_logs]
    data['dos_lin'] = np.array([dos_lin['dos'].curve(dos_log_lin.range, sigma=smear)[1] 
                            for dos_lin in dos_log_lin.ens]).reshape(npts,1)
    data['dos_cub'] = np.array([dos_cub['dos'].curve(dos_log_cub.range, sigma=smear)[1] 
                            for dos_cub in dos_log_cub.ens]).reshape(npts,1)    
    
    # hack to get fermi level for linear - in this case norb = norb_occ * 2
    data['fermi'] = [Ha2eV * log.evals[0][0][num_occ_ks(log)] for log in data['logfiles']]   
    #data['gap'] = [Ha2eV * (log.evals[0][0][ev+1] - log.evals[0][0][ev]) in data['logfiles']]
           
    data['fermi_lin'] = Ha2eV * log_lin.evals[0][0][num_occ_ks(log)]
    #data['gap_ref'] = Ha2eV * (log_ref.evals[0][0][ev+1] - log_ref.evals[0][0][ev])
    data['fermi_cub'] = dos_log_cub.ef
    # leave cubic gap unset for the mo
    #data['gap_ref'] = -1    
    
    return data



#hack for getting Fermi level, only works for Si substitutional, di-vacancy and pristine for 336/334/336 atoms
def num_occ_ks(log):
    if 'Si' in log.log['posinp']['positions'][0]:
        return 672 - 1
    elif log.nat==336:
	return 672 -1
    else:
        return 668 - 1

# function for making v_sim parameter file with atoms in the appropriate colours
def write_vsim_file(nfrag,nfrag_max):
    
    filename="v_sim.res"
    print "writing to ",filename
    ff=open(filename,'w')

    ff.write("#V_Sim resources file v3.0\n")
    ff.write("#====================\n")

    ff.write("# Set the background of the background ; four floating point values (0. <= v <= 1.)\n")
    ff.write("backgroundColor_color:\n")
    ff.write("    1.000 1.000 1.000 1.000\n")
    ff.write("# Control if a box is drawn around the rendering area ; boolean (0 or 1)\n")
    ff.write("box_is_on:\n")
    ff.write("    1\n")
    ff.write("# Define the color of the box ; three floating point values (0. <= v <= 1.)\n")
    ff.write("box_color:\n")
    ff.write("    0.000 0.000 0.000\n")
    ff.write("# Define the width of the lines of the box ; one integer (1. <= v <= 10.)\n")
    ff.write("box_line_width:\n")
    ff.write("       2\n")
    ff.write("# Control if the legend is drawn ; boolean (0 or 1)\n")
    ff.write("legend_is_on:\n")
    ff.write("    0\n")

    
    ff.write("# The radius of the element and its shape, a real > 0. & [Sphere Cube Elipsoid Point]\n")
    for nf in range(nfrag):
        ff.write("atomic_radius_shape:\n")
        ff.write("    C"+str(nf+1)+"  1.000 Sphere\n")   
    ff.write("atomic_radius_shape:\n")
    ff.write("    CA    0.95 Sphere\n")
    ff.write("atomic_radius_shape:\n")
    ff.write("    CB    0.95 Sphere\n")
    
    colors = colours.find_colours(nfrag_max)['indices']
    
    #print colors['indices']
    #print colors['html']
    
    # need to add if for Si
    
    ff.write("# Codes the main color in RedGreenBlueAlpha formatand the light effects on material, nine floats between 0. and 1.\n")
    for nf in range(nfrag):
        #print colors[nf]
        ff.write("element_color:\n")
        ff.write("    C"+str(nf+1)+" "+str(colors[nf][0]/255.0)+" "\
                 +str(colors[nf][1]/255.0)+" "\
                 +str(colors[nf][2]/255.0)+" "\
                 +" 1.000   0.25 0.25 0.25 0.25 0.25\n")
        #element_color:
        #    C1 0.502 0.000 0.000 1.000   0.25 0.25 0.25 0.25 0.25
        #element_color:
        #    C2 1.000 0.000 0.000 1.000   0.25 0.25 0.25 0.25 0.25
        
    #'#8f8f8f', '#4f4f4f'
    cab=[143, 79]
    ff.write("element_color:\n")
    ff.write("    CA "+str(cab[0]/255.0)+" "\
            +str(cab[0]/255.0)+" "\
            +str(cab[0]/255.0)+" "\
            +" 1.000   0.25 0.25 0.25 0.25 0.25\n")
    ff.write("element_color:\n")
    ff.write("    CB "+str(cab[1]/255.0)+" "\
            +str(cab[1]/255.0)+" "\
            +str(cab[1]/255.0)+" "\
            +" 1.000   0.25 0.25 0.25 0.25 0.25\n")
    
    
    ff.write("# This value is the width for all pairs drawn ; 0 < integer < 10\n")
    ff.write("pairWire_width:\n")
    ff.write("    2\n")
    ff.write("# Widths detail for each drawn link ; 0 < integer < 10\n")

    ff.write("# It chooses the colors of the cylinders according differents criterion ; 0 <= integer < 2\n")
    ff.write("cylinder_colorType:\n")
    ff.write("    1\n")
    ff.write("# This value is the default radius of the pairs drawn as cylinders ; 0 < real < 10\n")
    ff.write("pairCylinder_radius:\n")
    ff.write("    0.250000\n")
    ff.write("# This value is the radius for specific pairs drawn as cylinders ; element1 elemen2 0 < real < 10\n")

    ff.write("# Ask the opengl engine to draw pairs between elements ; boolean 0 or 1\n")
    ff.write("pairs_are_on:\n")
    ff.write("    1\n")
    ff.write("# Favorite method used to render files ; chain ('Wire pairs', 'Cylinder pairs')\n")
    ff.write("pairs_favoriteMethod:\n")
    ff.write("    Cylinder pairs\n")
    ff.write("# Draw a link between [ele1] [ele2] [0. <= dmin] [0. <= dmax]\n")
    ff.write("#                     [0. <= RGB <= 1.]x3 [bool: drawn] [bool: printLength] [string: method]\n")

    for nf in range(nfrag):
        ff.write("pair_link:\n")
        ff.write("    CA C"+str(nf+1)+" 2.500 3.000\n")
        ff.write("    1.000 0.600 0.200  1  0  Cylinder pairs\n")
        ff.write("pair_link:\n")
        ff.write("    CB C"+str(nf+1)+" 2.500 3.000\n")
        ff.write("    1.000 0.600 0.200  1  0  Cylinder pairs\n")
        for mf in range(nfrag):
            ff.write("pair_link:\n")
            ff.write("    C"+str(mf+1)+" C"+str(nf+1)+" 2.500 3.000\n")
            ff.write("    1.000 0.600 0.200  1  0  Cylinder pairs\n")
    
    ff.close()


def read_overlap_onsite(filename):
    f=open(filename,'r')

    oo={}
    oo['s']=list()
    oo['px']=list()
    oo['py']=list()
    oo['pz']=list()

    for l in f:
        field=l.split()
        # skip commented lines
        if field[0] == "#":
            continue
        # only read in if the tmbs are of the same type, i.e. s and s (assume no internal rotation between px,y,z)
        # store atom numbers and overlap, don't need tmb number
        if int(field[1])%4 == 1 and int(field[0])%4 == 1:
            oo['s'].append([field[3],field[4],field[2]])
        elif int(field[1])%4 == 2 and int(field[0])%4 == 2:
            oo['px'].append([field[3],field[4],field[2]])
        elif int(field[1])%4 == 3 and int(field[0])%4 == 3:
            oo['py'].append([field[3],field[4],field[2]])
        elif int(field[1])%4 == 0 and int(field[0])%4 == 0:
            oo['pz'].append([field[3],field[4],field[2]])

    f.close()

    return oo


















