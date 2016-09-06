
# coding: utf-8

# In[33]:

#Marcin Miklitz 13 Oct 2015 marcin.miklitz(at)gmail.com, m.miklitz14@imperial.ac.uk
import numpy as np
import glob
import operator as op
import mmap
import contextlib as ct
import sys
import time
import multiprocessing as mp
import matplotlib.pyplot as plt
import copy
import scipy
from mpl_toolkits.mplot3d import axes3d, Axes3D
from sklearn.neighbors import KDTree
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.metrics.pairwise import euclidean_distances
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

#Atom mass dictionary taken from www.ccdc.cam.ac.uk/Lists/ResourceFileList/Elemental_Radii.xlsx 13 Oct 2015
#Excluding unstable radioisotopes, dummy atom denoted X (and D) have atomic weight equal 1
atom_mass = {
             'AL': 26.982,  'SB': 121.76,  'AR': 39.948,  'AS': 74.922,  'BA': 137.327, 'BE': 9.012,   'BI': 208.98, 
             'B':  10.811,  'BR': 79.904,  'CD': 112.411, 'CS': 132.905, 'CA': 40.078,  'C':  12.011,  'CE': 140.116, 
             'CL': 35.453,  'CR': 51.996,  'CO': 58.933,  'CU': 63.546,  'DY': 162.5,   'ER': 167.26,  'EU': 151.964,
             'F':  18.998,  'GD': 157.25,  'GA': 69.723,  'GE': 72.61,   'AU': 196.967, 'HF': 178.49,  'HE': 4.003,
             'HO': 164.93,  'H':  1.008,   'IN': 114.818, 'I':  126.904, 'IR': 192.217, 'FE': 55.845,  'KR': 83.8,
             'LA': 138.906, 'PB': 207.2,   'LI': 6.941,   'LU': 174.967, 'MG': 24.305,  'MN': 54.938,  'HG': 200.59,
             'MO': 95.94,   'ND': 144.24,  'NE': 20.18,   'NI': 58.693,  'NB': 92.906,  'N':  14.007,  'OS': 190.23,
             'O':  15.999,  'PD': 106.42,  'P':  30.974,  'PT': 195.078, 'K':  39.098,  'PR': 140.908, 'PA': 231.036,
             'RE': 186.207, 'RH': 102.906, 'RB': 85.468,  'RU': 101.07,  'SM': 150.36,  'SC': 44.956,  'SE': 78.96,
             'SI': 28.086,  'AG': 107.868, 'NA': 22.991,  'SR': 87.62,   'S':  32.066,  'TA': 180.948, 'TE': 127.6,
             'TB': 158.925, 'TL': 204.383, 'TH': 232.038, 'TM': 168.934, 'SN': 118.71,  'TI': 47.867,  'W':  183.84,
             'U':  238.029, 'V':  50.942,  'XE': 131.29,  'YB': 173.04,  'Y':  88.906,  'ZN': 65.39,   'ZR': 91.224,
             'X':  1.00000, 'D':  1.0000
            }

#Atom vdW radii dictionary taken from www.ccdc.cam.ac.uk/Lists/ResourceFileList/Elemental_Radii.xlsx 13 Oct 2015
#Excluding unstable radioisotopes, dummy atom denoted X (and D) have atomic vdW radii equal 1
atom_vdw_radii = {
                  'AL': 2,    'SB': 2,    'AR': 1.88, 'AS': 1.85, 'BA': 2,    'BE': 2,    'BI': 2, 
                  'B':  2,    'BR': 1.85, 'CD': 1.58, 'CS': 2,    'CA': 2,    'C':  1.7,  'CE': 2, 
                  'CL': 1.75, 'CR': 2,    'CO': 2,    'CU': 1.4,  'DY': 2,    'ER': 2,    'EU': 2,
                  'F':  1.47, 'GD': 2,    'GA': 1.87, 'GE': 2,    'AU': 1.66, 'HF': 2,    'HE': 1.4,
                  'HO': 2,    'H':  1.09, 'IN': 1.93, 'I':  1.98, 'IR': 2,    'FE': 2,    'KR': 2.02,
                  'LA': 2,    'PB': 2.02, 'LI': 1.82, 'LU': 2,    'MG': 1.73, 'MN': 2,    'HG': 1.55,
                  'MO': 2,    'ND': 2,    'NE': 1.54, 'NI': 1.63, 'NB': 2,    'N':  1.55, 'OS': 2,
                  'O':  1.52, 'PD': 1.63, 'P':  1.8,  'PT': 1.72, 'K':  2.75, 'PR': 2,    'PA': 2,
                  'RE': 2,    'RH': 2,    'RB': 2,    'RU': 2,    'SM': 2,    'SC': 2,    'SE': 1.9,
                  'SI': 2.1,  'AG': 1.72, 'NA': 2.27, 'SR': 2,    'S':  1.8,  'TA': 2,    'TE': 2.06,
                  'TB': 2,    'TL': 1.96, 'TH': 2,    'TM': 2,    'SN': 2.17, 'TI': 2,    'W':  2,
                  'U':  1.86, 'V':  2,    'XE': 2.16, 'YB': 2,    'Y':  2,    'ZN': 1.29, 'ZR': 2,
                  'X':  1.0,  'D':  1.0
                 }
    
#Atom covalent radii dictionary taken from www.ccdc.cam.ac.uk/Lists/ResourceFileList/Elemental_Radii.xlsx 13 Oct 2015
#Excluding unstable radioisotopes, dummy atom denoted X (and D) have atomic cov radii equal 1
atom_cov_radii = {
                  'AL': 1.21, 'SB': 1.39, 'AR': 1.51, 'AS': 1.21, 'BA': 2.15, 'BE': 0.96, 'BI': 1.48, 
                  'B':  0.83, 'BR': 1.21, 'CD': 1.54, 'CS': 2.44, 'CA': 1.76, 'C':  0.68, 'CE': 2.04, 
                  'CL': 0.99, 'CR': 1.39, 'CO': 1.26, 'CU': 1.32, 'DY': 1.92, 'ER': 1.89, 'EU': 1.98,
                  'F':  0.64, 'GD': 1.96, 'GA': 1.22, 'GE': 1.17, 'AU': 1.36, 'HF': 1.75, 'HE': 1.5,
                  'HO': 1.92, 'H':  0.23, 'IN': 1.42, 'I':  1.4,  'IR': 1.41, 'FE': 1.52, 'KR': 1.5,
                  'LA': 2.07, 'PB': 1.46, 'LI': 1.28, 'LU': 1.87, 'MG': 1.41, 'MN': 1.61, 'HG': 1.32,
                  'MO': 1.54, 'ND': 2.01, 'NE': 1.5,  'NI': 1.24, 'NB': 1.64, 'N':  0.68, 'OS': 1.44,
                  'O':  0.68, 'PD': 1.39, 'P':  1.05, 'PT': 1.36, 'K':  2.03, 'PR': 2.03, 'PA': 2,
                  'RE': 1.51, 'RH': 1.42, 'RB': 2.2,  'RU': 1.46, 'SM': 1.98, 'SC': 1.7,  'SE': 1.22,
                  'SI': 1.2,  'AG': 1.45, 'NA': 1.66, 'SR': 1.95, 'S':  1.02, 'TA': 1.7,  'TE': 1.47,
                  'TB': 1.94, 'TL': 1.45, 'TH': 2.06, 'TM': 1.9,  'SN': 1.39, 'TI': 1.6,  'W': 1.62,
                  'U':  1.96, 'V':  1.53, 'XE': 1.5,  'YB': 1.87, 'Y':  1.9,  'ZN': 1.22, 'ZR': 1.75,
                  'X':  1.0,  'D':  1.0
                 }

#Atom keys for deciphering OPLS_2005 taken from DL_FIELD 3.3 dl_field.atom_type
opls_atom_keys = {'C': ['CA', 'CR3', 'CT', 'CD', 'C3T', 'C4T', 'CA5', 'CM', 'CT1', 'C', 'CZB', 'CZN', 'CO3', 'CG1',
                        'CO4', 'CALK', 'CQ', 'C5A', 'C5BC', 'CANI', 'CZ', 'CT1G', 'CTA', 'CTD', 'CTE', 'CTG', 'CTI',
                        'CTJ', 'CTM', 'CTP', 'CTQ', 'CTS', 'CTU', 'CP1', 'CTH', 'C5M', 'C5N', 'CN', 'CB', 'CTC', 
                        'CRA', 'C5X', 'CT1', 'C5B', 'CRA'],
                  'H': ['H', 'HE', 'HC', 'HAE', 'HANI', 'HMET', 'HA', 'HO', 'HS', 'H_OH'],
                  'O': ['O3T', 'O4T', 'OS', 'OVE', 'OHP', 'O', 'OH', 'OAL', 'O2Z', 'OES', 'OA', 'OM', 'ON', 'OY', 
                        'OZ'],
                  'N': ['NA', 'NA5B', 'NA5', 'N5B', 'N', 'NE', 'NI', 'NZB', 'NG', 'NP', 'NB', 'NZ', 'NOM', 'NO2', 
                        'NT', 'NE1'],
                  'F': ['F', 'FG'],
                  'S': ['S', 'SH', 'SY', 'SZ', 'SX6']}

def uniq(input):
    """
    This function removes duplicates from a list and returns it.
    """
    output = []
    for x in input:
        if x not in output:
            output.append(x)
    return(output)

def two_points_distance(point_a,point_b):
    """
    This is a function that takes XYZ positions for two points A and B and calculates the distance between them.
    This function assumes a list as an input of kind [x, y, z] and x, y, z should be floats.
    The point can be atoms XYZ coordinates, or center of mass XYZ coordinates, therefore you can calculate:
    atom-atom, com-atom, com-com distances. This equation is faster than invoking numpy.linalg
    Output: float
    """
    return(((point_a[0]-(point_b[0]))**2 + (point_a[1]-(point_b[1]))**2 + (point_a[2]-(point_b[2]))**2)**0.5)

def center_of_mass(atom_list):
    """
    This function calculates the centre of mass (COM) of a given list of atoms. It requires a list of type:
    [id, no, x, y, z], where id is the atom key identifier (It should be in uppercase and for force fielellipds
    atom keys, they should already be deciphered with an appropriate function. Periodic table notation is required.)
    no is an index number (not relevant for this function) and x, y, z are xyz atom coordinates, respectively,
    and should be floats.
    Output: list [a, b, c] where a, b, c are COM coordinates
    """
    mass = 0
    mass_x = 0
    mass_y = 0
    mass_z = 0
    for i in atom_list:
        if i[0] in atom_mass:
            mass += atom_mass[i[0]]                 #Total mass of the compound
            mass_x += atom_mass[i[0]] * i[2]        #Coordinate multiplied by the atom mass
            mass_y += atom_mass[i[0]] * i[3]        # -//-
            mass_z += atom_mass[i[0]] * i[4]        # -//-
        else:                                       #Atom key is not recognised and is missing from dictionary
            print('No mass for atom {1} in atom mass dictionary'.format(i[0]))
    return([mass_x/mass, mass_y/mass, mass_z/mass]) #Each mass*coordinate has to be divided by total mass of compound

def center_of_geometry(atom_list):
    """
    This function calculates the centre of geometry (COG) of a given list of atoms. It requires a list of type:
    [id, no, x, y, z], where id is the atom key identifier (It should be in uppercase and for force fields
    atom keys, they should already be deciphered with an appropriate function. Periodic table notation is required.)
    no is an index number (not relevant for this function) and x, y, z are xyz atom coordinates, respectively,
    and should be floats.
    Output: list [a, b, c] where a, b, c are COG coordinates
    """
    no_of_atoms = len(atom_list)                    #We need total number of molecules to take an avarage
    coordinate_a = 0
    coordinate_b = 0
    coordinate_c = 0
    for i in atom_list:
        coordinate_a += i[2]                        #We sum all coordinates for each component x,y,z
        coordinate_b += i[3]
        coordinate_c += i[4]                        #Each sum needs to be avaraged over all atoms
    return([coordinate_a/no_of_atoms, coordinate_b/no_of_atoms, coordinate_c/no_of_atoms])

def com2zero(atom_list):
    """
    This function first calculates the center of mass of the molecule and than translate all the coordinates
    So that new center of mass ends up in origin.
    """
    com = center_of_mass(atom_list) 
    return([[i[0],i[1],i[2]-com[0],i[3]-com[1],i[4]-com[2]] for i in atom_list])


def max_dim(atom_list, atom_coord, atom_vdw_vertical, atom_vdw_horizontal):
    """
    New version (01-12-15) of maximum dimension function. It now uses distance matrices, but still the data
    on which atom creates the maximum distance is retained and corrected by vdw radii (120x faster then previous)
    """
    dist_matrix = euclidean_distances(atom_coord, atom_coord)
    vdw_matrix = atom_vdw_vertical + atom_vdw_horizontal
    re_dist_matrix = dist_matrix + vdw_matrix
    final_matrix = np.triu(re_dist_matrix)
    i,j = np.unravel_index(final_matrix.argmax(), final_matrix.shape)
    answer = [atom_list[i][:2],atom_list[j][:2],final_matrix[i,j]]
    return(answer[2])

def void_diameter(atom_list, atom_coord, atom_vdw):
    """
    New version (01-12-15) of void diameter function. It now uses distance matrices, but still the data
    on which atom is closest to the COM is retained and corrected by vdw radii (230x faster then previous)
    """
    dist_matrix = euclidean_distances(atom_coord, [0,0,0])
    new_dist_matrix = dist_matrix - atom_vdw
    answer = [atom_list[np.argmin(new_dist_matrix)][:2],new_dist_matrix[np.argmin(new_dist_matrix)][0]*2]
    return(answer[1])

#Calculate normal vector for a provided set: COM(centre) of this plane and two vectors for two points on this plane 
def normal_vector(x,y):
    va = (y[0][2:])
    vb = (y[1][2:])
    v1 = np.array([va[0]-x[0],va[1]-x[1],va[2]-x[2]])
    v2 = np.array([vb[0]-x[0],vb[1]-x[1],vb[2]-x[2]])
    v3 = np.cross(v1,v2)
    return v3

#Calculates the angle between two normal vectors
def vec_angle(x,y):
    first_step = abs(x[0]*y[0] + x[1]*y[1] + x[2]*y[2]) / (np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)*                                                           np.sqrt(y[0]**2 + y[1]**2 + y[2]**2))
    first_step_prim = abs(x[0]*y[0] + x[1]*y[1] + x[2]*y[2]) /                       (x[0]**2 + x[1]**2 + x[2]**2)**0.5 * (y[0]**2 + y[1]**2 + y[2]**2)**0.5
    second_step = np.arccos(first_step)
    third_step = np.rad2deg(second_step)
    return(second_step)

#Calculates the angle between two planes (between guest benzene ring and all host benzene rings and windows)
def angle(x,y):
    CA_list_benzene = []
    for i in y:
        if i[0].upper() == 'CA':
            CA_list_benzene.append(i)
    com_benzene_guest = list(com(CA_list_benzene))
    nv_guest = normal_vector(com_benzene_guest,CA_list_benzene)
    CA_list_host = []
    for i in x:
        if i[0].upper() == 'CA':
            CA_list_host.append(i)
    benzene_rings = []
    itters = 0
    while len(CA_list_host) > 0:
        temp_ring_list = []
        for i in CA_list_host:
            if atom_distance(CA_list_host[0],i) < 3:
                temp_ring_list.append(i)
        for i in temp_ring_list:
            CA_list_host.remove(i)
        benzene_rings.append(temp_ring_list)
    com_list = []
    for i in benzene_rings:
        com_list.append(com(i))   
    nv_host_benzenes = []
    for i,j in zip(benzene_rings,com_list):
        nv_host_benzenes.append(normal_vector(j,i))
    angles = []
    for i in nv_host_benzenes:
        angles.append(nv_angle(nv_guest,i))
    return angles

def vector_analysis(vector, atom_list, atom_coord, atom_vdw, increment):
    """
    First part of this function calculates a set of points on a vector, in specified increments, 
    drawn in XYZ system starting at origin.
    Output: 'pathway' list of XYZ coordinates (floats) [[a1,b1,c1],[a2,b2,c2], ... ,[an,bn,cn]]
    Second part calculates the diameter of a biggest sphere that can be drown on a single point in XYZ space,
    with respect to the atom positions and their vdw radii.
    Output: float
    Third part of this function returns only these vectors, that have assigned sphere of positive radius value
    for the whole vector pathway
    Otherways it returns 'None'
    Output: list [[inc1, s_d1, a1, b1, c1, x, y, z], ... , [incn, s_dn, an, bn, cn, x, y, z]]
    inc1 - the increment on the vector assigning a point that was analysed
    s_d1 - the sphere diameter for this point
    a, b, c - XYZ coordinates of this point
    x, y, z - XYZ coordinates of analysed vector (unchanging value)
    UPDATE 30/11/15: added a treshold for the sphere size to be at least bigger than half minimum molecular
    diameter of hydrogen which is 2.18 so 1.09 -> lets make it 1A in diameter cutoff!
    UPDATE 2: No treshold should be used, not at this step at least, because it is only raw etimation!
    """
    pathway = []
    chunks = int(round((((vector[0])**2 + (vector[1])**2 + (vector[2])**2)**0.5) / increment, 0))
    for i in range(chunks+1):
        pathway.append(np.array([vector[0]/chunks*i, vector[1]/chunks*i, vector[2]/chunks*i]))
    increments = 0
    sphere_sizes = []
    values_mat = []
    sphere_sizes_mat = []
    for i in pathway:
        sphere_diameter_li = []
        dist_mat = euclidean_distances(atom_coord, i)
        dist_mat_red = (dist_mat - atom_vdw) * 2
        sphere_size_mat = np.amin(dist_mat_red)
        if sphere_size_mat > 0:
            values_mat.append([increments, sphere_size_mat, i[0], i[1], i[2], vector[0], vector[1], vector[2]])
            sphere_sizes_mat.append(sphere_size_mat)
        else:
            break
        increments += increment
    if len(values_mat) == len(pathway):
        return(values_mat[min(enumerate(sphere_sizes_mat), key=op.itemgetter(1))[0]])
    else:
        return(None)
    
def no_of_discrete_molecules(atom_list):
    """
    A function that iterates through atom list of type [id, no, x, y, z] and finds discrete molecules
    within single list and spits out list of lists containing each seperate molecule [[mol-1], ..., [mol-n]]
    Withinn a list the order is more or less random. The neighbours distance cutoff is 2 Angstroms and the
    allowed bond length cutoff is 1.6 Angstroms
    EDIT (15-10-15): For C-Br the distance is 1.91, need to be changed for 2 Angstroms
    Output: no of descrete molecules (int), list of lists
    """
    output = []
    alist = copy.deepcopy(atom_list)          #Deep copy of the original atom list not to modify it
    iteration = 0
    while len(alist) > 0:                     #Continue as long as there are unassigned atoms
        final_list = []                       #Final list for a single iteration contaning single molecule
        length_of_alist = len(alist)
        iteration += 1 
        final_list.append(alist[0])           #Append first atom from the pdb list, it's random which one is first
        alist.pop(0)                          #Delete this atom from the pdb list to avoid infinite loops
        for i in final_list:                  #Although it is for every i, it is usually the first one only 
            temp_list = []                    #Temp storage for all bonded neighbours to our i 
            for j in alist:                   #If a neighbour is found and the distance is ok, append it!
                if i[2]-2 < j[2] < i[2]+2 and i[3]-2 < j[3] < i[3]+2 and i[4]-2 < j[4] < i[4]+2:
                    if 0.1 < two_points_distance(i[2:],j[2:]) < 2:
                        temp_list.append(j)
        while len(temp_list) > 0:             #As long as there are candidates we look for neighbours for them!
            temp_list = uniq(temp_list)       #Because of 'cyclic' problems we delete replicates
            final_list.append(temp_list[0])   #Similiar to previous steps
            ref = temp_list.pop(0)         
            alist.remove(ref)                 #Here it works like a snow ball efect, especially with cages
            for j in alist:                   #The search spreads like a web and ends quickly
                if ref[2]-2 < j[2] < ref[2]+2 and ref[3]-2 < j[3] < ref[3]+2 and ref[4]-2 < j[4] < ref[4]+2:
                    if 0.1 < two_points_distance(ref[2:],j[2:]) < 2:
                        temp_list.append(j)
        output.append(final_list)             #Append results for SINGLE molecue to output list as a list
        no_of_atoms = length_of_alist - len(alist) 
    return(iteration,output) 

#Function that translates coordinates from HISTORY file
def coor_translate2(h,c):
    history_list = copy.deepcopy(h)
    zero_trans = []
    x_trans = []
    y_trans = []
    z_trans = []
    xy_trans = []
    xz_trans = []
    yz_trans = []
    xyz_trans = []
    c_a = c[0]
    c_b = c[1]
    c_c = c[2]
    for i in history_list:
        #When in position 1. No translation is necessary
        if 0 < i[2] <= c_a/2 and 0 < i[3] <= c_b/2 and 0 < i[4] <= c_c/2+0.01:
            zero_trans.append(i[1])
        #When in position 2. Need to add c_c
        elif 0 < i[2] <= c_a/2 and 0 < i[3] <= c_b/2 and i[4] <= 0:
            z_trans.append(i[1])
            i[4] = i[4] + c_c
        #When in position 3. Need to add c_a
        elif i[2] <= 0 and 0 < i[3] <= c_b/2 and 0 < i[4] <= c_c/2+0.01:
            x_trans.append(i[1])
            i[2] = i[2] + c_a
        #When in position 4. Need to add c_a and c_c
        elif i[2] <= 0 and 0 < i[3] <= c_b/2 and i[4] <= 0:
            xz_trans.append(i[1])
            i[2] = i[2] + c_a
            i[4] = i[4] + c_c
        #When in position 5. Need to add c_b
        elif 0 < i[2] <= c_a/2 and i[3] <= 0 and 0 < i[4] <= c_c/2+0.01:
            y_trans.append(i[1])
            i[3] = i[3] + c_b
        #When in position 6. Need to add c_b and c_c
        elif 0 < i[2] <= c_a/2 and i[3] <= 0 and i[4] <= 0:
            yz_trans.append(i[1])
            i[3] = i[3] + c_b
            i[4] = i[4] + c_c
        #When in position 7. Need to add c_a and c_b
        elif i[2] <= 0 and i[3] <= 0 and 0 < i[4] <= c_c/2+0.01:
            xy_trans.append(i[1])
            i[2] = i[2] + c_a
            i[3] = i[3] + c_b
        #When in position 8. Need to add c_a and c_b and c_c
        elif i[2] <= 0 and i[3] <= 0 and i[4] <= 0:
            xyz_trans.append(i[1])
            i[2] = i[2] + c_a
            i[3] = i[3] + c_b
            i[4] = i[4] + c_c
        else:
            print(i[1], i[2],i[3],i[4], c_a,c_b,c_c, ' couldnt allocate')
    return(history_list)

def coor_translate(history_coor,c,progress):
    """
    This is possibly the most important function of this script. It takes the atom coordinates
    within the unit cell and reconstructs (rebuilds) the cages into whole molecules
    """
    sys.stdout.flush()
    h = copy.deepcopy(history_coor)
    c_a = c[0]
    c_b = c[1]
    c_c = c[2]
    cell_pseudocenter = [c_a/4,c_b/4,c_c/4]
    final_list = []
    list_of_number = []
    list_of_number_temp = []
    first_length = 0
    second_length = 0
    itter = 0
    num_of_atoms = []
    num_of_atoms_temp = []
    list_of_lengths = []
    list_to_append = []
    number_of_atoms = 0
    inside_atoms = copy.deepcopy(h)
    halo_atoms = halo(inside_atoms,c)
    while second_length >= first_length:
        working_list_temp = []
        temp_list_to_append = []
        if second_length == first_length:
            list_of_number_temp = uniq(list_of_number_temp)
            for i in list_of_number_temp:
                list_of_number.append(i)
            list_of_number_temp = []
        if second_length == first_length and len(list_of_number) != len(h):
            if len(list_of_number) > 0:
                num_of_atoms.append(len(list_of_number)-num_of_atoms_temp[-1])
                for i in list_of_number[number_of_atoms:len(list_of_number)]:
                    temp_list_to_append.append(i)
                number_of_atoms = number_of_atoms + len(list_of_number)-num_of_atoms_temp[-1]
                list_of_lengths[-1].append(num_of_atoms[-1])
                list_to_append.append(temp_list_to_append)
            num_of_atoms_temp.append(len(list_of_number))
            itter += 1
            prog = round(len(list_of_number)*100/len(h),0) + 2
            sys.stdout.write("\rOverall progress: %d%% | Frame: 1 | Wrapping molecules progress: %d%%" % (progress,prog-1))
            sys.stdout.flush()
            new_inside_atoms = []
            for i in inside_atoms:
                if i[1] not in list_of_number:
                    new_inside_atoms.append(i)
            inside_atoms = new_inside_atoms
            working_list = []
            starting_point_temp = []
            for i in inside_atoms:
                if i[0] == 'CA':
                    starting_point_temp.append((i,two_points_distance(i[2:],cell_pseudocenter)))
            starting_point_temp = sorted(starting_point_temp, key=op.itemgetter(1))
            new_start = starting_point_temp[0][0]
            working_list.append(starting_point_temp[0][0])
            list_of_lengths.append([new_start[1]])
            inside_atom_set = []
            for i in inside_atoms:
                if two_points_distance(i[2:5],new_start[2:5]) < 16:
                    inside_atom_set.append(i)
            halo_atom_set = []
            for i in halo_atoms:
                if two_points_distance(i[2:5],new_start[2:5]) < 16:
                    halo_atom_set.append(i)
        if len(list_of_number) == len(h):
            num_of_atoms.append(len(list_of_number)-num_of_atoms_temp[-1])
            for i in list_of_number[number_of_atoms:len(list_of_number)]:
                temp_list_to_append.append(i)
            number_of_atoms = number_of_atoms + len(list_of_number)-num_of_atoms_temp[-1]
            list_to_append.append(temp_list_to_append)
            list_of_lengths[-1].append(num_of_atoms[-1])
            break
        first_length = len(working_list_temp)
        working_list = uniq(working_list)
        for i in working_list:
            if i[1] not in list_of_number_temp:
                if i[0] != 'HA' and i[0] != 'HC' and i[0] != 'HE':
                    for j in inside_atom_set:
                        if 0.1 < two_points_distance(i[2:5],j[2:5]) < 1.8:
                            working_list_temp.append(j)
                    if i[2] < 1.8 or i[3] < 1.8 or i[4] < 1.8 or i[2] > c_a-1.8 or i[3] > c_b-1.8 or i[4] > c_c-1.8:
                        for j in halo_atom_set:
                            if 0.1 < two_points_distance(i[2:5],j[2:5]) < 1.8:
                                working_list_temp.append(j)
        second_length = len(working_list_temp)
        for i in working_list:
            final_list.append(i)
            list_of_number_temp.append(i[1])
        working_list_temp = uniq(working_list_temp)
        for i in working_list_temp:
            working_list.append(i)

    final_list = uniq(final_list)
    atom_sequence = 0
    for i in final_list:
        atom_sequence += 1
        i.append(atom_sequence)
    sorted_working_list = sorted(final_list, key=op.itemgetter(1))
    num_of_descrete = len(num_of_atoms)
    num_of_distinct = len(unique(num_of_atoms))
    cell_description = []
    cell_description.append(num_of_descrete)
    cell_description.append(num_of_distinct)
    #list_of_lengths = sorted(list_of_lengths, key=itemgetter(0))
    list_of_len = []
    for i in list_of_lengths:
        list_of_len.append(i[1])
    cell_description.append(list_of_len)
    cell_description.append(list_to_append)
    return(sorted_working_list,cell_description)

#Optimisation with scipy over x,y plane, the returning value is negative of radius
def analyse_xy(position, *params):
    x_pos, y_pos = position
    z_pos, atom_list_trans = params
    return(-min([two_points_distance([x_pos,y_pos,z_pos],i[2:5])-atom_vdw_radii[i[0]]                                                        for i in atom_list_trans]))
def analyse_z(position, *params):
    z_pos = position
    x_pos, y_pos, atom_list_trans = params
    return(min([two_points_distance([x_pos,y_pos,z_pos],i[2:5])-atom_vdw_radii[i[0]]                                                        for i in atom_list_trans]))

def window_analysis(window, atom_list, atom_coor, atom_vdw):
    max_value = max(np.array(window)[:,1])
    number = [a for a, j in enumerate(window) if j[1] == max_value]
    vector2analyse = window[number[0]][5:8]
    vector_analysed = vector_analysis(vector2analyse, atom_list, atom_coor, atom_vdw, 0.1) 

    #UPDATE: Try rotation first, than translation
    vector_main = np.array([vector_analysed[5],vector_analysed[6],vector_analysed[7]])
    vec_a = [1, 0, 0]
    vec_b = [0, 1, 0]
    vec_c = [0, 0, 1]
    angle_2 = vec_angle(vector_main, vec_c)
    angle_1 = vec_angle([vector_main[0], vector_main[1], 0], vec_a)

    if vector_main[0] >= 0 and vector_main[1] >= 0 and vector_main[2] >= 0:
        angle_1 = -angle_1
        angle_2 = -angle_2                
    if vector_main[0] < 0 and vector_main[1] >= 0 and vector_main[2] >= 0:
        angle_1 = np.pi*2 + angle_1
        angle_2 = angle_2
    if vector_main[0] >= 0 and vector_main[1] < 0 and vector_main[2] >= 0:
        angle_1 = angle_1
        angle_2 = -angle_2
    if vector_main[0] < 0 and vector_main[1] < 0 and vector_main[2] >= 0:
        angle_1 = np.pi*2 -angle_1
    if vector_main[0] >= 0 and vector_main[1] >= 0 and vector_main[2] < 0:
        angle_1 = -angle_1
        angle_2 = np.pi + angle_2
    if vector_main[0] < 0 and vector_main[1] >= 0 and vector_main[2] < 0:
        angle_2 = np.pi - angle_2
    if vector_main[0] >= 0 and vector_main[1] < 0 and vector_main[2] < 0:
        angle_2 = angle_2 + np.pi
    if vector_main[0] < 0 and vector_main[1] < 0 and vector_main[2] < 0:
        angle_1 = -angle_1
        angle_2 = np.pi - angle_2 

    #First rotation around z-axis with angle_1

    rot_matrix_z = np.array([[np.cos(angle_1), -np.sin(angle_1),      0], 
                             [np.sin(angle_1),  np.cos(angle_1),      0],
                             [                0,                  0,      1]])

    resulting_vector = np.dot(rot_matrix_z, vector_main)
    atom_list_translated = [[i[0],i[1],np.dot(rot_matrix_z, i[2:])[0],
                            np.dot(rot_matrix_z, i[2:])[1],
                            np.dot(rot_matrix_z, i[2:])[2]] for i in atom_list]

    #Second rotation around y-axis with angle_2
    rot_matrix_y = np.array([[ np.cos(angle_2),         0,       np.sin(angle_2)],
                             [               0,          1,                     0],
                             [-np.sin(angle_2),         0,       np.cos(angle_2)]])

    resulting_vector2 = np.dot(rot_matrix_y, resulting_vector)
    atom_list_translated = [[i[0],i[1],np.dot(rot_matrix_y, i[2:])[0],
                            np.dot(rot_matrix_y, i[2:])[1],
                            np.dot(rot_matrix_y, i[2:])[2]] for i in atom_list_translated]

    #Third step is translation! We are now at approximetely [0,0,-z] 
    #We need to shift the origin into the point of the window
    #the value for z we know from the original vector analysis (it is the length on vector where
    #there was the biggest sphere (step - vector_analysed[0]) first value!)
    #We can reason that because you can see that the length of original vector and now the
    #z value for rotated point is the same! the length of vector is preserved
    new_z = vector_analysed[0]
    old_origin = np.array([0,0,-new_z])
    new_resulting_vector3 = np.add(resulting_vector2, old_origin)
    atom_list_translated2 = [[i[0],i[1],i[2],i[3],i[4]-new_z] for i in atom_list_translated]

    #!!!Here the xy and z sampling has to take place!!!
    #First sample the same point to check if nothing has changed
    #This point should represent 0,0,0 the origin
    
    distance_list = []
    for i in atom_list_translated2:
        distance_list.append(np.linalg.norm(np.array([i[2:5]]))-atom_vdw_radii[i[0]])

    iteration_opt = 0
    stdev_1 = 0
    stdev_2 = 0

    x_opt = 0
    y_opt = 0
    z_opt = new_z

    xyz_window = []

    xyz_window.append(2*min([two_points_distance([x_opt,y_opt, z_opt],i[2:5])-atom_vdw_radii[i[0]]                                                    for i in atom_list_translated]))

    parameters1 = (x_opt, y_opt, atom_list_translated)
    normal_optimisation = scipy.optimize.minimize(analyse_z, x0=z_opt, args=parameters1)
    z_opt = normal_optimisation.x[0]

    #BRUTE
    rranges = ((-max_value/2, max_value/2), (-max_value/2, max_value/2))
    parameters2 = (z_opt, atom_list_translated)
    brute_optimisation = scipy.optimize.brute(analyse_xy, rranges, args=parameters2, 
                                             full_output=True, finish=scipy.optimize.fmin)

    x_opt = brute_optimisation[0][0]
    y_opt = brute_optimisation[0][1]
    xyz_window.append(2*min([two_points_distance([x_opt,y_opt, z_opt],i[2:5])-atom_vdw_radii[i[0]]                                                for i in atom_list_translated]))
    
    #control_point_1 = np.array([x_opt, y_opt, z_opt])
    
    #Reverse translation step
    #rev_resulting_vector3 = np.add(new_resulting_vector3, [0,0,new_z])
    
    #Reversing the second rotation around axis y
    #angle_2_1 = - angle_2 
    #rev_matrix_y = np.array([[ np.cos(angle_2_1),         0,       np.sin(angle_2_1)],
    #                         [               0,          1,                     0],
    #                         [-np.sin(angle_2_1),         0,       np.cos(angle_2_1)]])

    #control_point_1 = np.dot(rev_matrix_y, control_point_1)
    #rev_resulting_vector4 = np.dot(rev_matrix_y, rev_resulting_vector3)

    #Reversing the first rotation around axis z
    #angle_1_1 = - angle_1 
    #rev_matrix_z = np.array([[np.cos(angle_1_1), -np.sin(angle_1_1),      0], 
    #                         [np.sin(angle_1_1),  np.cos(angle_1_1),      0],
    #                         [                0,                  0,      1]])

    #control_point_1 = np.dot(rev_matrix_z, control_point_1)
    #rev_resulting_vector5 = np.dot(rev_matrix_z, rev_resulting_vector4)
    return(xyz_window[1])

class XYZ:
    """
    This class load and open an XYZ file.
    self.path: Path of the file
    self.name: Name of the file
    self.body: list of string as in loaded file
    self.no_of_atoms: total number of atoms
    self.comment: Second (comment) line in XYZ file
    self.atom_list: Creates list of kind [[id1, no1, x1, y1, z1], ...,[idn, non, xn, yn, zn]]
    self.array: It also creates a numpy array containing only atom coordinates [[x1,y1,z1], ..., [xn,yn,zn]]
    self.com2zero(): Translates the atom coordinates so that COM is at origin
    """
    def __init__(self, xyz):
        self.path = xyz
        self.name = xyz.split('/')[-1][:-4]
        with open(xyz, 'r') as xyz_source:
            self.body = [i.split() for i in xyz_source.readlines()]
            self.no_of_atoms = int(self.body[0][0])
            self.comment = self.body[1][0]
        iteration = 1
        atom_list = []
        array = []
        for i in self.body[2:]:
            atom_list.append([i[0].upper(), iteration, float(i[1]), float(i[2]), float(i[3])])
            array.append([float(i[1]), float(i[2]), float(i[3])])
            iteration += 1
        self.atom_list = atom_list
        self.array = np.array(array)
    
    def com2zero(self):
        """
        This function first calculates the center of mass of the molecule and than translate all the coordinates
        So that new center of mass ends up in origin.
        It is done for both self.atom_list and self.xyz_array instances of this class
        """
        com = center_of_mass(self.atom_list)
        self.atom_list = [[i[0],i[1],i[2]-com[0],i[3]-com[1],i[4]-com[2]] for i in self.atom_list]
        self.array = np.array([[i[0]-com[0],i[1]-com[1],i[2]-com[2]] for i in self.array])
        
class MOL:
    """
    This class load and open an MOL file.
    self.path: Path of the file
    self.name: Name of the file
    self.body: list of string as in loaded file
    self.no_of_atoms: total number of atoms
    self.comment: Second (comment) line in XYZ file
    self.atom_list: Creates list of kind [[id1, no1, x1, y1, z1], ...,[idn, non, xn, yn, zn]]
    self.array: It also creates a numpy array containing only atom coordinates [[x1,y1,z1], ..., [xn,yn,zn]]
    """
    def __init__(self, mol):
        self.path = mol
        self.name = mol.split('/')[-1][:-4]
        with open(mol, 'r') as mol_source:
            self.body = [i.split() for i in mol_source.readlines()]
            #self.no_of_atoms = int(self.body[0][0])
            #self.comment = self.body[1][0]
        iteration = 1
        atom_list = []
        array = []
        flag = False
        for i in self.body:
            if len(i) >= 4:
                if i[2] == 'END' and i[3] == 'ATOM':
                    flag = False
                if flag == True:
                    atom_list.append([i[3].upper(), iteration, float(i[4]), float(i[5]), float(i[6])])
                    array.append([float(i[4]), float(i[5]), float(i[6])])
                    iteration += 1
                if i[2] == 'BEGIN' and i[3] == 'ATOM':
                    flag = True
        self.atom_list = atom_list
        self.array = np.array(array)
        
    def com2zero(self):
        """
        This function first calculates the center of mass of the molecule and than translate all the coordinates
        So that new center of mass ends up in origin.
        It is done for both self.atom_list and self.xyz_array instances of this class
        """
        com = center_of_mass(self.atom_list)
        self.atom_list = [[i[0],i[1],i[2]-com[0],i[3]-com[1],i[4]-com[2]] for i in self.atom_list]
        self.array = np.array([[i[0]-com[0],i[1]-com[1],i[2]-com[2]] for i in self.array])
        
class PDB:
    """
    This class load and open an PDB file.
    self.path: Path of the file
    self.name: Name of the file
    self.block: list of string as in loaded file
    self.remarks: All the REMARK lines in PDB
    self.crystal: Crystal structure parameters [a, b, c, alpha, beta, gamma]
    self.body: Properly extracted PDB data (see the PDB v3.3 manual)
    self.connect: Connectivity block in PDB file
    self.atom_list: Creates list of kind [[id1, no1, x1, y1, z1], ...,[idn, non, xn, yn, zn]]
    self.array: It also creates a numpy array containing only atom coordinates [[x1,y1,z1], ..., [xn,yn,zn]]
    """
    def __init__(self, pdb):
        self.path = pdb
        self.name = pdb.split('/')[-1][:-4]
        self.block = [i.strip('\n') for i in open(pdb, 'r').readlines()]
        self.remarks = [i for i in self.block if i[:6] == 'REMARK']
        self.crystal = [float(x) for i in self.block for x in [i[6:15], i[15:24], i[24:33], i[33:40],
                         i[40:47], i[47:54]] if i[:6] == 'CRYST1']
        self.body = [[i[0:6], i[6:11], i[12:16], i[16], i[17:20], i[21], i[22:26], i[26],
                     i[30:38], i[38:46], i[46:54], i[54:60], i[60:66], i[76:78], i[78:80]] \
                    for i in self.block if i[:6] == 'HETATM' or i[:6] == 'ATOM  ']
        self.conect = [i for i in self.block if i[:6] == 'CONECT']
        
        iteration = 1
        atom_list = []
        array = []
        for i in self.body:
            atom_list.append([i[-2].strip().upper(), iteration, float(i[8]), float(i[9]), float(i[10])])
            array.append([float(i[8]), float(i[9]), float(i[10])])
            iteration += 1
        self.atom_list = atom_list
        self.array = np.array(array)
        
    def com2zero(self):
        """
        This function first calculates the center of mass of the molecule and than translate all the coordinates
        So that new center of mass ends up in origin.
        It is done for both self.atom_list and self.xyz_array instances of this class
        """
        com = center_of_mass(self.atom_list)
        self.atom_list = [[i[0],i[1],i[2]-com[0],i[3]-com[1],i[4]-com[2]] for i in self.atom_list]
        self.array = np.array([[i[0]-com[0],i[1]-com[1],i[2]-com[2]] for i in self.array])


# In[35]:

#@profile
def sub_main(file_atom_list, verbose=False, figures=True, adjust=250, psd_output=False):
    """
    This function is responsible for the main analysis. It requires an atom list for a single molecule, crystal
    or as in case of trajectory file - single frame.
    If verbose = True: print output in the terminal
    scale: by each the number of samples will change with each iteration
    psd_output: if the txt file ready to plot the PSD should be saved at the end
    """
    #current_milli_time = int(round(time.time() * 1000))
    #print('start mili function :', current_milli_time)
    output_list = []
    discrete = no_of_discrete_molecules(file_atom_list) #Each discrete molecule has to be analysed seperately 
    if verbose == True:
        print("Number of discrete molecules found: {0}".format(discrete[0]))
    iteration = 0
    data = []
    for i in discrete[1]:  #For each found discrete molecule
        iteration += 1
        atom_list = i
        com = center_of_mass(atom_list)           #Center of mass
        cog = center_of_geometry(atom_list)       #Center of geometry (mass distribution can be uneven)
        atom_list_old = atom_list
        atom_list = com2zero(atom_list)           #Shift COM to origin for the sake of simplicity
        atom_coor = np.array([[x[2],x[3],x[4]] for x in atom_list])
        atom_vdw = np.array([[atom_vdw_radii[x[0]]] for x in atom_list])
        atom_vdw_vertical = np.matrix([[atom_vdw_radii[x[0]]] for x in atom_list])
        atom_vdw_horizontal = np.matrix([atom_vdw_radii[x[0]] for x in atom_list])
        com_zero = center_of_mass(atom_list)      #Double check that it is shifted to the origin
        maximum_dimension = max_dim(atom_list, atom_coor, atom_vdw_vertical, atom_vdw_horizontal)
        void_diam = void_diameter(atom_list, atom_coor, atom_vdw)
        output_list.append(void_diam)
        output_list.append(maximum_dimension)
        if verbose == True:
            print("\n    Analysing molecule {0}".format(iteration))
            print("    COM: {0:.2f} {1:.2f} {2:.2f}".format(*com))
            print("    COG: {0:.2f} {1:.2f} {2:.2f}".format(*cog))
            print("    COM_0: {0:.2f} {1:.2f} {2:.2f}".format(*com_zero))
            print("    Maximum dimension {0:.2f}".format(maximum_dimension))
            print("    Void diameter {0:.2f}".format(void_diam))
            print("    Kr volume % of void volume: {0:.2f}".format(3.69/void_diam*100))
            print("    Xe volume % of void volume: {0:.2f}".format(4.10/void_diam*100))

        output_list2 = []
        window_coms = []
        R = maximum_dimension/2 #Finding out the sphere radius as a half of maximum cage dimension
        sphere_surface_area = 4 * np.pi * R * R #Calculating surface area of a halo sphere
        number_of_points = int(round(np.log10(sphere_surface_area)*adjust, 0)) #quadratic_to_logaritmic
        output_list2.append(number_of_points)
        points_per_1A = number_of_points/sphere_surface_area
        points_per_5A_round = int(round(points_per_1A*5, 0))
        if verbose == True:
            print("    Number of sampling points: {0}".format(number_of_points))
            print("    Number of sampling points per A^2: {0:.2f}".format(points_per_1A))
            print("    Number of sampling points per 5A^2: {0:.2f}".format(points_per_5A_round))

        #Estimating golden angle between the points on the sphere and creating a halo points sphere
        golden_angle = np.pi * (3 - np.sqrt(5))
        theta = golden_angle * np.arange(number_of_points)
        z = np.linspace(1 - 1.0 / number_of_points, 1.0 / number_of_points - 1.0, number_of_points)
        radius = np.sqrt(1 - z * z)

        points = np.zeros((number_of_points, 3))
        points[:,0] = radius * np.cos(theta) * R
        points[:,1] = radius * np.sin(theta) * R
        points[:,2] = z * R

#        print('\nSKLERN MATRIX')
        values = []
        tree = KDTree(points)
        for i in points:
            dist, ind = tree.query(i, k=10)
            values.append(dist[0][1])
            values.append(dist[0][2])
            values.append(dist[0][3])

        eps = np.mean(values) 
        eps_sqrt = eps + eps**0.5
#        print('EPS_KDTree_sqrt ', eps_sqrt)

        #Multiprocessing    
#        pool = mp.Pool(processes=8)
#        parallel = [pool.apply_async(vector_analysis, 
#                    args=(x, atom_list, atom_coor, atom_vdw, 2.0)) for x in points]
#        results = [p.get() for p in parallel if p.get() is not None]
#        pool.terminate()
#        dataset = [x[5:8] for x in results]
#        all_basins = [x[1] for x in results]

        #Use this part instead if you don't want parallel
        results = []
        for x in points:
            results.append((vector_analysis(x, atom_list, atom_coor, atom_vdw, 2.0)))
        results = [x for x in results if x is not None]
        dataset = [x[5:8] for x in results if x is not None]
        all_basins = [x[1] for x in results if x is not None]

        if len(results) == 0:
            return(None)
#            print('No windows were found. Finishing the analysis')
        else:
            biggest_window = max(all_basins)
            output_list2.append(biggest_window)

            X = np.array(dataset)

            if points_per_5A_round/2 > 5:
                points_per_5A_round = 10
#            print("    Number of sampling points per 5A^2 used: {0:.2f}".format(points_per_5A_round))
            # Compute DBSCAN
            db = DBSCAN(eps=eps_sqrt).fit(X)
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_

            # Number of clusters in labels, ignoring noise if present.
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            if verbose == True:
                print("    Estimated number of clusters: {0}".format( n_clusters_))
            data.append(n_clusters_)

            clusters = []
            for i,j in zip(X, db.labels_):
                clusters.append([i,j])

            clustered_results = {}
            for i in range(n_clusters_):
                no_of_samples = 0
                clustered_results[i] = []
                for j,k in zip(clusters,results):
                    if j[1] == i:
                        clustered_results[i].append(k)
            
            #Here starts the multiprocessing for the window size analysis and minimasation
#            pool = mp.Pool(processes=8)
#            window_parallel = [pool.apply_async(window_analysis, args=(clustered_results[x], atom_list, atom_coor, 
#                                        atom_vdw,)) for x in clustered_results]
#            window_results = [p.get() for p in window_parallel if p.get() is not None]
#            pool.terminate()
            
            #Use this part instead for serial analysis
            window_results = []
            for i in clustered_results:
                window_results.append(window_analysis(clustered_results[i], atom_list, atom_coor, atom_vdw))
                
            #for i in window_results:
            #    print(i)
            #    print(np.add(i[1], com))

    #current_milli_time1 = int(round(time.time() * 1000))
    #print('end mili function :', current_milli_time1- current_milli_time)
    return(window_results)
        

def main():
    """
    This function handles the recognition of file types (XYZ,PDB,MOL and 'HISTORY' (trajectory files) are allowed)
    It than passes the data to the sub_main() function that is responsible for the main analysis
    HISTORY trajectory files as much more complicated requires a lot more steps of loding it,
    extracting each frame and than passing this frame to the sub_main() function. Therefore much of this
    operation has to be done here, as I want sub_main() function to work with each of the allowed inputs.
    """
    print('Output generated on {0} {1}\n'.format(time.strftime("%H:%M:%S"), time.strftime("%d/%m/%Y")))
    
    for i in glob.glob('molecules/*'):   #ENRICO! Create an input dir in your working dir and run the script :)
        if i.split('/')[-1][-3:] == 'xyz':               #If it is XYZ file start the analysis rout for XYZ files
            file = XYZ(i)
            print('\n                    Analysing {0}'.format(file.path))
            print('Start time {0}\n'.format(time.strftime("%H:%M:%S")))
            file_atom_list = file.atom_list
            sub_main(file_atom_list) 
        elif i.split('/')[-1][-3:] == 'pdb':             #If it is PDB file start the analysis rout for PDB files
            file = PDB(i)
            print('\n                    Analysing {0}'.format(file.path))
            print('Start time {0}\n'.format(time.strftime("%H:%M:%S")))
            file_atom_list = file.atom_list
            sub_main(file_atom_list)
        elif i.split('/')[-1][-3:] == 'mol':             #If it is MOL file start the analysis rout for MOL files
            file = MOL(i)
#            print('\n                    Analysing {0}'.format(file.path))
#            print('Start time {0}\n'.format(time.strftime("%H:%M:%S")))
            file_atom_list = file.atom_list
            print(sub_main(file_atom_list))
        elif i.split('/')[-1].split('_')[0] == 'HISTORY' or i.split('/')[-1].split('-')[0] == 'HISTORY':
            #History file need special means of loding it, than extracting each frame and each frame
            #Can than be analysed with the sub_main() function as examples above
            timestep_flag = False
            reconstruction_flag = False
            end_flag = False
            global_line = 0
            frame = 0
            frame_line = 0
            binary_step = 0
            progress = 0
            frame_data = []
            frame_data_line = []
            cell_parameters = []
            
            with open(i, 'r') as history_file:
                with ct.closing(mmap.mmap(history_file.fileno(), 0, access=mmap.ACCESS_READ)) as history_binary_map:
                    while binary_step <= len(history_binary_map):
                        global_line += 1
                        binary_line = history_binary_map.readline()
                        if len(binary_line) == 0:
                            end_flag = True
                        binary_step = binary_step + len(binary_line)
                        progress_old = progress
                        progress = round(binary_step*100/len(history_binary_map),0)
                        if progress_old < progress:
                            sys.stdout.write("\rOverall progress: {:d}% | Frame: {:d}".format(int(progress),frame))
                            sys.stdout.flush()
                        string_line = binary_line.decode("utf-8").strip('\n').split()
                        if global_line == 2:
                            periodic_boundry_key = int(string_line[1])
                            # 0 is for non-periodic 1 is for cubic 2 is for orthorombic 3 is for parallelepiped
                            # 4-7 too complicated not covered here
                            if periodic_boundry_key == 0:
                                reconstruction_flag = False
                            elif periodic_boundry_key == 1 or periodic_boundry_key == 2 or periodic_boundry_key == 3:
                                reconstruction_flag = True
                            else:
                                print("The periodic boundry condition of this trajectory file is not covered by " + 
                                     "this script. This file will not be analysed. \nTerminate\n")
                                break
                            number_of_atoms = int(string_line[2])
                        #THSI IS TO BE DONE!!!!!!!!!!!!!!!!!
                        """if timestep_flag == True and reconstruction_flag == True:
                            frame_line += 1
                            if frame_line == 1:
                                cell_parameters.append(float(string_line[0]))
                                cell_matrix.append([float(string_line[0]),float(string_line[1]),float(string_line[2])])
                            if frame_line == 2:
                                cell_parameters.append(float(string_line[1]))
                                cell_matrix.append([float(string_line[0]),float(string_line[1]),float(string_line[2])])
                            if frame_line == 3:
                                cell_parameters.append(float(string_line[2]))
                                cell_matrix.append([float(string_line[0]),float(string_line[1]),float(string_line[2])])
                                cell_matrix = np.array(cell_matrix)
                            if frame_line > 3 and frame_line % 2 == 0:
                                #Need to translate the OPLS atom types here!!!
                                if string_line[0].upper() in opls_atom_keys['C']:
                                    frame_data_line.append('C')
                                elif string_line[0].upper() in opls_atom_keys['H']:
                                    frame_data_line.append('H')
                                elif string_line[0].upper() in opls_atom_keys['O']:
                                    frame_data_line.append('O')
                                elif string_line[0].upper() in opls_atom_keys['S']:
                                    frame_data_line.append('S')
                                elif string_line[0].upper() in opls_atom_keys['N']:
                                    frame_data_line.append('N')
                                elif string_line[0].upper() in opls_atom_keys['F']:
                                    frame_data_line.append('F')
                                else:
                                    if string_line[0].upper() != 'TIMESTEP':
                                        print('Atom key {0} has not been found in the '.format(string_line[0].upper())
                                             + 'OPLS atom keys library. I have left it unchanged, although this might '
                                             + 'cause the script to breake')
                                frame_data_line.append(int(string_line[1]))
                            if frame_line > 3 and frame_line % 2 != 0:
                                frame_data_line.append(float(string_line[0]))
                                frame_data_line.append(float(string_line[1]))
                                frame_data_line.append(float(string_line[2]))
                                frame_data.append(frame_data_line)
                                frame_data_line = []
                            if string_line[0] == 'timestep':
                                print('UNIT CELL')
                                print('START')
                                if frame == 1:
                                    if number_of_atoms != len(frame_data):
                                        print("Number of atoms does not equal number of entries in frame data")
                                        break
                                    # Here we extract a trajectory frame 
                                    host_xyz_original = copy.deepcopy(frame_data)
                                    
                                    # We reconstruct the solid-state structure from the trajectory file for the host
                                    # First we translate atom positions from a 0,0,0 centered unit cell to a normalised
                                    # unit cell where 0,0,0 is the beginning of the unit cell and c_a/2,c_b/2,c_c/2 
                                    host_xyz_step_1 = coor_translate2(host_xyz_original, cell_parameters)

                                    # This is the most important part. Now the cage molecules are reconstructed
                                    # So that each molecule is a whole. For this to happen, they has to be allowed
                                    # to penetrate beyond unit cell boundaries. Number of molecules per cage
                                    # cell paramaeters with numbers assigned to molecules are also produced
                                    # However!!! The resulting list is sorted with the original atom numbers, not the good ones
                                    # Therefore we sort it later once again. But to create matrix we need the original numbering
                                    host_xyz_step_2 = coor_translate(host_xyz_step_1, cell_parameters, progress)
                                    frame_to_pdb(host_xyz_step_2[0], frame, cell_parameters, history_file_name, 'step3')

                                    # Some cage molecules might be in a wrong positions in respect to the unit cell.
                                    # The idea is that for the rebuild unit cell to have everything in good positions
                                    # A COM of a cage molecule has to be as close as possible to 1/4,1/4,1/4 of the unit cell
                                    host_xyz_step_3 = position_swap(host_xyz_step_2[0], cell_parameters, host_xyz_step_2[1])
                                    frame_to_pdb(host_xyz_step_3, frame, cell_parameters, history_file_name, 'step4')

                                    # We create translational matrix used for the next frame it contains:
                                    # Atom id (neccessery so that we do not assign HE every time)
                                    # Original atom no. (it has to be sorted this way so that we can do list comprehension between
                                    # and the next frames. In other way it takes a really long time to itterate)
                                    # Good atom no. (the atom no. that was assigned during the cage reconstruction, we will use it
                                    # as an additional reference for atoms and sort it this way each time before triads recognition)
                                    # A/B/C translational operations for the xyz
                                    matrix = translational_matrix(host_xyz_original, host_xyz_step_4, progress)

                                    # We sort the list with coordinates so that it is in the good order assigned during cage reconstruction
                                    host_xyz_step_5 = sorted(host_xyz_step_4, key=itemgetter(5))
                                    print(host_xyz_step_5[:5])

                                if frame > 1:
                                    # Deepcopy of a reference matrix from the 1st frame
                                    reference_matrix = copy.deepcopy(matrix)

                                    # Deepcopy of a host atoms positions from the trajectory frame no. x
                                    host_xyz = copy.deepcopy(frame_data)

                                    # Updating the matrix, because some of the atoms could change position over periodic boundry
                                    # Therefore we have to compare the positions of atoms from the 1st frame from unchanged trajectory
                                    # and for the unchanged trajectory from this frame and make corrections to the matrix
                                    new_matrix = structure_check(host_xyz,host_xyz_original,cell_parameters,reference_matrix)

                                    # We update the coordinates according to the matrix. This will automatically cover the step for
                                    # the cage reconstruction which is almost the longest. It will also change atom is from HC to HE
                                    # when neccessery and will append the 'good' itteration number into the list this will also help 
                                    # in the next steps
                                    host_xyz_updated = update_coor(host_xyz,new_matrix,cell_parameters)


                                    # We sort according to the second 'good' number position on the list allowing to calculate triads etc.
                                    host_xyz_sorted = sorted(host_xyz_updated, key=itemgetter(5))

                                timestep_condition = False
                                frame_line = 0"""
                                
                        if timestep_flag == True and reconstruction_flag == False:
                            frame_line += 1
                            if frame_line > 0 and frame_line % 2 != 0 and len(string_line) > 0 and end_flag == False:
                                #Need to translate the OPLS atom types here!!!
                                if string_line[0].upper() in opls_atom_keys['C']:
                                    frame_data_line.append('C')
                                elif string_line[0].upper() in opls_atom_keys['H']:
                                    frame_data_line.append('H')
                                elif string_line[0].upper() in opls_atom_keys['O']:
                                    frame_data_line.append('O')
                                elif string_line[0].upper() in opls_atom_keys['S']:
                                    frame_data_line.append('S')
                                elif string_line[0].upper() in opls_atom_keys['N']:
                                    frame_data_line.append('N')
                                elif string_line[0].upper() in opls_atom_keys['F']:
                                    frame_data_line.append('F')
                                else:
                                    if string_line[0].upper() != 'TIMESTEP':
                                        print('Atom key {0} has not been found in the '.format(string_line[0].upper())
                                             + 'OPLS atom keys library. I have left it unchanged, although this might '
                                             + 'cause the script to breake')
                                    frame_data_line.append(str(string_line[0]).upper())
                                frame_data_line.append(int(string_line[1]))
                            if frame_line > 0 and frame_line % 2 == 0 and len(string_line) > 0 and end_flag == False:
                                frame_data_line.append(float(string_line[0]))
                                frame_data_line.append(float(string_line[1]))
                                frame_data_line.append(float(string_line[2]))
                                frame_data.append(frame_data_line)
                                frame_data_line = []
                            
                            if len(binary_line) > 0:
                                if string_line[0] == 'timestep':
                                    if frame >= 1:
                                        if number_of_atoms != len(frame_data):
                                            print("Number of atoms does not equal number of entries in frame data")
                                            break
                                        # Here we extract a trajectory frame 
                                        host_xyz_original = copy.deepcopy(frame_data)
                                        output = sub_main(host_xyz_original)
                                        #if frame == 1:
                                            #line_to_append = ['Frame','Iteration','Void_d','Max_dim','Num_points',
                                            #                 'Biggest_W','Num_of_W']
                                            #for i in range(int((len(output[2])-3)/4)):
                                            #    line_to_append.append('WD{0}'.format(i+1))
                                            #for i in range(int((len(output[2])-3)/4)):
                                            #    line_to_append.append("WD{0}'".format(i+1))
                                            #for i in range(int((len(output[2])-3)/4)):
                                            #    line_to_append.append('WD{0}_COM'.format(i+1))
                                            #for i in range(int((len(output[2])-3)/4)):
                                            #    line_to_append.append("WD{0}'_COM".format(i+1))
                                            #txt_output.append(line_to_append)
                                        second_line_to_append = []
                                        itter = 0
                                        for i in output[2:]:
                                            itter += 1
                                            second_line_to_append.append(frame)
                                            second_line_to_append.append(itter)
                                            second_line_to_append.append(output[0])
                                            second_line_to_append.append(output[1])
                                            for j in i:
                                                if type(j) != list:
                                                    second_line_to_append.append(round(float(j), 4))
                                                if type(j) == list:
                                                    for k in j:
                                                        second_line_to_append.append(round(float(k), 4))
                                            txt_output.append(second_line_to_append)
                                            second_line_to_append = []
                                        

                                    timestep_flag = False
                                    frame_line = 0
                                    
                            if len(binary_line) == 0 and end_flag == True:
                                if number_of_atoms != len(frame_data):
                                    print("Number of atoms does not equal number of entries in frame data")
                                    break
                                # Here we extract a trajectory frame 
                                host_xyz_original = copy.deepcopy(frame_data)
                                output = sub_main(host_xyz_original)
                                second_line_to_append = []
                                itter = 0
                                for i in output[2:]:
                                    itter += 1
                                    second_line_to_append.append(frame)
                                    second_line_to_append.append(itter)
                                    second_line_to_append.append(output[0])
                                    second_line_to_append.append(output[1])
                                    for j in i:
                                        if type(j) != list:
                                            second_line_to_append.append(round(float(j), 4))
                                        if type(j) == list:
                                            for k in j:
                                                second_line_to_append.append(round(float(k), 4))
                                    txt_output.append(second_line_to_append)
                                    second_line_to_append = []
                                print('Analysis is complete. Break.')
                                break
                                
                        if len(binary_line) > 0:
                            if string_line[0] == 'timestep':
                                frame += 1
                                frame_data = []
                                cell_parameters = []
                                cell_matrix = []
                                frame_data_line = []
                                timestep_flag = True
        else:
            pass
        
def window_sizes(mol_file):
    return sub_main(MOL(mol_file).atom_list)