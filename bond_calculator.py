#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
import numpy as np
import math

# Load the PDB file
file_name = "nvt"
file = file_name + ".gro"
u = mda.Universe(file)

# Setting the parameter
atom_number = len(u.atoms)
frame_number = len(u.trajectory)
dt = u.trajectory.dt
output = []
CO2_bond_default = 1.16 #Angstrom
CO2_angle_default = 180 #Degree
H2O_bond_default = 0.9572 #Angstrom
H2O_angle_default = 104.52 #Degree

# F4 calculation between two oxygen atoms of the water molecule which in the cutoff range(0.35nm)
def angle( a1, a2, a3, bond_default, angle_default, box):        
    # Get the vector and angle
    a2 -= a1
    a3 -= a1
    a1 = [0, 0, 0]
    vector1 = a1 - a2
    vector2 = a1 - a3
    bond1 = np.linalg.norm(vector1)
    bond2 = np.linalg.norm(vector2)
    if bond1 > (10*bond_default):
        for i in range(0, 3):
            if abs(a2[i]-a1[i]) > (box[i]/2): 
                if a2[i]-a1[i] > 0: a2[i] -= box[i]
                else: a2[i] += box[i]
                vector1 = a1 - a2
                bond1 = np.linalg.norm(vector1)
    if bond2 > (10*bond_default):
        for i in range(0, 3):
            if abs(a3[i]-a1[i]) > (box[i]/2): 
                if a3[i]-a1[i] > 0: a3[i] -= box[i]
                else: a3[i] += box[i]
                vector2 = a1 - a3
                bond2 = np.linalg.norm(vector2)
    cos_theta = np.dot(vector1, vector2) / (bond1*bond2)
    
    bond1 = round(bond1, 4)
    bond2 = round(bond2, 4)
    bond1_error = round(((bond1-bond_default)/bond_default),4)
    bond2_error = round(((bond2-bond_default)/bond_default),4) 
    cos_theta = round(cos_theta, 4)
    theta = math.acos(cos_theta)*180.0/math.pi
    theta = round(theta, 4)
    theta_error = round(((theta-angle_default)/angle_default),4)

    return (bond1, bond2, theta, bond1_error, bond2_error, theta_error)


# F4 calculation of one frame
def one_frame_calcu( f, output_file):
    print("the calculating frame is",f,"\n")
    # Access information for each frame
    coordinates = u.trajectory[f].positions
    atom_names = u.atoms.names
    residue_names = u.atoms.resnames
    box = [u.trajectory[f].dimensions[0], u.trajectory[f].dimensions[1], u.trajectory[f].dimensions[2]]

   
    # Initialize the parameters
    CO2_cluster = []
    H2O_cluster = []
    
    # Record the molecular index 
    for i in range(0, atom_number):
        # Select the carbon atoms of the CO2 molecule
        if residue_names[i] == 'CO2' and atom_names[i] == 'C':
            CO2_cluster.append([i,i+1,i+2])
    
    # Calculate CO2 angle
    for p in range(0, len(CO2_cluster)):
        i = CO2_cluster[p][0]
        j = CO2_cluster[p][1]
        k = CO2_cluster[p][2]
        (b1, b2, theta, b1er, b2er, ther) = angle( coordinates[i], coordinates[j], coordinates[k], CO2_bond_default, CO2_angle_default, box)
        output_file.append("For CO2 molecule(index = {:<5d}, {:<5d}, {:<5d}), bond1 = {:.5f}(Angstrom) and error = {:+.3%}, \
bond2 = {:.5f}(Angstrom) and error = {:+.3%}, angle = {:.5f}(degree) and error = {:+.3%}.".format(i+1,j+1,k+1,b1,b1er,b2,b2er,theta,ther))
    
    
    for i in range(0, atom_number):
        # Select the oxygen atoms of the H2O molecule
        if residue_names[i] == 'SOL' and atom_names[i] == 'OW':
            H2O_cluster.append([i,i+1,i+2])
    
    # Calculate H2O angle
    for p in range(0, len(H2O_cluster)):
        i = H2O_cluster[p][0]
        j = H2O_cluster[p][1]
        k = H2O_cluster[p][2]
        (b1, b2, theta, b1er, b2er, ther) = angle( coordinates[i], coordinates[j], coordinates[k], H2O_bond_default, H2O_angle_default, box)
        output_file.append("For H2O molecule(index = {:<5d}, {:<5d}, {:<5d}), bond1 = {:.5f}(Angstrom) and error = {:+.3%}, \
bond2 = {:.5f}(Angstrom) and error = {:+.3%}, angle = {:.5f}(degree) and error = {:+.3%}.".format(i+1,j+1,k+1,b1,b1er,b2,b2er,theta,ther))
    
        


    
# Generate the parameters to run the multi-process
for f in range(0,len(u.trajectory)):
    output.append("Frame = {}, CO2: bond_default = 1.16(Angstrom), angle_default = 180(Degree), and H2O: bond_default = 0.9572(Angstrom), angle_default = 104.52(Degree)".format(f))
    one_frame_calcu( f, output)

    
# Output the F4-time array in the file
np.savetxt(file_name + "_angle.txt", output, fmt='%s')
