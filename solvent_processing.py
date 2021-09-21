#! /usr/bin/env python

import mdtraj as mdt
import numpy as np
import subprocess as sp
import os, sys
from optparse import OptionParser
import time

def ParseInput(ArgsIn):
    UseMsg = '''
    Usage: python [script] [options] [pdb_file] [target_dir] [dim_file]
    '''
    parser = OptionParser(usage=UseMsg)
    parser.add_option('--start_idx', dest='start_idx', action='store', type='int', default=0, help='index of the starting frame (default 0)')
    parser.add_option('--cutoff_values', dest='cutoff_values', action='callback', type='string', default=None, callback=string_sp_callback, help='specify multiple cutoff values (default: 0A and 7A)')
    parser.add_option('--stride', dest='stride', type='int', action='store', default=50, help='interval between two selected frames')
    parser.add_option('--water', dest='do_water', action='store_true', default=False, help='set true when solvent is water')
    parser.add_option('--solute', dest='solute', type='string', action='store', default='CCX', help='the solute name in pdb file (default: CCX)')
    parser.add_option('--lig_name', dest='lig_name', action='store', default='LIG', type='string', help='specify the name of ligand')
    parser.add_option('--xtc', dest='xtc', action='store_true', default=False, help='generate the xtc file for the strided trajectory')

    options, args = parser.parse_args(ArgsIn)
    if len(args) < 3:
        parser.print_help()
        sys.exit(0)
    return options, args

def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))

def compute_solv_solute_distances(traj_location, dims_location, stride = 50, do_water = False, solute_name = 'CCX', lig_name = 'LIG', generate_xtc = False):
    
    """
    This function takes a pdb trajectory file of a solute molecule
    in solvent and calculate the distance between each solute atom
    and the COM of each solvent molecule
    
    Input
    traj_location: file path for the input trajectory. Needs to be 
    a pdb file
    do_water: solvent water or not
    
    Output
    traj: directly from mdtraj load
    distances: regularized distances between each solute atom and solvent COMs
    [nframes, n_solute_atoms, n_solvents]

    Example:
    compute_solv_solute_distances("CCX_water_100_frames.pdb", do_water = True)
    
    """
    

    # load the trajectory
    print("Loading trajectory ...")
    start = time.time()
    traj    = mdt.load(traj_location, stride=stride)
    nframes = len(traj)
    print("Done loading %d frames in %.1f minutes" %(nframes, (time.time()-start)/60))
    #print(traj.top.atom(0).element.symbol)
    print(traj.xyz.shape)

    if generate_xtc:
        print ("Save strided trajectory to xtc file")
        xtc_filename = traj_location[:-3]+'xtc'
        traj.save_xtc(xtc_filename)    
    
    # unit cell dimensions
    if nframes > 1:
        cell_dims = np.loadtxt(dims_location)[::stride,:] * 0.1 # nanometers
    else:
        cell_dims = np.zeros((1,3))
        cell_dims[0, :] = np.loadtxt(dims_location) * 0.1 
    traj.unitcell_lengths = cell_dims #overwrite this member variable of traj
    assert(cell_dims.shape == (nframes, 3))
    
    # pick solvent and solute indices
    if do_water:
        solvent_indices = traj.top.select('water')
    else:
        solvent_indices = traj.top.select('resname ' + lig_name)
    #solute_indices  = traj.top.select('resname ' + solute_name)
    solute_indices  = traj.top.select('resid == 0')
    print ("number of solute atoms: %d" %len(solute_indices))

    # number of atoms each solvent molecule contains
    natoms_per_solvent = len(traj.top.select('resid == 1'))
    print("natoms per solvent: %d" %natoms_per_solvent)
    nsolvent_molecules = int(len(solvent_indices)/natoms_per_solvent)
    print("total number of solvent molecules: %d" %(len(solvent_indices)/natoms_per_solvent))  # number of solvent molecules
    
    # start by calculating the COM of solvents using the "angle" method. This fixes the PBC problem
    solvent_com = np.empty((nframes, nsolvent_molecules, 3), dtype=float)
    ai = solvent_indices[0]

    while ai < solvent_indices[-1]:

        # sum up the contribution from each atom
        sum_M  = 0.0
        sum_ME = 0.0
        sum_MC = 0.0
        
        for iatm in range(natoms_per_solvent):
            mass_i  = traj.top.atom(ai+iatm).element.mass
            theta_i = (traj.xyz[:, ai+iatm, :] / traj.unitcell_lengths) * 2 * np.pi # shape nframes * 3
            Ei = np.cos(theta_i)
            Ci = np.sin(theta_i)
            sum_M += mass_i
            sum_ME += Ei * mass_i
            sum_MC += Ci * mass_i
        
        mean_E = sum_ME / sum_M
        mean_C = sum_MC /sum_M
        
        mean_theta = np.arctan2(-1*mean_C, -1*mean_E) + np.pi
        
        com = traj.unitcell_lengths * mean_theta / (2 * np.pi)
        
        #index in the solvent_com array
        index_com = int((ai - solvent_indices[0]) / natoms_per_solvent)
        
        solvent_com[:, index_com, :] = com
        ai += natoms_per_solvent
    
    # now calculate the displacements
    displacements = np.empty((nframes, len(solute_indices), nsolvent_molecules, 3), dtype=float)
    
    for i in range(len(solute_indices)):
        solute_index = solute_indices[i]
        # calculate the raw difference
        solute_locations = traj.xyz[:, solute_index, :]
        displacements[:, i, :, :] = solvent_com - solute_locations[:, np.newaxis, :]
    
    # apply the minimum image convention
    displacements = displacements - (np.rint(displacements / traj.unitcell_lengths[:, np.newaxis, np.newaxis, :]) * 
                                    traj.unitcell_lengths[:, np.newaxis, np.newaxis, :])
        
    # now pick values below a certain threshold
    distances = np.sqrt(np.sum(np.square(displacements), axis=3))
    return traj, distances

def write_solv_solute_xyz(traj, distances, sphere_radius, target_dir, do_water = False, start_idx = 0, solute_name = 'CCX', lig_name = 'LIG'):
    """
    This function removes all solvent molecules that don't meet a 
    cutoff radius threshold, and write the processed solute-solvent coordinates
    to xyz files
    
    Input
    traj: previously loaded mdtraj data
    distances: calculated distances between each solute atom and solvent COMs
    sphere_radius: radius criterion for solvent molecules. Needs to be in 
    angstrom
    target_dir: the place to store xyz files
    do_water: solvent water or not
    
    Output
    xyz files containing each frame of the original trajectory, with excess water molecules
    removed

    """
    #This part is repeated work
    nframes = len(traj)
    # pick solvent and solute indices
    if do_water:
        solvent_indices = traj.top.select('water')
    else:
        solvent_indices = traj.top.select('resname ' + lig_name)
    #solute_indices  = traj.top.select('resname ' + solute_name)
    solute_indices  = traj.top.select('resid == 0')

    # number of atoms each solvent molecule contains
    natoms_per_solvent = len(traj.top.select('resid == 1'))

    indices_d = np.argwhere(distances < (.1 * sphere_radius))  # 0.1 to convert to nanometers
    
    # start writing the out files
    if not os.path.exists(target_dir):
        sp.call(['mkdir', target_dir])
    current_frame = 0
    current_index = 0
    
    while current_frame < nframes:
        
        # collect the solvent indices for this frame
        all_solvents = []
        
        while current_index < len(indices_d) and indices_d[current_index][0] == current_frame:
            all_solvents.append(indices_d[current_index][2])
            current_index += 1
            
        all_solvents  = list(set(all_solvents))
        total_atoms = len(all_solvents) * natoms_per_solvent + len(solute_indices)
        cell_dims   = traj.unitcell_lengths[current_frame]

        if current_frame == 0:
            print ("n solvent molecules: %d" %len(all_solvents))
            print ("natoms_total: %d" %total_atoms)
        
        # open an xyz file for this particular frame
        with open(target_dir+"/r{:.2f}_frame_{:04n}.xyz".format(sphere_radius, current_frame+start_idx), "w") as out_file:
            
            # opening lines for the xyz file
            out_file.write("{:n}".format(total_atoms))
            out_file.write("\n")
            
            out_file.write("cell dimensions: {:.4f} {:.4f} {:.4f}".format(cell_dims[0],
                                                                          cell_dims[1],
                                                                          cell_dims[2]))
            out_file.write("\n")
        
            # first write down the solute molecule:
            for atom in solute_indices:
                out_file.write("{:s}  {:.5f}  {:.5f}  {:.5f}".format(traj.top.atom(atom).element.symbol,
                                                                 10 * traj.xyz[current_frame, atom, 0],
                                                                 10 * traj.xyz[current_frame, atom, 1],
                                                                 10 * traj.xyz[current_frame, atom, 2]))
                out_file.write("\n")
            
            # then the solvent molecules
            for solvent in all_solvents:
                
                rel_index = solvent * natoms_per_solvent + solvent_indices[0]
                
                for iatm in range(natoms_per_solvent):
                    #if current_frame == 0:
                    #    print (traj.top.atom(rel_index+iatm).element.symbol)
                    out_file.write("{:s}  {:.5f}  {:.5f}  {:.5f}".format(traj.top.atom(rel_index+iatm).element.symbol,
                                                                     10 * traj.xyz[current_frame, rel_index+iatm, 0],
                                                                     10 * traj.xyz[current_frame, rel_index+iatm, 1],
                                                                     10 * traj.xyz[current_frame, rel_index+iatm, 2]))
                    out_file.write("\n")
            
        current_frame += 1

options, args = ParseInput(sys.argv)
pdb_file = args[1]
target_dir = args[2]
dims_location = args[3]

if target_dir[-1] == '/':
    target_dir = target_dir[:-1]
if not os.path.exists(target_dir):
    sp.call(['mkdir', target_dir])
if options.cutoff_values == None:
    cutoff_values = [0.0, 7.0]
else:
    cutoff_values = []
    for cutoff in options.cutoff_values:
        cutoff_values.append(float(cutoff))

traj, distances = compute_solv_solute_distances(pdb_file, dims_location, options.stride, options.do_water, options.solute, options.lig_name, options.xtc) 
n_frames = len(traj)
for cutoff in cutoff_values:
    xyz_dir = target_dir + '/' + "{:d}A_{:d}".format(int(cutoff), n_frames-1)
    write_solv_solute_xyz(traj, distances, cutoff, xyz_dir, options.do_water, options.start_idx, options.solute, options.lig_name)
