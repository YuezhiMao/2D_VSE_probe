#! /usr/bin/env python

import os, sys, re, glob
import numpy as np
import xyzgeom
import subprocess as sp
from optparse import OptionParser

def ParseInput(ArgsIn):
    UseMsg = '''
    Usage: python [script] [options] [xyz_dir] [box_size_file]
    '''
    parser = OptionParser(usage=UseMsg)
    parser.add_option('-a','--all',dest='all',action='store_true',default=False,help='processing all the directories under the xyz_path')
    parser.add_option('-k','--keyword',dest='keyword',action='store',type='string',default=None,help='processing directories containing a given keyword')
    parser.add_option('-t','--target',dest='target',action='store_true',default=False,help='processing a specific xyz directory')
    parser.add_option('--stride',dest='stride',action='store',type='int',default=50,help='stride in the box size file (default 50)')
    parser.add_option('--ref_atom',dest='ref_atom',action='store',type='int',default=2,help='the reference atom position on solute')
    parser.add_option('--com',dest='com',action='store_true',default=False,help='shift based on solute COM (calculated using the angle method)')
    parser.add_option('--natoms_solv',dest='natoms_solv',action='store',default=0,type='int',help='number of atoms per solvent')

    options, args = parser.parse_args(ArgsIn)
    if len(args) < 3:
        parser.print_help()
        sys.exit(0)

    if not options.all and options.target==None and options.keyword==None:
        print("The target directory must be specified: one or some or all")
        parser.print_help()
        sys.exit(0)

    if options.com and options.natoms_solv == 0:
        print("Specify number of atoms per solvent if shifting atoms based on COM")
        parser.print_help()
        sys.exit(0)

    return options, args

def get_mass(atom_symbol):
    if atom_symbol.upper() == 'C':
        return 12.00
    elif atom_symbol.upper() == 'H':
        return 1.00783
    elif atom_symbol.upper() == 'O':
        return 15.99491
    elif atom_symbol.upper() == 'N':
        return 14.00307
    else:
        print ("Unsupported atomic symbol")
        sys.exit(0)


def compute_com(AtomList, coords, natoms, box_size):
    #angle method for COM under PBC
    sum_M  = 0.0
    sum_ME = np.zeros(3)
    sum_MC = np.zeros(3)
    for i in range(natoms):
        M_i = get_mass(AtomList[i])
        theta_i = coords[i] / box_size * 2 * np.pi
        Ei = np.cos(theta_i)
        Ci = np.sin(theta_i)
        sum_M += M_i
        sum_ME += Ei * M_i
        sum_MC += Ci * M_i
    mean_E = sum_ME / sum_M
    mean_C = sum_MC / sum_M
    mean_theta = np.arctan2(-1*mean_C, -1*mean_E) + np.pi
    com = box_size * mean_theta / (2 * np.pi) 
    return com

def get_frame_index(xyzfile):
    idx = re.search('frame_(\d+).xyz', xyzfile).group(1)
    #print "frame " + idx
    return int(idx)

def apply_pbc_shift(xyzfile, natoms_solute, natoms_per_solvent, box_size):
    AtomList, coords = xyzgeom.parse_xyz_file(xyzfile)
    com = compute_com(AtomList, coords, natoms_solute, box_size)
    dist = coords - com
    natoms_total = len(AtomList)
    #process solute
    for d in range(0, 3):
        d_max = np.max(dist[0:natoms_solute, d]) 
        d_min = np.min(dist[0:natoms_solute, d])
        if d_max * d_min < 0 and (d_max - d_min) > box_size / 2.0:
            for j in range(natoms_solute):
                if abs(d_min) > abs(d_max):
                    if dist[j][d] < 0.0:
                        coords[j][d] += box_size
                else:
                    if dist[j][d] > 0.0:
                        coords[j][d] -= box_size
    #process solvent
    si = natoms_solute
    while si < natoms_total:
        for d in range(0, 3):
            d_max = np.max(dist[si:si+natoms_per_solvent, d])
            d_min = np.min(dist[si:si+natoms_per_solvent, d])
            if d_max * d_min < 0 and (d_max - d_min) > box_size / 2.0:
                for j in range(natoms_per_solvent):
                    if abs(d_min) > abs(d_max):
                        if dist[si+j][d] < 0.0:
                            coords[si+j][d] += box_size
                    else:
                        if dist[si+j][d] > 0.0:
                            coords[si+j][d] -= box_size
        si += natoms_per_solvent 

    #shift everything
    for i in range(natoms_total):
        for d in range(3):
            if ((coords[i][d] - com[d]) > 0.5 * box_size):
                coords[i][d] -= box_size
            elif ((coords[i][d] - com[d]) < -0.5 * box_size):
                coords[i][d] += box_size
    #write the shifted coordinates
    xyzgeom.write_xyz_file(xyzfile, AtomList, coords)

def apply_pbc_shift_atomref(xyzfile, ref_atom, box_size):
    AtomList, coords = xyzgeom.parse_xyz_file(xyzfile)
    natoms_total = len(AtomList)
    for i in range(natoms_total):
        for d in range(3):
            if ((coords[i][d] - coords[ref_atom-1][d]) > 0.5 * box_size):
                coords[i][d] -= box_size
            elif ((coords[i][d] - coords[ref_atom-1][d]) < -0.5 * box_size):
                coords[i][d] += box_size
    #write the shifted coordinates
    xyzgeom.write_xyz_file(xyzfile, AtomList, coords)

options, args = ParseInput(sys.argv)
xyzdir = args[1]
boxsize_file = args[2]
all_box_sizes = np.loadtxt(boxsize_file)
if all_box_sizes.ndim == 1:
   box_size_list = np.zeros((1,1))
   box_size_list[0] = all_box_sizes[0]
else:
   box_size_list = all_box_sizes[::options.stride, 0] #already in Angstrom

if xyzdir[-1] == "/":
    xyzdir = xyzdir[:-1]
xyzdir_list = []
if options.all:
    xyzdir_list = glob.glob(xyzdir+'/*/')
elif options.keyword:
    xyzdir_list = glob.glob(xyzdir+'/*'+options.keyword+'*/')
elif options.target:
    xyzdir_list.append(xyzdir)

print (xyzdir_list)
for target_dir in xyzdir_list:
    curdir = os.getcwd()
    os.chdir(target_dir)
    xyzfile_list = glob.glob('*.xyz')
    for xyzfile in sorted(xyzfile_list):
        box_size = box_size_list[get_frame_index(xyzfile)]
        if options.com:
            natoms_solute = 22
            apply_pbc_shift(xyzfile, natoms_solute, options.natoms_solv, box_size) 
        else:
            ref_atom = options.ref_atom
            apply_pbc_shift_atomref(xyzfile, ref_atom, box_size)    
    os.chdir(curdir)
