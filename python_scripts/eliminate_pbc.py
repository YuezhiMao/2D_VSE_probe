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
    parser.add_option('--ref_atom',dest='ref_atom',action='store',type='int',default=2,help='the reference atom position on solute (default: 2)')

    options, args = parser.parse_args(ArgsIn)
    if len(args) < 3:
        parser.print_help()
        sys.exit(0)

    if not options.all and options.target==None and options.keyword==None:
        print("The target directory must be specified: one or some or all")
        parser.print_help()
        sys.exit(0)

    return options, args


def get_frame_index(xyzfile):
    idx = re.search('frame_(\d+).xyz', xyzfile).group(1)
    #print "frame " + idx
    return int(idx)

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
        ref_atom = options.ref_atom
        apply_pbc_shift_atomref(xyzfile, ref_atom, box_size)    
    os.chdir(curdir)
