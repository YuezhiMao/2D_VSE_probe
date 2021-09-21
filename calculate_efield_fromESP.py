#! /usr/bin/env python

import os, glob, sys, re, math
import numpy as np
import xyzgeom
import pandas as pd
from optparse import OptionParser

def ParseInput(ArgsIn):
    UseMsg = '''
    Usage: python [script] [options] [target_dir] [xyz_dir]
    '''
    parser = OptionParser(usage=UseMsg)
    parser.add_option('--cutoff_values', dest='cutoff_values', action='callback', type='string', default=None, callback=string_sp_callback, help='specify multiple cutoff values (default: 7A)')
    parser.add_option('--probe_atoms', dest='probe_atoms', action='callback', type='string', default=None, callback=string_sp_callback, help='the indices of C, D, O atoms (default 20, 21, 22)')

    options, args = parser.parse_args(ArgsIn)
    if len(args) < 3:
        parser.print_help()
        sys.exit(0)
    if options.probe_atoms != None and len(options.probe_atoms) != 3:
        print "Has to specify the index of three probe atoms (C, D, O)"
        parser.print_help()
        sys.exit(0)

    return options, args

def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))

def parse_esp(esp_file, tmpfile):
   fr = open(tmpfile, 'r')
   raw_data = []
   count = 0
   for line in fr.readlines():
      l_sp = line.split()
      if l_sp[0] == 'TotESP':
         count += 1
      if len(l_sp) == 2 and re.search('(\d+)', l_sp[0]) != None:
         raw_data.append(float(l_sp[1]))
   fr.close()
   if count == 2:
      esp = np.array(raw_data[len(raw_data)/2 : ])
   else:
      esp = np.array(raw_data)
   np.savetxt(esp_file, esp, fmt="%.7f")

def generate_all_esp_files():
   outfile_list = glob.glob('*.out')
   for outfile in outfile_list:
      tmpfile = "tmp.txt"
      command = "sed -n " + "'/TotESP/,/Ele EField/p' " + outfile + " > " + tmpfile
      os.system(command)
      esp_file = outfile + ".esp"
      parse_esp(esp_file, tmpfile)
   os.system('rm ' + tmpfile)

def get_bond_length(xyzfile, atom_idx1, atom_idx2):
    atomlist, coords = xyzgeom.parse_xyz_file(xyzfile)
    bond_length = xyzgeom.compute_distance(coords, atom_idx1, atom_idx2)
    return bond_length

def get_bond_length_all_frames(xyz_path, atom_idx_list):
    bond_length = {}
    xyzfile_list = glob.glob(xyz_path+'/'+'r*.xyz')
    for xyzfile in xyzfile_list:
        frame = int(re.search("frame_([^_]+).xyz", xyzfile).group(1))
        bond_length[frame] = {}
        bond_length[frame]["CD"] = get_bond_length(xyzfile, atom_idx_list[0], atom_idx_list[1]) * 1.88973  #convert to au
        bond_length[frame]["CO"] = get_bond_length(xyzfile, atom_idx_list[0], atom_idx_list[2]) * 1.88973
    return bond_length

def get_esp_single_cutoff(data_esp, target_dir, cutoff, atom_idx_list):
    esp_file_list = glob.glob(target_dir+"/r{:.2f}_frame*.esp".format(cutoff))
    data_esp[cutoff] = {}
    for esp_file in esp_file_list:
        data = np.loadtxt(esp_file)
        frame = int(re.search("frame_([^_]+)_", esp_file).group(1))
        data_esp[cutoff][frame] = {}
        for idx in atom_idx_list:
            data_esp[cutoff][frame][idx] = data[idx-1] #unlike e-field, esp should just be a single number
    return data_esp

def calculate_efield_from_esp(data_esp, bond_length, target_dir, atom_idx_list):
   au_to_MVcm = 5.142E+3
   idx_C, idx_D, idx_O = atom_idx_list[0], atom_idx_list[1], atom_idx_list[2]
   data_efield_stat = {}
   for cutoff in data_esp:
      efield_data_file = target_dir + '/' + "efield_per_frame_"+str(cutoff)+"A.csv"
      fw = open(efield_data_file, 'w')
      fw.write("frame,E_CO,E_CD\n")
      data_efield_stat[cutoff] = {}
      sum_E_CO = sum2_E_CO = 0.0
      sum_E_CD = sum2_E_CD = 0.0
      for frame in sorted(data_esp[cutoff]):
         r_CO = bond_length[frame]["CO"]
         r_CD = bond_length[frame]["CD"]
         V_C = data_esp[cutoff][frame][idx_C]
         V_D = data_esp[cutoff][frame][idx_D]
         V_O = data_esp[cutoff][frame][idx_O] 
         E_CO = (V_C - V_O) / r_CO * au_to_MVcm
         E_CD = (V_C - V_D) / r_CD * au_to_MVcm
         fw.write("%d,%.3f,%.3f\n" %(frame, E_CO, E_CD))
         sum_E_CO += E_CO
         sum2_E_CO += E_CO ** 2
         sum_E_CD += E_CD
         sum2_E_CD += E_CD ** 2
      n_frames = len(data_esp[cutoff])
      print "Cutoff = %.1f A, Number of frames: %d" %(cutoff, n_frames)
      fw.close()

      data_efield_stat[cutoff]["avg_E_CO"] = sum_E_CO 
      data_efield_stat[cutoff]["avg_E_CO"] = sum_E_CO / float(n_frames)
      data_efield_stat[cutoff]["std_E_CO"] = math.sqrt(sum2_E_CO / float(n_frames) - data_efield_stat[cutoff]["avg_E_CO"]**2)
      data_efield_stat[cutoff]["se_E_CO"] = data_efield_stat[cutoff]["std_E_CO"] / math.sqrt(float(n_frames))
      data_efield_stat[cutoff]["avg_E_CD"] = sum_E_CD / float(n_frames)
      data_efield_stat[cutoff]["std_E_CD"] = math.sqrt(sum2_E_CD / float(n_frames) - data_efield_stat[cutoff]["avg_E_CD"]**2)
      data_efield_stat[cutoff]["se_E_CD"] = data_efield_stat[cutoff]["std_E_CD"] / math.sqrt(float(n_frames))

   return data_efield_stat


options, args = ParseInput(sys.argv) 
curdir = os.getcwd()
print "Generating all ESP files from Q-Chem output..."
target_dir = args[1]
os.chdir(target_dir)
generate_all_esp_files()
os.chdir(curdir)

#probe atoms
if options.probe_atoms != None:
   atom_idx_list = []
   for iatm in options.probe_atoms:
      atom_idx_list.append(int(iatm))
else:
   atom_idx_list = [20, 21, 22]

print "Calculating the C-O and C-D directions in each snapshot"
xyz_path = args[2]
if xyz_path[:-1] == '/':
    xyz_path = xyz_path[:-1]
bond_length = get_bond_length_all_frames(xyz_path, atom_idx_list)

print "Getting the ESP values on relevant atoms"
data_esp = {}
cutoff_values = []
if options.cutoff_values == None:
   cutoff_values = [7]
else:
   [cutoff_values.append(int(val)) for val in options.cutoff_values]

for cutoff in cutoff_values:
   get_esp_single_cutoff(data_esp, target_dir, cutoff, atom_idx_list)

print "Calculating E-field along bond directions"
data_efield_stat = calculate_efield_from_esp(data_esp, bond_length, target_dir, atom_idx_list)

#print E-field on chemical bonds
df_efield = pd.DataFrame.from_dict(data_efield_stat, orient="index")
df_efield.index.name = 'cutoff'
for cutoff in sorted(cutoff_values):
   if cutoff == 0:
      continue
   print "E-field along C=O bond (cutoff = %dA):" %cutoff
   print "Mean: %.3f MV/cm, Std Err: %.3f MV/cm" %(df_efield["avg_E_CO"][cutoff], df_efield["se_E_CO"][cutoff])

   print "E-field along C-D bond (cutoff = %dA):" %cutoff
   print "Mean: %.3f MV/cm, Std Err: %.3f MV/cm" %(df_efield["avg_E_CD"][cutoff], df_efield["se_E_CD"][cutoff])
df_efield = df_efield.reindex(sorted(df_efield.columns), axis=1)
df_efield.to_csv(target_dir + '/efield_fromESP_bond.csv')
