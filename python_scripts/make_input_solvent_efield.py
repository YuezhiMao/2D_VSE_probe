#! /usr/bin/env python

import os, glob, re, sys
import subprocess as sp
import numpy as np
from optparse import OptionParser
import qrems, qmol

def ParseInput(ArgsIn):
   UseMsg = '''
   make_input_sp_geom [options] [xyz_path]
   '''
   parser = OptionParser(usage=UseMsg)
   parser.add_option('-a','--all',dest='all',action='store_true',default=False,help='run all the xyz files under the xyz_path')
   parser.add_option('-k','--keyword',dest='keyword',action='store',type='string',default=None,help='create inputs for xyz files containing the keyword')
   parser.add_option('-t','--target',dest='target',action='callback',callback=string_sp_callback,type='string',default=None,help='create input for certain xyz files')
   parser.add_option('-b','--basis',dest='basis',action='callback',callback=string_sp_callback, type='string', default='aug-cc-pvtz',help='The target basis (default is aug-cc-pVTZ)')
   parser.add_option('-m','--method',dest='method',action='callback',type='string', default=None, callback=string_sp_callback,help='The methods (density functionals) to use')
   parser.add_option('-i','--input_path',dest='input_path',action='store',type='string', default='input/',help='The directory storing the generated inputs')
   parser.add_option('--almo',dest='almo',action='store_true',default=False,help='Use ALMO method to compute E-field')
   parser.add_option('--charge',dest='charge',action='store',type='int',default=0,help='total charge of the system')
   parser.add_option('--mult',dest='mult',action='store',type='int',default=1,help='total multiplicity of the system') 
   parser.add_option('--coarse',dest='coarse',action='store_true',default=False,help='use less tight integral thresh and less fine grid for SCF calculations')
   parser.add_option('--efield',dest='efield',action='store_true',default=False,help='calculate electric field')
   parser.add_option('--both_field',dest='both_field',action='store_true',default=False,help='calculate electric field in both jobs')
   parser.add_option('--solute',dest='solute',action='store',type='int',default=22,help='number of atoms on solute')

   options, args = parser.parse_args(ArgsIn)
   if len(args) < 2 and (options.all or options.keyword!=None):
      parser.print_help()
      sys.exit(0)
   if options.both_field:
      options.efield = True
   if not options.all and options.target==None and options.keyword==None:
      print("The target directory must be specified: one or some or all")
      parser.print_help()
      sys.exit(1)
   return options, args

def string_sp_callback(option, opt, value, parser):
   setattr(parser.values, option.dest, value.split(','))

def generate_FRGM(XYZ, natoms_frgm1, xyz_path):
   natoms_total = XYZ.NAtom
   natoms_frgm2 = natoms_total - natoms_frgm1
   frgm_file = xyz_path + XYZ.Name + '.frgm'
   fw = open(frgm_file, 'w')
   fw.write('0\n1\n0 0\n1 1\n')
   fw.write('%d %d' %(natoms_frgm1, natoms_frgm2)) 
   fw.close()
   FRGM = qmol.FRGM(frgm_file)
   return FRGM

def Append_Dev_Section(fw, options):
   fw.write('\n')
   fw.write('$development\n')
   if options.almo:
      fw.write('almo_efield  1\n')
   else:
      fw.write('embedding_early_stop  1\n')
   fw.write('$end\n\n')

def Append_Rem_Frgm(fw):
   fw.write('$rem_frgm\n')
   fw.write('SCF_CONVERGENCE  4\n')
   fw.write('$end\n\n')

def XYZ_to_Input(fw, parsed_XYZ, FRGM, myrems_1, myrems_2, options):
   #job 1
   qmol.WriteMolecule_Frgm(fw, parsed_XYZ, FRGM)
   qrems.AppendRem(fw, myrems_1)
   Append_Dev_Section(fw, options)
   if options.almo:
      Append_Rem_Frgm(fw)
   fw.write("@@@\n\n")
   #job 2
   qmol.WriteMolecule_supersub(fw, parsed_XYZ, FRGM, monomer=2, ghost=True)
   qrems.AppendRem(fw, myrems_2)
   if options.almo:
      Append_Dev_Section(fw, options)

options, args = ParseInput(sys.argv)
xyz_path =''
if len(args) > 1:
   xyz_path = args[1]
   if xyz_path[-1:] != '/':
      xyz_path += '/'

#determine xyz_file
curdir = os.getcwd()
xyzfile_list = []
if options.all:
   xyzfile_list = glob.glob(xyz_path+'*.xyz')
elif options.keyword:
   xyzfile_list = glob.glob(xyz_path+'*'+options.keyword+'*.xyz')
if options.target!=None:
   for xyz_file in options.target:
      if xyz_file not in xyzfile_list:
         xyzfile_list.append(xyz_file)

#check the input_path
input_path = options.input_path
if input_path[-1:] != '/':
   input_path += '/'
if not os.path.exists(input_path):
   sp.call(['mkdir', input_path])

rem_file = 'rem_stdscf'
myrems_1 = qrems.ParseRems(rem_file)
myrems_2 = qrems.ParseRems(rem_file)

for method in options.method:
   for basis in options.basis:
      qrems.set_rems_common(myrems_1, method, basis, options.coarse)
      qrems.set_rems_common(myrems_2, method, basis, options.coarse)
      #job 1 specific
      if options.almo:
         qrems.ModRem('FRGM_METHOD', 'STOLL', myrems_1)
         qrems.ModRem('SCFMI_MODE', '1', myrems_1)
         qrems.ModRem('SCF_GUESS', 'FRAGMO', myrems_1)
      else:
         qrems.ModRem('ENV_METHOD', method, myrems_1)
         qrems.ModRem('GEN_SCFMAN_EMBED', 'TRUE', myrems_1)
         qrems.ModRem('SPADE_PARTITION', 'TRUE', myrems_1)
         if options.both_field:
            qrems.ModRem('ESP_GRID', '0', myrems_1)
            qrems.ModRem('ESP_EFIELD', '1', myrems_1)
      #job 2 specific
      if options.almo:
         qrems.ModRem('SCF_GUESS', 'READ_DEN', myrems_2)
      else:
         qrems.ModRem('SCF_GUESS', 'READ', myrems_2)
      qrems.ModRem('GEN_SCFMAN', 'FALSE', myrems_2)
      qrems.ModRem('SKIP_SCFMAN', 'TRUE', myrems_2)
      if options.efield:
         qrems.ModRem('ESP_GRID', '0', myrems_2)
         qrems.ModRem('ESP_EFIELD', '1', myrems_2)
      for xyz_file in xyzfile_list:
         parsed_XYZ = qmol.XYZ(xyz_file)
         FRGM = generate_FRGM(parsed_XYZ, options.solute, xyz_path)
         if options.almo:
            inputfile = input_path+parsed_XYZ.Name+'_almo_'+method+'_'+qrems.basis_abbr(basis)+'.in' 
         else:
            inputfile = input_path+parsed_XYZ.Name+'_frozen_orb_'+method+'_'+qrems.basis_abbr(basis)+'.in' 
         fw = open(inputfile, 'w')
         XYZ_to_Input(fw, parsed_XYZ, FRGM, myrems_1, myrems_2, options)
         fw.close()
