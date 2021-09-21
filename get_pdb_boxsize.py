#! /usr/bin/env python

import os, sys, re, glob
import subprocess as sp
from optparse import OptionParser

if len(sys.argv) < 3:
   print ("Usage: python [script] [pdb_file] [dim_file]")
   sys.exit()

pdb_file = sys.argv[1]
dim_file = sys.argv[2]
command = 'grep \"CRYST1\" ' + pdb_file + ' | awk \'{print $2, $3, $4}\' > ' + dim_file
print (command) 
os.system(command)
