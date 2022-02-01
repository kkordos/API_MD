# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 13:24:45 2021

@author: Kostis
"""

#import mupif as mp
import os
import numpy as np

def write_xscfile(cell_vectors, filename):
  
  celline = cell_vectors.reshape(9).astype(str)
  
  name = os.path.splitext(os.path.abspath(filename))[0]
  with open(str(name) + '.xsc', 'w' ) as f:
    f.write('# NAMD extended system configuration output file\n')
    f.write('#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w\n')
    f.write("0 {} 0 0 0 0 0 0 0 0 0".format(' '.join(celline) ))
  
  return None

def write_pdbfile(boxlines, filename):
  
  with open(filename+'.pdb', 'w') as f:
    f.write(''.join(boxlines))
  
  return None

def write_pdbline(ia, atname, im, atpos):
  
  # this is the standard pdb format for generating the host input file
  # https://cupnet.net/pdb-format/
  atx, aty, atz = atpos[0], atpos[1], atpos[2]
  
  line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format('ATOM', ia+1, atname, ' ', 'MOL', ' ', im+1, ' ', atx,aty,atz, 1.00, 0.00, ' ', ' ')
  
  return line

#cell_vectors = np.array([ [100.,0.,0.], [0.,100.,0.], [0.,0.,100.]])

### Prototype API implementation ### (Part of a class with MuPIF datatype specifiactions)

self.grainState = mp.heavydata.HeavyDataHandle(id=mp.dataid.MiscID.ID_GrainState) # we need a specification of which host grain
 
# mode probably must be changed here 
# The host grain object must be loaded here
hostgrain = self.grainState.getData(mode='create',schemaName='grain',schemasJson=mp.heavydata.sampleSchemas_json)
 
# The cell dimensions and (identity) are needed here
grainCell = hostgrain.getTopology().getCell() 
grainName = hostgrain.getIdentity().getMaterial()
 
# Input generation for replace.py
pdblines = []
for im,m in enumerate(hostgrain.getMolecules()):
  for ia,a in enumerate(hostgrain.getAtoms()):
     
     atype  = a.getTopology().getType()
     atname = a.getTopology().getName()
     atpos  = a.getTopology().getPosition()
     
     pdbline = write_pdbline(ia, atname, im, atpos) # generates one single line in pdbfile
     pdblines.append(pdbline+'\n')

# This function generates a pdb file needed for replace.py    
write_pdbfile(pdblines, grainName)

# This function generates an xsc file needed for replace.py 
write_xscfile(grainCell, '{}'.format(grainName))  

#- Mol2 file ?

parmchk2 < .mol2 > .frcmod

1mol: in.tleap ( .frcmod, .mol2, leaprc.gaff )
1molbulk: in.tleap ( 1.frcmod, 1.mol2, bulk.pdb, leaprc.gaff )  
2molbulk: in.tleap ( 1.frcmod, 2.frcmod, 1.mol2, 2.mol2, mix.pdb, leaprc.gaff )

tleap -f in.tleap > .prmtop 

namd < namd.in

























# These are completely arbitrary. Just relying on Model1 and Model2 of Example11 
# to give you an idea of the sequence of services to create the necessary inputs for models 
  
# Here the dopant molecule must be loaded from database 
self.moleculeState = mp.heavydata.HeavyDataHandle(id=mp.dataid.MiscID.ID_MoleculeState) # we need a specification of which dopant Molecule
dopantmol = self.moleculeState.getData(mode='create',schemaName='molecule',schemasJson=mp.heavydata.sampleSchemas_json)

dopantmolName = dopantmol.getIdentity().getChemicalName()

pdblines = [] 
for ia,a in enumerate(dopantmol.getAtoms()):
  
  datype  = a.getTopology().getType()
  datname = a.getTopology().getName()
  datpos  = a.getTopology().getPosition()
  
  # this is the standard pdb format for generating the dopant input file
  pdbline = write_pdbline(ia, atname, 1, atpos)
  pdblines.append(pdbline+'\n')
  
write_pdbfile(pdblines, dopantmolName)

# *********************** up to this point code creates 2 pdb files and an xsc file *****************

#################### this must be part of the workflow ###################

while True:
   
  os.system('python3 replace.py > xx') #  <---- Specify the output of the program here
  os.system('python3 validate.py > yy') # <---- Specify the output of the program here
  
  # Assuming that validate.py returns True/False
  with open('yy', 'r') as valout:
    valout = valout.readline()
  
  if valout[0] == True:
    break
##########################################################################

# ************************** SET NEW GRAIN HERE ***************************

  
  
  
  
  