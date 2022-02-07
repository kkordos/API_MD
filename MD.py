# -*- coding: utf-8 -*-

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

def write_pdbline(indexatom, atomname, indexmolecule, atomcoords):
  
  # this is the standard pdb format for generating the host input file
  # https://cupnet.net/pdb-format/
  atx, aty, atz = atpos[0], atpos[1], atpos[2]
  
  line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format('ATOM', ia+1, atname, ' ', 'MOL', ' ', im+1, ' ', atx,aty,atz, 1.00, 0.00, ' ', ' ')
  
  return line

#cell_vectors = np.array([ [100.,0.,0.], [0.,100.,0.], [0.,0.,100.]])

def write_tleapfile(lines)
    with open ('in.tleap', 'w') as f:
    # leaprc.gaff is included in amber tools when it is installed
    f.write ('source leaprc.gaff')
    # tleap needs the BOX pdb file to run and extract the prmtop file with the parameters for MD 
    f.write ('bulk=loadPdb BOX.pdb')
    f.write(''.join(lines))
            
    return None

def write_tleapline(frcmod, mol2):
    line='loadamberparams', frcmod, '\n MOL=',mol2
    
    return line 
    

grainState = mp.heavydata.HeavyDataHandle(id=mp.dataid.MiscID.ID_GrainState) # we need a specification of which host grain
#hostgrain = grainState.getData(mode='create',schemaName='grain',schemasJson=mp.heavydata.sampleSchemas_json)

# The cell dimensions and (identity) are needed here
grainCell = hostgrain.getStructure().getCell() 
grainName = hostgrain.getIdentity().getMaterial()
 
# Input generation
pdblines = []

# create pdb from GRAIN ----------------------------------------------
for indexmolecule,molecule in enumerate(hostgrain.getMolecules()):
  for indexatom,atom in enumerate(hostgrain.getAtoms()):
     
     atomtype    = atom.getStructure().getType()
     atomname    = atom.getStructure().getName()
     atomcoords  = atom.getStructure().getPosition()
     
     pdbline = write_pdbline(indexatom, atomname, indexmolecule, atomcoords) # generates one single line in pdbfile
     pdblines.append(pdbline+'\n')

# This function generates a pdb file    
write_pdbfile(pdblines, grainName)

# This function generates an xsc file 
write_xscfile(grainCell, '{}'.format(grainName))  


# create mol2 from GRAIN ? -------------------------------------------------
for indexmolecule,molecule in enumerate(hostgrain.getMolecules()):
  for indexmolecule,atom in enumerate(hostgrain.getAtoms()):
     
     atomtype    = atom.getStructure().getType()
     atomname    = atom.getStructure().getName()
     atomcoords  = atom.getStructure().getPosition()
     
     partialcharge  = atom.getProperties().getPartialCharge()
     
     bond     = atom.getProperties().getbond()
     bondType = atom.getProperties().getbondType()   
         
     ATOMSmol2line = write_mol2line(indexatom, atomname, atomcoords, atomtype, indexmolecule, partialcharge ) 
     ATOMSmol2lines.append(mol2line+'\n')
     
     write_mol2file(ATOMSmol2lines, grainName)
     
     # we need to make mol2 bonds in a proper way
     BONDSmol2line = write_mol2line(indexbond, bond, bondType)
     BONDSmol2lines.append(BONDSmol2line+'\n')


# run tleap ---------------------------------------
# mol2 file from database ?

# take molecule species from database
    moleculeSpecies1frcmod = get_Molecule_frcmod_from_DB() 
    moleculeSpecies1mol2   = get_Molecule_mol2_from_DB() 
    moleculeSpecies2frcmod = get_Molecule_frcmod_from_DB() 
    moleculeSpecies2mol2   = get_Molecule_mol2_from_DB() 

# write the in.tleap file to run tleap as preprocessor to extract prmtop file
    tleapline = write_tleapline(moleculeSpecies1frcmod, moleculeSpecies1mol2)
    tleaplines.append(tleapline+'\n')
    
    tleapline = write_tleapline(moleculeSpecies2frcmod, moleculeSpecies2mol2)
    tleaplines.append(tleapline+'\n')
    
# run tleap
    os.system('tleap -f in.tleap') # > .prmtop 


# run namd ---------------------------------------------
# we need a service to get in.namd input file from database to run namd software for MD
  
# get namd input from database    
    namdIn  = get_namdInput_from_DB()   
    
# charmrun is used to run namd parallel at {N} cores 
    os.system("charmrun +p{} namd2 namdin > out.namd".format(NUMBER_OF_CPUS))
