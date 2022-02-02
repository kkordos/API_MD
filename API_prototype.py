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

# read pdb ----------------------------------------------
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


# read mol2 files -------------------------------------------------
for im,m in enumerate(hostgrain.getMolecules()):
  for ia,a in enumerate(hostgrain.getAtoms()):
     
     atype  = a.getTopology().getType()
     atname = a.getTopology().getName()
     atpos  = a.getTopology().getPosition()
     
     pc= a.getProperties().getPartialCharge()
     
     bond     = a.getProperties().getbond()
     bondType = a.getProperties().getbondType()   
         
     ATOMSmol2line = write_mol2line(ia, atname, atpos, atype, im, pc ) 
     ATOMSmol2lines.append(mol2line+'\n')
     write_mol2file(ATOMSmol2lines, grainName)
     BONDSmol2line = write_mol2line(ia, bond, bondType )
     BONDSmol2lines.append(BONDSmol2line+'\n')



# run tleap ---------------------------------------     
    tleaplines.append('source leaprc.gaff')
 
  for im,m in enumerate (hostgrain.getMolecules()):
    frcmod = m.getProperties().getFrcmod()
    mol2   = m.getIdentity().getMol2()
    
    tleapline = write_tleapline(frcmod, mol2)
    tleaplines.append(tleapline+'\n')
    
    tleaplines.append('bulk = loadPdb mix.pdb')

    os.system('tleap -f in.tleap') # > .prmtop 


# run namd ---------------------------------------------
    TEMP=%temp
    NUMSTEPS=%num
    MARGIN=%margin 
    
    
    namdlines.append('extendedSystem XSC')
    namdlines.append('parmfile PRMTOP')
    namdlines.append('coordinates mix.pdb')

    os.system ("charmrun +pN namd2 in.namd > out.namd")



  
  
  
  