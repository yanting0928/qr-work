import os
import iotbx.pdb

pdb_inp = iotbx.pdb.input(file_name = "1bdw.pdb")
ph = pdb_inp.construct_hierarchy()
for rg in ph.residue_groups():
  #n_confs = len(rg.conformers())
  #if (n_confs > 1):
  ###  Is  contain conformers  ###
  if rg.have_conformers():
    conformer_container=[]
    for atom_group in rg.atom_groups():
      if (atom_group.occupancy() < 1):
        conformer_container.append(atom_group)
        print atom_group.occupancy() 
        print atom_group.resname
        print atom_group.confid()
