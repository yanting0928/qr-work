import os
import iotbx.pdb
from cctbx.array_family import flex

def have_conformers(ph):
  for model in ph.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        if residue_group.have_conformers():
          return True
  return False

if __name__ == "__main__":
  pdb_inp = iotbx.pdb.input(file_name = "1yjp.pdb")
  ph = pdb_inp.construct_hierarchy()
  if have_conformers(ph=ph):
    ph.remove_alt_confs(always_keep_one_conformer=True) 
    ph.write_pdb_file(file_name="tst.pdb")
  else:
    print "no"
