import iotbx.pdb
from mmtbx.monomer_library import server
from mmtbx.monomer_library import idealized_aa
from mmtbx.rotamer import rotamer_eval
from scitbx.array_family import flex

mon_lib_server = server.server()


def run(pdb_hierarchy):
  ideal_dict = idealized_aa.residue_dict()
  pdb_atoms = pdb_hierarchy.atoms()

  selection = flex.bool(pdb_atoms.size(), True)
  partial_sidechains = []
  for chain in pdb_hierarchy.only_model().chains():
    for residue_group in chain.residue_groups():
      if(residue_group.atom_groups_size() != 1):continue
      for residue in residue_group.atom_groups():
        i_seqs = residue.atoms().extract_i_seq()
        residue_sel = selection.select(i_seqs)
        if(not residue.resname.lower() in ideal_dict.keys()): continue
        missing_atoms = rotamer_eval.eval_residue_completeness(
          residue          = residue,
          mon_lib_srv      = mon_lib_server,
          ignore_hydrogens = True)
        if(len(missing_atoms) > 0):
          return True
  return False 

#n_changed = extend_sidechains.extend_protein_model(
#      pdb_hierarchy,
#      mon_lib_server,
#      add_hydrogens=False,
#)
#  pdb_hierarchy.write_pdb_file(file_name = "is_complete.pdb")
#print n_changed
if __name__ == '__main__':
  pdb_inp = iotbx.pdb.input(file_name = "1kyc.pdb")
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  print   run(pdb_hierarchy = pdb_hierarchy) 
