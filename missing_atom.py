from mmtbx.building import extend_sidechains
import os
import iotbx.pdb
from mmtbx.monomer_library import server

mon_lib_server = server.server()

pdb_inp = iotbx.pdb.input(file_name = "1kyc.pdb")
pdb_hierarchy = pdb_inp.construct_hierarchy()

n_changed = extend_sidechains.extend_protein_model(
      pdb_hierarchy,
      mon_lib_server,
      add_hydrogens=False,
)
pdb_hierarchy.write_pdb_file(file_name = "is_complete.pdb")
print n_changed
