from __future__ import division
from mmtbx import model
import mmtbx
import iotbx.pdb
import os
import sys
from libtbx.utils import null_out
from scitbx.array_family import flex
from cctbx import geometry_restraints

def get_model(pdb_file_name,selection):
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  params =  mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.flip_symmetric_amino_acids = True
  model = mmtbx.model.manager(
    model_input               = pdb_inp,
    pdb_interpretation_params = params,
    stop_for_unknowns         = False,
    log                       = null_out())
  sel = model.selection(string = selection)
  return model.select(sel) 

def run(pdb_file_name):
  pdb_code = os.path.basename(pdb_file_name)[:4]
  resname = "ARG"
  atom_name = ["CZ", "NE", "NH1", "NH2", "CD"]
  selection = "resname " + resname + " and name " + " or name ".join(i for i in atom_name)
  print selection
  model = get_model(pdb_file_name = pdb_file_name, selection = selection)
  pdb_hierarchy = model.get_hierarchy()
  pdb_hierarchy.write_pdb_file(file_name = "test.pdb")
  grm = model.get_restraints_manager().geometry
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  atoms = pdb_hierarchy.atoms()
  for proxy in grm.angle_proxies:
    angle_proxy = proxy.i_seqs
    atom0 = atoms[angle_proxy[0]]
    atom1 = atoms[angle_proxy[1]]
    atom2 = atoms[angle_proxy[2]]
    angle = geometry_restraints.angle(
      sites_cart = sites_cart,
      proxy      = proxy)
    print "%s-%s-%s angle:"%(atom0.name,atom1.name,atom2.name),angle.angle_model
  for proxy in grm.dihedral_proxies:
    dihedral_proxy = proxy.i_seqs
    atom0 = atoms[dihedral_proxy[0]]
    atom1 = atoms[dihedral_proxy[1]]
    atom2 = atoms[dihedral_proxy[2]]
    atom3 = atoms[dihedral_proxy[3]]
    dihedral = geometry_restraints.dihedral(
      sites_cart = sites_cart,
      proxy = proxy)
    print "%s-%s-%s-%s dihedral:"%(atom0.name,atom1.name,atom2.name,atom3.name),dihedral.angle_model

if __name__ == '__main__':
#  exercise()
  if 0:
    path = "/net/anaconda/raid1/afonine/work/high_res_survey/high_res_pdb_mtz/"
    for pdb_file in os.listdir(path):
      if(pdb_file.endswith("_update.pdb")):
        code = pdb_file[:-4]
        #
        #if code != "1etm": continue # For debugging
        #
        pdb_file = "%s%s.pdb"%(path, code)
        assert os.path.isfile(pdb_file)
        pdb_file_result = "%s_analyze.pdb"%code
        if(os.path.isfile(pdb_file_result)):
          print "SKIP (already processed):", pdb_file_result
          continue # Skip already done case
        print code, "-"*75
        run(pdb_file_name=pdb_file)
  else:
    run(pdb_file_name="1l37.pdb")
