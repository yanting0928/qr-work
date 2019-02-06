from __future__ import division
from mmtbx import utils
from mmtbx import model
import iotbx.pdb
from iotbx import reflection_file_utils
from cStringIO import StringIO
import mmtbx
import time
import os
import sys
import libtbx
from scitbx.math import superpose
from libtbx.utils import null_out
from iotbx.file_reader import any_file
from libtbx import group_args

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def get_fmodel(crystal_symmetry,reflection_files,model):
  hkl_file = any_file(file_name = reflection_files)
  rfs = reflection_file_server(
        crystal_symmetry = crystal_symmetry,
        reflection_files = [hkl_file.file_content])
  determine_data_and_flags_result = utils.determine_data_and_flags(
        reflection_file_server  = rfs,
        keep_going              = True,
        log                     = null_out())
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs,
      r_free_flags   = r_free_flags,
      xray_structure = model.get_xray_structure())
  fmodel.update_all_scales()
  return fmodel

def get_model(crystal_symmetry,pdb_hierarchy):
  model = mmtbx.model.manager(
    model_input               = None,
    pdb_hierarchy             = pdb_hierarchy,
    crystal_symmetry          = crystal_symmetry,
    log                       = null_out())
  return model

def move(ref,rg):
  reference_sites = ref.atoms().extract_xyz()
  moving = rg.atoms().extract_xyz()
  ls = superpose.least_squares_fit(
    reference_sites = reference_sites,
    other_sites     = moving)
  moved = ls.r.elems * moving + ls.t.elems
  rmsd = reference_sites.rms_difference(moved)
  print rmsd
  rg.atoms().set_xyz(moved) 
  
def model_select(pdb_hierarchy,select):
  asc = pdb_hierarchy.atom_selection_cache()
  sel = asc.selection(select)
  return pdb_hierarchy.select(sel)

def get_map(fmodel):
  f_map = fmodel.electron_density_map().map_coefficients(
            map_type = "mFo-DFc",
            isotropize   = True,
            fill_missing = False)
  fft_map = f_map.fft_map(resolution_factor=1/5)
  return fft_map.real_map_unpadded()

def run(pdb_file_name,data_file_name):
#  pdb_code = os.path.basename(pdb_file_name)[3:7]
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  resolution = pdb_inp.resolution()
  d_min = 0.9
  target = [' CD ',' NE ',' CZ ',' NH1',' NH2']
  if (resolution is not None and resolution <= d_min):
    print resolution
    crystal_symmetry = pdb_inp.crystal_symmetry()
    pdb_hierarchy = pdb_inp.construct_hierarchy()
#    model = get_model(crystal_symmetry = crystal_symmetry,
#                      pdb_hierarchy    = pdb_hierarchy)
#    fmodel = get_fmodel(crystal_symmetry = crystal_symmetry,
#                        reflection_files = data_file_name,
#                        model            = model)
#    print  "Initial r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(),
#                                            fmodel.r_free())
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          for conformer  in residue_group.conformers():
            for residue in conformer.residues():
              if (residue.resname == "ARG"):
                if not conformer.altloc: 
                  select = "not (chain " + chain.id +" and resseq" + residue.resseq + ")"
                else:
                  select = "not (chain " + chain.id +" and resseq " + residue.resseq + " and altid " + conformer.altloc + ")"
                print select
                pdb_hierarchy_new = model_select(pdb_hierarchy = pdb_hierarchy,
                                                        select = select) 
                model = get_model(crystal_symmetry = crystal_symmetry,
                                  pdb_hierarchy    = pdb_hierarchy_new)
                fmodel = get_fmodel(crystal_symmetry = crystal_symmetry,
                                    reflection_files = data_file_name,
                                    model            = model)
                print  "Initial r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(),
                                            fmodel.r_free())
                             

if __name__ == '__main__':
  if 0:
    path = "/home/pdb/pdb/"
    dpath = "/home/pdb/structure_factors"
    of = open("".join([path,"INDEX"]),"r")
    files = ["".join([path,f]).strip() for f in of.readlines()]
    of.close()
    of = open("/".join([dpath,"INDEX"]),"r")
    dfiles = [
      os.path.basename("".join([path,f]).strip())[1:5] for f in of.readlines()]
    of.close()
    for f in files:
      pdb_code = os.path.basename(f)[3:7]
      if (pdb_code in dfiles):
        d = "/".join([dpath,pdb_code[1:3],"r"+pdb_code+"sf.ent.gz"])
        run(pdb_file_name=f,data_file_name=d)
  else:
    run(pdb_file_name="1us0.pdb",data_file_name="1us0.mtz")
