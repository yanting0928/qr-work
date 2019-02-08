from __future__ import division
from mmtbx import utils
from mmtbx import model
from mmtbx import map_tools
import iotbx.pdb
from iotbx import reflection_file_utils
from cStringIO import StringIO
import mmtbx
import os
import sys
from scitbx.math import superpose
from libtbx.utils import null_out
from iotbx.file_reader import any_file
from cctbx import miller
from scitbx.array_family import flex

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def get_fmodel(crystal_symmetry, reflection_files, xray_structure):
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
    xray_structure = xray_structure)
  fmodel.update_all_scales()
  return fmodel

def get_map(fmodel, resolution_factor=1./5):
  f_map = fmodel.electron_density_map().map_coefficients(
    map_type = "mFo-DFc",
    isotropize   = True,
    fill_missing = False)
  fft_map = f_map.fft_map(resolution_factor=resolution_factor)
  return fft_map.real_map_unpadded()

def move_residue_atoms(map_data, residue, crystal_symmetry):
  for atom in residue.atoms():
    sites_cart = atom.xyz
    map_best = None
    xyz_best = []
    x_site = -0.1 
    while x_site < 0.1:
      x_site += 0.01
      y_site = -0.1
      while y_site < 0.1:
        y_site += 0.01
        z_site = -0.1
        while z_site < 0.1:
          z_site += 0.01
          site_cart = [
            sites_cart[0]+x_site,
            sites_cart[1]+y_site,
            sites_cart[2]+z_site] 
          site_frac = crystal_symmetry.unit_cell().fractionalize(site_cart)
          map_value = map_data.tricubic_interpolation(site_frac)
          if (map_value > map_best):
            map_best = map_value           
            xyz_best = site_cart
    print xyz_best               
    atom.set_xyz(xyz_best)

def run(pdb_file_name, data_file_name):
  pdb_code = os.path.basename(pdb_file_name)[:4]
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  crystal_symmetry = pdb_inp.crystal_symmetry()
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry = crystal_symmetry)
  fmodel_ini = get_fmodel(
    crystal_symmetry = crystal_symmetry,
    reflection_files = data_file_name,
    xray_structure   = xray_structure)
  print  "Initial r_work=%6.4f r_free=%6.4f" % (fmodel_ini.r_work(),
    fmodel_ini.r_free())
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for conformer  in residue_group.conformers():
          for residue in conformer.residues():
            sel_int = residue.atoms().extract_i_seq()
            n_atoms = xray_structure.scatterers().size()
            sel_bool = flex.bool(n_atoms, sel_int)
            xrs_sel = xray_structure.select(~sel_bool)
            fmodel = fmodel_ini.deep_copy()
            fmodel.update_xray_structure(xray_structure = xrs_sel, 
              update_f_calc=True)
            print "r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free())
            map_data = get_map(fmodel = fmodel)
            move_residue_atoms(
              map_data         = map_data, 
              residue          = residue,
              crystal_symmetry = crystal_symmetry)
  pdb_hierarchy.write_pdb_file(file_name="%s_update.pdb"%pdb_code,
                             crystal_symmetry = crystal_symmetry)
                   
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
#    run(pdb_file_name="data_perturbed.pdb", data_file_name="data.mtz")
    run(pdb_file_name="p1-1akg_shift.pdb", data_file_name="p1-1akg.mtz")
