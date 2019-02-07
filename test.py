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

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def get_fmodel(crystal_symmetry,reflection_files,xray_structure):
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
  fmodel.update_all_scales(remove_outliers=False)
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
  fft_map = f_map.fft_map(resolution_factor=1/5.)
  return fft_map.real_map_unpadded()

def run(pdb_file_name,data_file_name):
  pdb_code = os.path.basename(pdb_file_name)[3:7]
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  resolution = pdb_inp.resolution()
  d_min = 0.9
  target = [' CD ',' NE ',' CZ ',' NH1',' NH2']
  if (resolution is not None and resolution <= d_min):
    print resolution
    crystal_symmetry = pdb_inp.crystal_symmetry()
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    xray_structure = pdb_hierarchy.extract_xray_structure(crystal_symmetry = crystal_symmetry)
    fmodel_ini = get_fmodel(crystal_symmetry = crystal_symmetry,
                        reflection_files = data_file_name,
                        xray_structure   = xray_structure)
    print  "Initial r_work=%6.4f r_free=%6.4f" % (fmodel_ini.r_work(),
                                            fmodel_ini.r_free())
    prefix = 0
    hierarchy = iotbx.pdb.hierarchy.root()
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
#                pdb_hierarchy_new.write_pdb_file(file_name = "update_pdb_%d.pdb"%prefix)
#                pdb_hierarchy.write_pdb_file(file_name = "pdb_%d.pdb"%prefix)
                xray_structure_new = pdb_hierarchy_new.extract_xray_structure(crystal_symmetry = crystal_symmetry)
                fmodel = fmodel_ini.deep_copy()
                prefix += 1
                fmodel.update_xray_structure(xray_structure = xray_structure_new)
                fmodel.update_all_scales(remove_outliers=False)
                map_data = get_map(fmodel = fmodel)
                rg = iotbx.pdb.hierarchy.residue_group()
                ag = iotbx.pdb.hierarchy.atom_group()
                ag.resname = "ARG"
                for atom in residue.atoms():
                  if atom.name in target:
                    print atom.name
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
                          site_cart = [sites_cart[0]+x_site,
                                        sites_cart[1]+y_site,
                                        sites_cart[2]+z_site] 
                          site_frac = xray_structure.unit_cell().fractionalize(site_cart)
                          map_value = map_data.tricubic_interpolation(site_frac)
                          if (map_value > map_best):
                            atom.set_xyz(site_cart)
                            map_best = map_value           
                            xyz_best = site_cart
                    print sites_cart
                    print map_best
                    print xyz_best               
                    print atom.xyz
                    ag.append_atom(atom.detached_copy()) 
                rg.resseq=str(prefix)
                rg.append_atom_group(ag)
                model = iotbx.pdb.hierarchy.model()
                chain = iotbx.pdb.hierarchy.chain()
                chain.id = "A"
                model.append_chain(chain)
                hierarchy.append_model(model)
                chain.append_residue_group(rg)
    hierarchy.write_pdb_file(file_name="%s.pdb"%pdb_code)
                   
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
