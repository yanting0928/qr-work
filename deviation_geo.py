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
from libtbx.easy_mp import pool_map

#def ccp4_map(cg, file_name, map_data):
#  from iotbx import mrcfile
#  mrcfile.write_ccp4_map(
#    file_name=file_name,
#    unit_cell=cg.unit_cell(),
#    space_group=cg.space_group(),
#    #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
#    #gridding_last=n_real,  # This causes a bug (map gets shifted)
#    map_data=map_data,
#    labels=flex.std_string([""]))
#
#def exercise():
#  pdb_str = """
#CRYST1   13.123   11.874   14.671  93.00  95.00  98.00 P 1
#SCALE1      0.076202  0.010710  0.007382        0.00000
#SCALE2      0.000000  0.085045  0.005569        0.00000
#SCALE3      0.000000  0.000000  0.068569        0.00000
#ATOM      1  CA  GLY A   1       5.000   5.000   5.000  1.00  5.00           C
#ATOM      2  C   GLY A   1       5.428   6.025   6.024  1.00  5.00           C
#END
#  """
#  from cctbx import maptbx
#  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
#  h = pdb_inp.construct_hierarchy()
#  residue = h.only_residue()
#  cs = pdb_inp.crystal_symmetry()
#  xrs = h.extract_xray_structure(crystal_symmetry = cs)
#  fc = xrs.structure_factors(d_min=0.8).f_calc()
#  cg = fc.crystal_gridding(
#    d_min             = fc.d_min(),
#    symmetry_flags    = maptbx.use_space_group_symmetry,
#    resolution_factor = 1./8)
#  fft_map = miller.fft_map(
#    crystal_gridding     = cg,
#    fourier_coefficients = fc)
#  fft_map.apply_sigma_scaling()
#  map_data = fft_map.real_map_unpadded()
#  #
#  ccp4_map(cg=cg, file_name="map.map", map_data=map_data)
#  h.write_pdb_file(file_name="two_atoms.pdb", crystal_symmetry = cs)
#  #
##  pdb_str = """
##CRYST1   13.123   11.874   14.671  93.00  95.00  98.00 P 1
##SCALE1      0.076202  0.010710  0.007382        0.00000
##SCALE2      0.000000  0.085045  0.005569        0.00000
##SCALE3      0.000000  0.000000  0.068569        0.00000
##ATOM      1  CA  GLY A   1       5.023   5.020   5.050  1.00  5.00           C
##ATOM      2  C   GLY A   1       5.400   6.017   6.000  1.00  5.00           C
##END
##  """
##  from cctbx import maptbx
##  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
##  h = pdb_inp.construct_hierarchy()
##  residue = h.only_residue()
#  #
#  move_residue_atoms(map_data=map_data, residue=residue, crystal_symmetry=cs)
#  h.write_pdb_file(file_name="two_atoms_updated.pdb", crystal_symmetry = cs)

def move_residue_atoms(map_data, atom, crystal_symmetry):
  sites_cart = atom.xyz
  map_best = -9999
  xyz_best = []
  inc = 0.001
  x_site = -0.1 
  while x_site < 0.1:
    x_site += inc
    y_site = -0.1
    while y_site < 0.1:
      y_site += inc
      z_site = -0.1
      while z_site < 0.1:
        z_site += inc
        site_cart = [
          sites_cart[0]+x_site,
          sites_cart[1]+y_site,
          sites_cart[2]+z_site] 
        site_frac = crystal_symmetry.unit_cell().fractionalize(site_cart)
        map_value = map_data.tricubic_interpolation(site_frac)
        if (map_value > map_best):
          map_best = map_value           
          xyz_best = site_cart[:]
  diff = list(flex.vec3_double([sites_cart])-flex.vec3_double([xyz_best]))[0]
  atom.set_xyz(xyz_best)
  return " ".join(["%8.3f"%i for i in diff]),xyz_best
    
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

def get_map(fmodel, resolution_factor=1./10):
  f_map = fmodel.electron_density_map().map_coefficients(
    map_type = "mFo-DFc",
    isotropize   = True,
    fill_missing = False)
  fft_map = f_map.fft_map(resolution_factor=resolution_factor)
  return fft_map.real_map_unpadded()

def map_peak_coordinate(args):
  atom_i_seq,pdb_hierarchy,xray_structure,crystal_symmetry,fmodel_ini = args
  atom = pdb_hierarchy.atoms()[atom_i_seq]
  sel_int = flex.size_t([atom.i_seq])
  n_atoms = xray_structure.scatterers().size()
  sel_bool = flex.bool(n_atoms, sel_int)
  xrs_sel = xray_structure.select(~sel_bool)
  fmodel = fmodel_ini.deep_copy()
  fmodel.update_xray_structure(xray_structure = xrs_sel,
    update_f_calc=True)
  map_data = get_map(fmodel = fmodel)
  shift, xyz_best = move_residue_atoms(
    map_data         = map_data,
    atom             = atom,
    crystal_symmetry = crystal_symmetry)
  print "r_work=%6.4f r_free=%6.4f"%(
    fmodel.r_work(), fmodel.r_free()), shift
  return  xyz_best

def run(pdb_file_name, data_file_name, nproc=8):
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
  get_class = iotbx.pdb.common_residue_names_get_class
  atom_seq=[]
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for conformer  in residue_group.conformers():
          for residue in conformer.residues():
            if(get_class(name=residue.resname) == "common_water"): continue
            for atom in residue.atoms():
              if(atom.element.strip().upper()=="H"): continue
              atom_seq.append(atom.i_seq)
  args = [(atom_seq[i],pdb_hierarchy,xray_structure,crystal_symmetry,
    fmodel_ini)  for i in range(len(atom_seq))] 
  results = pool_map( 
    func      = map_peak_coordinate,
    iterable  = args,
    processes = nproc)
  site_carts = flex.vec3_double(pdb_hierarchy.atoms().extract_xyz())
  select = flex.size_t([i for i in atom_seq])
  site_cart_shifted = flex.vec3_double(results) 
  site_cart_new = site_carts.set_selected(select,site_cart_shifted)
  pdb_hierarchy.atoms().set_xyz(site_cart_new)
  pdb_hierarchy.write_pdb_file(file_name="%s_update.pdb"%pdb_code,
                             crystal_symmetry = crystal_symmetry)
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry = crystal_symmetry)
  fmodel = get_fmodel(
    crystal_symmetry = crystal_symmetry,
    reflection_files = data_file_name,
    xray_structure   = xray_structure)
  print  "Final r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(),
    fmodel.r_free())
                   
if __name__ == '__main__':
#  exercise()
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
    run(pdb_file_name="1akg.pdb", data_file_name="1akg.mtz")
