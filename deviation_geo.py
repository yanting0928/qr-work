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
from libtbx.utils import null_out
from cctbx import maptbx

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
  xyz_best = maptbx.fit_point_3d_grid_search(site_cart=sites_cart, 
    map_data=map_data, unit_cell=crystal_symmetry.unit_cell(), increment=0.001)
  diff = list(flex.vec3_double([sites_cart])-flex.vec3_double([xyz_best]))[0]
  return " ".join(["%8.3f"%i for i in diff]), xyz_best

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

def get_map(fmodel, resolution_factor=1./10, map_type="mFo-DFc", scale=False):
  f_map = fmodel.electron_density_map().map_coefficients(
    map_type     = map_type,
    isotropize   = True,
    fill_missing = False)
  fft_map = f_map.fft_map(resolution_factor=resolution_factor)
  if(scale): fft_map.apply_sigma_scaling()
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
  if 0: # Good for debugging!
    print "r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()), shift
  return xyz_best

def show(model, fmodel, prefix):
  sel = model.selection(string="protein")
  model_ = model.select(sel)
  r = model_.geometry_statistics(use_hydrogens=False).result()
  f="%s r_work=%6.4f r_free=%6.4f bond: %7.4f angle: %7.3f"
  print f%(prefix, fmodel.r_work(), fmodel.r_free(), r.bond.mean, r.angle.mean)
  sys.stdout.flush()

def run(pdb_file_name, data_file_name, nproc):
  pdb_code = os.path.basename(pdb_file_name)[:4]
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  model = mmtbx.model.manager(
    model_input       = pdb_inp,
    build_grm         = True,
    stop_for_unknowns = False,
    log               = null_out())
  crystal_symmetry = model.crystal_symmetry()
  pdb_hierarchy = model.get_hierarchy()
  xray_structure = model.get_xray_structure()
  fmodel_ini = get_fmodel(
    crystal_symmetry = crystal_symmetry,
    reflection_files = data_file_name,
    xray_structure   = xray_structure)
  cmpl = fmodel_ini.f_obs().completeness()
  #if(cmpl<0.9):
  #  print "Low completeness:", cmpl
  #  return # SKIP
  if(fmodel_ini.r_work()*100.>20):
    print "  skip: too large starting Rwork, inputs must be bad."
    return # SKIP
  unit_cell = crystal_symmetry.unit_cell()
  tfofc = get_map(
    fmodel=fmodel_ini, resolution_factor=1./4, map_type="2mFo-DFc", scale=True)
  show(model = model, fmodel = fmodel_ini, prefix = "Start:")
  get_class = iotbx.pdb.common_residue_names_get_class
  atom_seq=[]
  for model_ in pdb_hierarchy.models():
    for chain in model_.chains():
      for residue_group in chain.residue_groups():
        for conformer  in residue_group.conformers():
          for residue in conformer.residues():
            if(get_class(name=residue.resname) == "common_water"): continue
            for atom in residue.atoms():
              if(atom.element.strip().upper()=="H"): continue
              map_at_center = tfofc.tricubic_interpolation(
                unit_cell.fractionalize(atom.xyz))
              if(map_at_center<1.0): continue
              atom_seq.append(atom.i_seq)
  args = [(atom_seq[i], pdb_hierarchy, xray_structure, crystal_symmetry,
    fmodel_ini) for i in range(len(atom_seq))]
  results = pool_map(
    func      = map_peak_coordinate,
    iterable  = args,
    processes = nproc)
  site_carts = flex.vec3_double(pdb_hierarchy.atoms().extract_xyz())
  select = flex.size_t([i for i in atom_seq])
  site_cart_shifted = flex.vec3_double(results)
  site_cart_new = site_carts.set_selected(select,site_cart_shifted)
  pdb_hierarchy.atoms().set_xyz(site_cart_new)
  pdb_hierarchy.write_pdb_file(file_name="%s_updated.pdb"%pdb_code,
    crystal_symmetry = crystal_symmetry)
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry = crystal_symmetry)
  fmodel = get_fmodel(
    crystal_symmetry = crystal_symmetry,
    reflection_files = data_file_name,
    xray_structure   = xray_structure)
  model.set_sites_cart(sites_cart = xray_structure.sites_cart())
  show(model = model, fmodel = fmodel, prefix = "Final:")

if __name__ == '__main__':
#  exercise()
  if 1:
    path = "/net/anaconda/raid1/afonine/work/high_res_survey/high_res_pdb_mtz/"
    pdbs  = flex.std_string()
    mtzs  = flex.std_string()
    codes = flex.std_string()
    sizes = flex.double()
    for pdb_file in os.listdir(path):
      if(pdb_file.endswith(".pdb")):
        code = pdb_file[:-4]
        #
        #if code != "1etm": continue # For debugging
        #
        pdb_file = "%s%s.pdb"%(path, code)
        mtz_file = "%s%s.mtz"%(path, code)
        assert os.path.isfile(mtz_file)
        assert os.path.isfile(pdb_file)
        pdb_file_result = "%s_updated.pdb"%code
        if(os.path.isfile(pdb_file_result)):
          print "SKIP (already processed):", pdb_file_result
          continue # Skip already done case
        hierarchy = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy()
        size = hierarchy.atoms().size()
        if(hierarchy.models_size()>1): continue # Skip multi-model files
        #print pdb_file, mtz_file, size
        pdbs .append(pdb_file)
        mtzs .append(mtz_file)
        codes.append(code)
        sizes.append(size)
    sel = flex.sort_permutation(sizes) # Order by size, from smallest to largest
    pdbs  = pdbs .select(sel)
    mtzs  = mtzs .select(sel)
    codes = codes.select(sel)
    sizes = sizes.select(sel)
    for pdb_file, mtz_file, code in zip(pdbs, mtzs, codes):
      print code, "-"*75
      run(pdb_file_name=pdb_file, data_file_name=mtz_file, nproc=50)
  else:
    run(pdb_file_name="3nir_refine_001.pdb",
        data_file_name="3nir_refine_data.mtz", nproc=90)
