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
import mmtbx.hydrogens
from mmtbx.utils import run_reduce_with_timeout
from libtbx import group_args
from libtbx import adopt_init_args
from cctbx import xray
from phenix.refinement import weight_xray_chem
import scitbx

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
    map_data=map_data, unit_cell=crystal_symmetry.unit_cell(), amplitude=0.1,
    increment=0.001)
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
    target_name    = "ls_wunit_k1",
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xray_structure)
  fmodel.update_all_scales()
  return fmodel

def get_map(fmodel, step, map_type="mFo-DFc", scale=False):
  f_map = fmodel.electron_density_map().map_coefficients(
    map_type     = map_type,
    isotropize   = True,
    fill_missing = False)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = f_map.unit_cell(),
    space_group_info = f_map.crystal_symmetry().space_group_info(),
    step             = step,
    symmetry_flags   = maptbx.use_space_group_symmetry)
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = f_map)
  if(scale): fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def map_peak_coordinate(args):
  atom_i_seq,pdb_hierarchy,xray_structure,crystal_symmetry,fmodel_ini,step=args
  atom = pdb_hierarchy.atoms()[atom_i_seq]
  sel_int = flex.size_t([atom.i_seq])
  n_atoms = xray_structure.scatterers().size()
  sel_bool = flex.bool(n_atoms, sel_int)
  xrs_sel = xray_structure.select(~sel_bool)
  fmodel = fmodel_ini.deep_copy()
  fmodel.update_xray_structure(xray_structure = xrs_sel,
    update_f_calc=True)
  map_data = get_map(fmodel = fmodel, step = step)
  shift, xyz_best = move_residue_atoms(
    map_data         = map_data,
    atom             = atom,
    crystal_symmetry = crystal_symmetry)
  if 0 and (shift.split()[0].strip() in ["0.100","-0.100"] or
            shift.split()[1].strip() in ["0.100","-0.100"] or
            shift.split()[2].strip() in ["0.100","-0.100"]):
    print "r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()), shift
  return xyz_best

def show(model, fmodel, prefix):
  sel = model.selection(string="protein")
  model_ = model.select(sel)
  r = model_.geometry_statistics(use_hydrogens=False).result()
  f="%s r_work=%6.4f r_free=%6.4f bond: %7.4f angle: %7.3f"
  rw = fmodel.r_work()
  rf = fmodel.r_free()
  string=f%(prefix, rw, rf, r.bond.mean, r.angle.mean)
  return group_args(string = string, r_work = rw, r_free = rf,
    bonds = r.bond.mean, angles = r.angle.mean)

def get_model(pdb_file_name):
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  model = mmtbx.model.manager(
    model_input       = pdb_inp,
    build_grm         = False,
    stop_for_unknowns = False,
    log               = null_out())
  hd_sel = model.get_hd_selection()
  model = model.select(~hd_sel)
  rr = run_reduce_with_timeout(
    parameters  = " -quiet -",
    stdin_lines = model.model_as_pdb())
  pdb_inp = iotbx.pdb.input(source_info = None, lines = rr.stdout_lines)
  model = mmtbx.model.manager(
    model_input       = pdb_inp,
    build_grm         = True,
    stop_for_unknowns = False,
    log               = null_out())
  sel = model.selection(string = "element H and not protein")
  return model.select(~sel)

class lbfgs(object):
  def __init__(self, fmodel,
                     number_of_iterations = 100,
                     gradient_only = True,
                     stpmax = 0.25):
    adopt_init_args(self, locals())
    self.gradient_only = gradient_only
    self.xray_structure = self.fmodel.xray_structure
    self.correct_special_position_tolerance = 1.0
    self.x = flex.double(self.xray_structure.n_parameters(), 0)
    self._scatterers_start = self.xray_structure.scatterers()
    self.target_functor = self.fmodel.target_functor()
    lbfgs_core_params = scitbx.lbfgs.core_parameters(
      stpmin = 1.e-9,
      stpmax = stpmax)
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = number_of_iterations)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator          = self,
      gradient_only             = gradient_only,
      line_search               = True,
      core_params               = lbfgs_core_params,
      termination_params        = lbfgs_termination_params,
      use_fortran               = False,
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True))
    self.apply_shifts()

  def apply_shifts(self):
    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = self.xray_structure.unit_cell(),
      scatterers     = self._scatterers_start,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    self.xray_structure.replace_scatterers(scatterers = scatterers_shifted)

  def compute_functional_and_gradients(self):
    self.apply_shifts()
    self.fmodel.update_xray_structure(
      xray_structure = self.xray_structure,
      update_f_calc  = True)
    tg = self.target_functor(compute_gradients=True)
    self.f = tg.target_work() # not really used if gradient_only=True
    self.g = tg.gradients_wrt_atomic_parameters().packed()
    return self.f, self.g

def run(pdb_file_name, data_file_name, step, nproc, use_lbfgs=True,
        use_map_match=True):
  pdb_code = os.path.basename(pdb_file_name)[:4]
  model = get_model(pdb_file_name = pdb_file_name)
  model.idealize_h_riding()
  crystal_symmetry = model.crystal_symmetry()
  pdb_hierarchy = model.get_hierarchy()
  xray_structure = model.get_xray_structure()
  fmodel_ini = get_fmodel(
    crystal_symmetry = crystal_symmetry,
    reflection_files = data_file_name,
    xray_structure   = xray_structure)
  if(fmodel_ini.r_work()*100.>20):
    print "  skip: too large starting Rwork, inputs must be bad."
    return # SKIP
  unit_cell = crystal_symmetry.unit_cell()
  tfofc = get_map(
    fmodel=fmodel_ini, step=0.2, map_type="2mFo-DFc", scale=True)
  stat_start = show(model = model, fmodel = fmodel_ini, prefix = "Start:  ")
  if(use_map_match):
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
      fmodel_ini,step) for i in range(len(atom_seq))]
    results = pool_map(
      func      = map_peak_coordinate,
      iterable  = args,
      processes = nproc)
    site_carts = flex.vec3_double(pdb_hierarchy.atoms().extract_xyz())
    select = flex.size_t([i for i in atom_seq])
    site_cart_shifted = flex.vec3_double(results)
    site_cart_new = site_carts.set_selected(select,site_cart_shifted)
    model.set_sites_cart(sites_cart = site_cart_new)
    model.idealize_h_riding()
    fmodel_ini.update_xray_structure(
      xray_structure = model.get_xray_structure(),
      update_f_calc  = True,
      update_f_mask  = True)
    fmodel_ini.update_all_scales()
    stat_final = show(model = model, fmodel = fmodel_ini, prefix = "Final 1:")
  if(use_lbfgs):
    fmodel_ini.xray_structure.scatterers().flags_set_grad_site(
      iselection = xray_structure.all_selection().iselection())
    for cycle in xrange(10):
      if(cycle%2==0): fmodel_ini.set_target_name(target_name="ml")
      else:           fmodel_ini.set_target_name(target_name="ls")
      minimized = lbfgs(fmodel = fmodel_ini)
      model.set_sites_cart(sites_cart = minimized.xray_structure.sites_cart())
      model.idealize_h_riding()
    fmodel_ini.update_xray_structure(
      xray_structure = model.get_xray_structure(),
      update_f_calc  = True,
      update_f_mask  = True)
    fmodel_ini.update_all_scales()
  stat_final2 = show(model = model, fmodel = fmodel_ini, prefix = "Final 2:")
  # Write output model
  fo = fmodel_ini.f_obs()
  d_min = fo.d_min()
  cmpl  = fo.completeness()
  dtpr_h = fo.data().size()/(model.size()*3)
  dtpr   = fo.data().size()/(model.select(~model.get_hd_selection()).size()*3)
  of = open("%s_refined.pdb"%code, "w")
  print >> of, "REMARK %s"%stat_start.string
  print >> of, "REMARK %s"%stat_final.string
  print >> of, "REMARK %s"%stat_final2.string
  print >> of, "REMARK d_min = %6.4f completeness = %6.4f"%(d_min, cmpl)
  print >> of, "REMARK dtpr(H) = %3.1f dtpr(no H) = %3.1f"%(dtpr_h, dtpr)
  print >> of, model.model_as_pdb()
  #
  return group_args(
    d_min       = d_min,
    cmpl        = cmpl,
    stat_start  = stat_start,
    stat_final  = stat_final,
    stat_final2 = stat_final2)

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
        #if code != "4u9h": continue # For debugging
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
    #steps = [i/100 for i in range(4,10)]+[i/100 for i in range(10,31)]
    for pdb_file, mtz_file, code in zip(pdbs, mtzs, codes):
      print code, "-"*75
      r = run(pdb_file_name=pdb_file, data_file_name=mtz_file, step=0.125,
        nproc=50, use_lbfgs=True, use_map_match=True)
      if(r is None): continue
      print r.stat_start.string, "%6.4f %6.4f"%(r.d_min, r.cmpl)
      print r.stat_final.string
      print r.stat_final2.string
      sys.stdout.flush()
  else:
    run(pdb_file_name="3nir_refine_001.pdb",
        data_file_name="3nir_refine_data.mtz", nproc=90)
