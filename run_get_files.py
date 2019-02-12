import os, sys
import iotbx.pdb
from iotbx import crystal_symmetry_from_any
from mmtbx.command_line import cif_as_mtz
from cctbx import crystal
from libtbx.utils import null_out
from scitbx.array_family import flex
from libtbx import easy_mp
import mmtbx.model
from libtbx import group_args
from cctbx import miller
import itertools
from mmtbx.utils import run_reduce_with_timeout
from libtbx.utils import null_out

log = sys.stdout
pdb_files = "/net/cci/pdb_mirror/pdb/"
hkl_files = "/net/cci/pdb_mirror/structure_factors/"

pdb_codes = [
"1AA5", "1AB1", "1ALZ", "1BYZ", "1CBN", "1D8G", "1DPL", "1DY5", "1EJG", "1EM0", "1ENN", "1ET1", "1ETL", 
"1ETM", "1ETN", "1FN8", "1FY4", "1FY5", "1G66", "1G6X", "1GCI", "1GDN", "1GWE", "1HHU", "1HHY", "1HJE", 
"1I0T", "1I1W", "1IUA", "1IX9", "1IXB", "1J8G", "1JXT", "1JXW", "1JXX", "1JXY", "1M24", "1M40", "1MC2", 
"1MUW", "1N55", "1N9B", "1NWZ", "1O56", "1OB6", "1OB7", "1OEW", "1P9G", "1PJX", "1PQ5", "1PQ7", "1R6J", 
"1SK5", "1SSX", "1UCS", "1US0", "1V6P", "1VYR", "1W0N", "1X6Z", "1X8P", "1X8Q", "1XVO", "1YK4", "1YWA", 
"2B97", "2BW4", "2DCG", "2DDX", "2DSX", "2F01", "2FMA", "2H5C", "2H5D", "2HS1", "2I16", "2I17", "2IDQ", 
"2IXT", "2IZQ", "2J8T", "2JFR", "2O7A", "2O9S", "2OL9", "2OV0", "2P74", "2PEV", "2PF8", "2PFH", "2PVE", 
"2PWA", "2PYA", "2QDV", "2QXW", "2VB1", "2VXN", "2WFI", "2WFJ", "2WUR", "2XU3", "2YGI", "2YGJ", "2YKZ", 
"2YL7", "3A2O", "3A38", "3A39", "3AL1", "3BCJ", "3D43", "3FIL", "3G21", "3GGK", "3IP0", "3K34", "3KLR", 
"3KS3", "3MFJ", "3MI4", "3NIR", "3NJW", "3O4P", "3P4J", "3PYP", "3Q8J", "3QPA", "3SBN", "3U7T", "3UI4", 
"3UI6", "3VOR", "3W5H", "3WDN", "3WVM", "3X2L", "3X2M", "3X32", "3X34", "3ZOJ", "3ZQV", "3ZR8", "3ZSJ", 
"3ZSK", "3ZTM", "4AYO", "4AYP", "4CE8", "4D5M", "4EA9", "4EIC", "4F1V", "4G13", "4GCA", "4HIF", "4HIG", 
"4HP2", "4HS1", "4I8G", "4I8H", "4I8J", "4I8K", "4I8L", "4IGS", "4K8M", "4LAU", "4LAZ", "4LB3", "4LB4", 
"4LBR", "4LBS", "4M7G", "4NAG", "4NPD", "4NSV", "4O6U", "4O8H", "4OCB", "4OY5", "4PSS", "4PSY", "4PTH", 
"4QXI", "4REK", "4RKV", "4TJZ", "4TKB", "4TKJ", "4TUT", "4U9H", "4UA6", "4UA7", "4UA9", "4UAA", "4UYR", 
"4WEE", "4X1A", "4XXG", "4Y9V", "4Y9W", "4ZM7", "5AKR", "5AL6", "5AVD", "5AVH", "5B28", "5D8V", "5D9E", 
"5DJ7", "5EMB", "5GJI", "5GV7", "5GV8", "5HB7", "5HBS", "5I5B", "5I6A", "5JQF", "5JZQ", "5KWM", "5L87", 
"5LP9", "5MAS", "5MN1", "5MNG", "5MNK", "5MNN", "5MRX", "5NFM", "5NQI", "5NW3", "5O99", "5OME", "5OQZ", 
"5TDA", "5VLE", "5WQQ", "5WQR", "5X9L", "5XUF", "5XVT", "5Y2S", "5YCE", "5YOK", "6ANM", "6ANN", "6B00", 
"6DIY", "6EIO", "6ETK", "6ETL", "6EVH", "6FMC"
]

pdb_codes = [c.lower() for c in pdb_codes]

aa_codes = [
"resname ALA",
"resname ARG",
"resname ASN",
"resname ASP",
"resname CYS",
"resname GLN",
"resname GLU",
"resname GLY",
"resname HIS",
"resname ILE",
"resname LEU",
"resname LYS",
"resname MET",
"resname MSE",
"resname PHE",
"resname PRO",
"resname SER",
"resname THR",
"resname TRP",
"resname TYR",
"resname VAL"
]

def get_stat(model):
  h = model.get_hierarchy()
  asc = h.atom_selection_cache()
  ss = " or ".join(aa_codes)
  ss = "(%s) and not element H"%ss
  sel = asc.selection(ss)
  model = model.select(sel)
  model.get_restraints_manager()
  ms = model.geometry_statistics().result()
  return ms

def process_model(pdb):
  pdb_inp = iotbx.pdb.input(file_name=pdb)
  # EXPDTA
  expdata = pdb_inp.get_experiment_type()
  if(expdata != "X-RAY DIFFRACTION"): return None
  # SPLIT
  if("SPLIT" in pdb_inp.record_type_counts().keys()): return None
  # MODEL
  xrs = None
  model = None
  try:
    model = mmtbx.model.manager(model_input = pdb_inp, log=null_out())
    xrs = model.get_xray_structure()
  except: return None
  # IS PROTEIN
  h = pdb_inp.construct_hierarchy()
  asc = h.atom_selection_cache()
  sel = asc.selection("water or single_atom_residue")
  h = h.select(~sel)
  #
  asc = h.atom_selection_cache()
  ss = " or ".join(aa_codes)
  sel = asc.selection(ss)
  frac_aa = sel.count(True)*100./sel.size()
  if(frac_aa<80.): return
  #
  r = pdb_inp.get_r_rfree_sigma()
  r_work = r.r_work
  r_free = r.r_free
  if(r_work is None and len(r.r_works)>0):
    r_work = min(r.r_works)
  if(r_free is None and len(r.r_frees)>0):
    r_free = min(r.r_frees)
  #
  hd_selection = xrs.hd_selection()
  frac_h = hd_selection.count(True)*100./hd_selection.size()
  #
  return group_args(
    xrs = xrs, 
    r_work=r_work, 
    r_free=r_free, 
    d_min=r.resolution,
    frac_h=frac_h, 
    frac_aa=frac_aa, 
    model = model,
    #ms=get_stat(model), 
    sz = model.get_number_of_atoms())
      
def process_data(hkl, crystal_symmetry):
  mas = iotbx.cif.reader(file_path=hkl).build_miller_arrays()
  labs = []
  F = ["FOBS","F"]
  I = ["IOBS","I"]
  R = ["R-free-flags"]
  FI = F+I
  accepted = F + I + R
  found_F = 0
  found_I = 0
  found_R = 0
  for (data_name, mas_) in mas.iteritems():
    sub_labs = []
    for ma in mas_.values():
      l=cif_as_mtz.get_label(miller_array=ma, output_r_free_label="R-free-flags")
      if(l is None): continue
      if(not l in accepted): continue
      ma.set_info(info = miller.array_info(labels=[l]))
      if(l in F): found_F += 1
      if(l in I): found_I += 1
      if(l in R): found_R += 1
      sub_labs.append(l)
    labs.append(sub_labs)
  # Filter out flag-less loops (only in case there are loops with flags)
  if(found_R>0):
    labs_ = []
    for sub_labs in labs:
      if(not R[0] in sub_labs): continue
      labs_.append(sub_labs)
    labs = labs_[:]
  # Filter out ambiguous (more than 1 same-kind item in the loop present)
  labs_ = []
  for sub_labs in labs:
    if(sub_labs.count(F[0])>1): continue
    if(sub_labs.count(F[1])>1): continue
    if(sub_labs.count(I[0])>1): continue
    if(sub_labs.count(I[1])>1): continue
    if(sub_labs.count(R[0])>1): continue
    labs_.append(sub_labs)
  labs = labs_[:]
  # Flatten the list
  labs = list(itertools.chain(*labs))
  # Select highest resolution F or I array
  d_min = 9999
  data_array = None
  for (data_name, mas_) in mas.iteritems():
    for ma in mas_.values():
      if(not ma.info().label_string() in FI): continue
      if(not crystal_symmetry.is_similar_symmetry(ma.crystal_symmetry())):
        return None
      if(ma.crystal_symmetry().unit_cell() is None): return None
      if(ma.d_min() < d_min):
        d_min = ma.d_min()
        dinfo = ma.info()
        ma = ma.resolution_filter(d_min=0.4, d_max=100.) # safeguard
        data_array = ma.deep_copy()
  ##
  ot = data_array.deep_copy()
  data_array = data_array.as_non_anomalous_array()
  data_array = data_array.merge_equivalents().array().customized_copy(
    crystal_symmetry = data_array)
  data_array = data_array.map_to_asu()
  data_array = data_array.remove_systematic_absences()
  data_array = data_array.set_info(dinfo)
  data_array = data_array.set_observation_type(ot)
  ##
  if(not is_writable(data=data_array, flags=None)): return None
  #
  return data_array

def is_writable(data, flags):
  result = True
  try:
    mtz_dataset = data.as_mtz_dataset(
      column_root_label=data.info().labels[0])
    if(flags is not None):
      mtz_dataset.add_miller_array(miller_array=flags,
        column_root_label=flags.info().labels[0])
    mtz_object = mtz_dataset.mtz_object()
    result = True
  except Exception, e:
    result = False
  return result

def write_mtz_file(data, flags, pdb_code):
  try:
    mtz_dataset = data.as_mtz_dataset(
      column_root_label=data.info().labels[0])
    if(flags is not None):
      mtz_dataset.add_miller_array(miller_array=flags,
        column_root_label=flags.info().labels[0])
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = "%s.mtz"%pdb_code)
  except Exception, e:
    print "Cannot write mtz:", pdb_code, str(e)

def file_name_pairs_sorted():
  ifn_p = open("/".join([pdb_files,"INDEX"]),"r")
  ifn_r = open("/".join([hkl_files,"INDEX"]),"r")
  cntr_matches = 0
  cntr_all = 0
  pdb_file_names = flex.std_string()
  hkl_file_names = flex.std_string()
  codes          = flex.std_string()
  sizes     = flex.double()
  for lp in ifn_p.readlines():
    cntr_all += 1
    #print cntr_all
    lp = lp.strip()
    pdb_file_name = "".join([pdb_files,lp])
    assert os.path.isfile(pdb_file_name)
    pdb_code = lp[-11:-7]
    #print pdb_code
    #XXXXXXXXXXXXXXXXXCCCCCCCCCCED
    #if cntr_all == 100: break
    if not pdb_code in pdb_codes: continue
    print pdb_code
    #XXXXXXXXXXXXXXXXXCCCCCCCCCCED
    lr = lp.replace("pdb","r")
    lr = lr.replace(".ent.gz","sf.ent.gz")
    hkl_file_name = "".join([hkl_files,lr])
    if(os.path.isfile(hkl_file_name)):
      cntr_matches+=1
      pdb_file_names.append(pdb_file_name)
      hkl_file_names.append(hkl_file_name)
      codes.append(pdb_code)
      sizes.append(os.path.getsize(hkl_file_name))
  sel = flex.sort_permutation(sizes)
  pdb_file_names = pdb_file_names.select(sel)
  hkl_file_names = hkl_file_names.select(sel)
  codes     = codes.select(sel)
  sizes     = sizes.select(sel)
  #for a,b,c in zip(pdb_file_names, hkl_file_names, sizes):
  #  print a,b,c
  print >> log, "cntr_matches:", cntr_matches
  print >> log, "cntr_all    :", cntr_all
  return pdb_file_names, hkl_file_names, codes

def run_one(arg):
  pdb, hkl, code = arg
  #
  m = process_model(pdb=pdb)
  if(m is None): return
  sys.stdout.flush()
  #
  mas = process_data(hkl=hkl, crystal_symmetry=m.xrs.crystal_symmetry())
  if(mas is None): return
  #
  of = open("%s.pdb"%code,"w")
  of.write(m.model.model_as_pdb())
  of.close()
  #
  diff = abs(m.d_min-mas.d_min())
  f = ""
  if(diff>0.01): f="<<<<<<<<"
  fmt = "%s %3.2f %3.2f %5.1f %5.1f %d"
  print fmt%(code, m.d_min, mas.d_min(), m.frac_h, m.frac_aa,
    m.sz), mas.info().label_string(), f
  write_mtz_file(data=mas, flags=None, pdb_code=code)

def run(NPROC=1):
  pdbs, hkls, codes = file_name_pairs_sorted()
  args = []
  for a,b,c in zip(pdbs, hkls, codes):
    arg = [a,b,c]
    args.append(arg)
  #
  if(NPROC>1):
    stdout_and_results = easy_mp.pool_map(
      processes    = NPROC,
      fixed_func   = run_one,
      args         = args,
      func_wrapper = "buffer_stdout_stderr")
  else:
    for arg in args:
      run_one(arg)
  return True

if (__name__ == "__main__"):
  run()
