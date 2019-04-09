import os
from libtbx import easy_run
import iotbx.pdb
import mmtbx.model
from cctbx import uctbx
import libtbx
from mmtbx.monomer_library import server
import qrefine.clustering as clustering
from qrefine.utils import yoink_utils
from qrefine.plugin.yoink.pyoink import PYoink
from qrefine.fragment import fragments
from qrefine.charges import charges_class
import qrefine.completion as completion
from mmtbx.monomer_library import server
from mmtbx.monomer_library import idealized_aa
from mmtbx.rotamer import rotamer_eval
from scitbx.array_family import flex
from libtbx import group_args
from mmtbx import model
from libtbx.utils import null_out

qrefine = libtbx.env.find_in_repositories("qrefine")
mon_lib_server = server.server()

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

def keep_protein_only(pdb_hierarchy):
  asc = pdb_hierarchy.atom_selection_cache()
  ss = " or ".join(aa_codes)
  ss = "(%s) and not (element H or element D)"%ss
  pdb_hierarchy =  pdb_hierarchy.select(asc.selection(ss))
  if pdb_hierarchy.atoms_size() == 0: return None
  else                              : return pdb_hierarchy

def have_conformers(pdb_hierarchy):
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        if residue_group.have_conformers():
          return True
  return False

def box_pdb(pdb_hierarchy):
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = sites_cart,
    buffer_layer = 10)
  pdb_hierarchy.atoms().set_xyz(box.sites_cart)
  return group_args(
    pdb_hierarchy = pdb_hierarchy,
    cs            = box.crystal_symmetry())

def clusters(pdb_hierarchy):
  yoink_utils.write_yoink_infiles("cluster.xml",
                                  "qmmm.xml",
                                  pdb_hierarchy,
                                  os.path.join(qrefine,"plugin","yoink","dat"))
  pyoink=PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"),
                os.path.join(qrefine,"plugin","yoink","dat"),
                "cluster.xml")
  interaction_list,weight = pyoink.get_interactions_list()
  for e in interaction_list:
    e.sort()
  interaction_list.sort()
  cc=clustering.betweenness_centrality_clustering(interaction_list,maxnum_residues_in_cluster=3)
  print cc.get_clusters()
  return  cc.get_clusters()

def check_missing_atom(pdb_filename):
  pdb_inp = iotbx.pdb.input(file_name = pdb_filename)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
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

def add_hydrogens_using_ReadySet(pdb_hierarchy):
  from elbow.command_line.ready_set import run_though_all_the_options
  pdb_lines = pdb_hierarchy.as_pdb_string()
  rc = run_though_all_the_options(
    pdb_lines,
    [], # args
    hydrogens=True,
    reduce_args=['-BUILD', '-NUClear'],
    ligands=False,
    add_h_to_water=True,
    metals=False,
    output_file_name="reasyset.pdb",
    verbose=False,
    silent=True,
    )
  return rc['model_hierarchy']

def validate(pdb_str, threshold_bonds=0.02*4, threshold_angles=2.5*4):
  pdb_inp = iotbx.pdb.input(source_info = None, lines = pdb_str)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  model = mmtbx.model.manager(
    model_input               = pdb_inp,
    build_grm                 = True,
    stop_for_unknowns         = True,#False,
    pdb_interpretation_params = params,
    log                       = null_out())
  grm = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart()
  b_deltas = flex.abs(
    grm.get_all_bond_proxies()[0].deltas(sites_cart=sites_cart))
  b_outl = b_deltas.select(b_deltas>threshold_bonds)
  if(b_outl.size()>0): return None
  a_deltas = flex.abs(grm.get_all_angle_proxies().deltas(sites_cart=sites_cart))
  a_outl = a_deltas.select(a_deltas>threshold_angles)
  if(a_outl.size()>0): return None
  return pdb_str

def run(file_name,
        d_min=0.8,
        maxnum_residues_in_cluster=3,
        filter_non_protein_and_hd=True):
  prefix = os.path.basename(file_name)[3:7]
  print file_name, prefix
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  resolution = pdb_inp.resolution()
  data_type = pdb_inp.get_experiment_type()
  if(resolution is not None and resolution <= d_min):
    if(data_type in ["X-RAY DIFFRACTION","NEUTRON DIFFRACTION"]):
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      xray_structure = pdb_hierarchy.extract_xray_structure()
      xray_structure.convert_to_isotropic()
      pdb_hierarchy.adopt_xray_structure(xray_structure)
      if(filter_non_protein_and_hd):
        pdb_hierarchy = keep_protein_only(pdb_hierarchy = pdb_hierarchy)
      if(pdb_hierarchy is not None):
        if(have_conformers(pdb_hierarchy=pdb_hierarchy)):
          pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=True)
        pdb_hierarchy = add_hydrogens_using_ReadySet(
          pdb_hierarchy=pdb_hierarchy)
        box = box_pdb(pdb_hierarchy = pdb_hierarchy)
        fq = fragments(
          pdb_hierarchy              = box.pdb_hierarchy,
          maxnum_residues_in_cluster = maxnum_residues_in_cluster,
          debug                      = True,
          crystal_symmetry           = box.cs)
        print "clusters:", fq.clusters
        os.mkdir(prefix)
        cntr=0
        for fname in os.listdir(os.getcwd()):
          if fname.endswith("_cluster.pdb"):
            if check_missing_atom(pdb_filename = fname):
              os.remove(fname)
            else:
              cl_fn = "%s_%d.pdb"%(prefix, cntr)
              fo = open(cl_fn,"w")
              ph_i = iotbx.pdb.input(file_name=fname).construct_hierarchy()
              ph_i = completion.run(
                pdb_hierarchy    = ph_i,
                crystal_symmetry = box.cs,
                model_completion = False)
              pdb_str = ph_i.as_pdb_string(
                crystal_symmetry = box.cs, append_end=True)
              pdb_str = validate(pdb_str = pdb_str)
              if(pdb_str is None): continue
              charge = charges_class(raw_records = pdb_str).get_total_charge()
              fo.write("REMARK charge= %s \n"%str(charge))
              fo.write(pdb_str)
              fo.close()
              cntr+=1
              libtbx.easy_run.fully_buffered("mv %s %s"%(cl_fn, prefix))
        libtbx.easy_run.fully_buffered("rm -rf *.pdb ase/")

if __name__ == '__main__':
  if 1:
    result = run(file_name = "/home/yanting/qr-work/5nfm.pdb")
    print result
  else:
    result = run(file_name = "/home/pdb/pdb/p4/pdb3p4j.ent.gz")
