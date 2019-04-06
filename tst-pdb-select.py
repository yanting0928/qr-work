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

qrefine = libtbx.env.find_in_repositories("qrefine")
mon_lib_server = server.server()

def keep_protein_only(pdb_hierarchy):
  selection = flex.size_t()
  get_class = iotbx.pdb.common_residue_names_get_class
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if(get_class(ag.resname) == "common_amino_acid"):
            selection.extend(ag.atoms().extract_i_seq())
  if(selection.size()==0): return None
  else:                    return pdb_hierarchy.select(selection)

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

def run(file_name, 
        d_min=0.9, 
        maxnum_residues_in_cluster=3, 
        filter_non_protein=True):
  prefix = os.path.basename(file_name)[:4]
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
      if(filter_non_protein):
        pdb_hierarchy = keep_protein_only(pdb_hierarchy = pdb_hierarchy)
      if(pdb_hierarchy is not None):
        if(have_conformers(pdb_hierarchy=pdb_hierarchy)):
          pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=True)
        # Remove ori H, using ReasySet add H
        hd_selection = pdb_hierarchy.atom_selection_cache().selection("element H")
        pdb_hierarchy = pdb_hierarchy.select(~hd_selection)
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
              charge = charges_class(raw_records = pdb_str).get_total_charge()
              fo.write("REMARK charge= %s \n"%str(charge))
              fo.write(pdb_str)
              fo.close()
              cntr+=1
              libtbx.easy_run.fully_buffered("mv %s %s"%(cl_fn, prefix))
        libtbx.easy_run.fully_buffered("rm -rf *.pdb ase/")

if __name__ == '__main__':
  result = run(file_name = "/home/yanting/QR/ANI/qr-work/4oy5.pdb")
