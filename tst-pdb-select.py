import os
from libtbx import easy_run
import iotbx.pdb
from iotbx.pdb import remark_2_interpretation
import mmtbx.model
from cctbx import uctbx
import libtbx
from mmtbx.building import extend_sidechains
from mmtbx.monomer_library import server
import qrefine.clustering as clustering
from qrefine.utils import yoink_utils
from qrefine.plugin.yoink.pyoink import PYoink
from qrefine.fragment import fragments
import qrefine.completion as completion
from mmtbx.monomer_library import server
from mmtbx.monomer_library import idealized_aa
from mmtbx.rotamer import rotamer_eval
from scitbx.array_family import flex


qrefine = libtbx.env.find_in_repositories("qrefine")
mon_lib_server = server.server()

def get_resolution(pdb_inp):
  resolution = None
  resolutions = iotbx.pdb.remark_2_interpretation.extract_resolution(
    pdb_inp.extract_remark_iii_records(2))
  if(resolutions is not None):
    resolution = resolutions[0]
  return resolution

def not_only_protein(pdb_hierarchy):
  get_class = iotbx.pdb.common_residue_names_get_class
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if  not (get_class(ag.resname) == "common_amino_acid" or  get_class(ag.resname) == "modified_amino_acid"):
            return True 
  return False

def have_conformers(pdb_hierarchy):
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        if residue_group.have_conformers():
          return True
  return False

def box_pdb(pdb_inp,filename):
  model = mmtbx.model.manager(model_input = pdb_inp)
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=model.get_sites_cart(),
    buffer_layer=10)
  model.set_sites_cart(box.sites_cart)
  cs  = box.crystal_symmetry()
  box_file=open(filename+"_box.pdb",'w')
  box_file.write(model.model_as_pdb())
  box_file.close()
  return cs

def remove_rna_dna(pdb_hierarchy):
  get_class = iotbx.pdb.common_residue_names_get_class
  not_protein_resname=[]
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if  not (get_class(ag.resname) == "common_amino_acid" or  get_class(ag.resname) == "modified_amino_acid"):
            not_protein_resname.append(ag.resname.strip())
  not_protein_resname_h=list(set(not_protein_resname))
  selection=" and ".join("not resname %s"%i for i in not_protein_resname_h)
  asc = pdb_hierarchy.atom_selection_cache()
  sel = asc.selection(selection)
  pdb_hierarchy = pdb_hierarchy.select(sel)

  return pdb_hierarchy 

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
 
def complete_missing_atom(pdb_hierarchy,mon_lib_server):
  n_changed = extend_sidechains.extend_protein_model(
        pdb_hierarchy = pdb_hierarchy,
        mon_lib_srv = mon_lib_server,
        add_hydrogens=False,
        )
  return n_changed,pdb_hierarchy 

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
#  pdb_lines = open(pdb_filename, 'rb').read()
  pdb_lines = pdb_hierarchy.as_pdb_string()
  output_file_name =  "reduce.pdb"
  rc = run_though_all_the_options(
    pdb_lines,
    [], # args
    hydrogens=True,
    reduce_args=['-BUILD'],
    ligands=False,
    add_h_to_water=True,
    metals=False,
    output_file_name=output_file_name+".pdb", 
    )
  return rc['model_hierarchy']  

def run(file_name):
  filename=os.path.basename(file_name)[:4]
  print filename
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  resolution = get_resolution(pdb_inp = pdb_inp)
  data_type = pdb_inp.get_experiment_type()
  if resolution > 0.9: 
    print resolution
    if data_type=="X-RAY DIFFRACTION" or  data_type=="NEUTRON DIFFRACTION":
      print data_type
      pdb_hierarchy = pdb_inp.construct_hierarchy()
#      n_changed,pdb_hierarchy = complete_missing_atom(pdb_hierarchy = pdb_hierarchy,mon_lib_server = mon_lib_server)
#      if n_changed:
#        pdb_hierarchy.write_pdb_file(file_name = filename+"_complete.pdb")
#        file_name = filename+"_complete"
      if not_only_protein(pdb_hierarchy=pdb_hierarchy): 
        print "PDB file not only protein"
        pdb_hierarchy  = remove_rna_dna(pdb_hierarchy=pdb_hierarchy)
      pdb_hierarchy = add_hydrogens_using_ReadySet(pdb_hierarchy=pdb_hierarchy)
      if have_conformers(pdb_hierarchy=pdb_hierarchy):
        pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=True)
      file_name = file_name+"_update"
      pdb_hierarchy.write_pdb_file(file_name = file_name+".pdb")
      pdb_inp = iotbx.pdb.input(file_name = file_name+".pdb")
      cs = box_pdb(pdb_inp = pdb_inp,filename=file_name)
      fq = fragments(pdb_hierarchy= iotbx.pdb.input(file_name = file_name+"_box.pdb").construct_hierarchy(),
                     maxnum_residues_in_cluster=5,
                     debug=True,
                     crystal_symmetry=cs)
      print fq.clusters
      for fname in os.listdir(os.getcwd()):
        if fname.endswith("_cluster.pdb"):
          if check_missing_atom(pdb_filename = fname):
            os.remove(fname)
          else:
            ph = completion.run(pdb_filename = os.path.join(os.getcwd(),fname),
                      crystal_symmetry=cs,
                      model_completion=False)
      os.mkdir(filename)
      libtbx.easy_run.fully_buffered("mv %s_*  *_cluster* %s/"%(filename,filename))
      libtbx.easy_run.fully_buffered("rm -rf *.pdb ase/")
      
if __name__ == '__main__':
  result = run(file_name = "/home/yanting/QR/ANI/qr-work/1kyc.pdb")
