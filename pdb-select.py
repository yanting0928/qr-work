import os
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
from qrefine.charges import charges_class
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

def clusters(pdb_hierarchy,filename):
  yoink_utils.write_yoink_infiles("cluster.xml",
                                  "qmmm.xml",
                                  pdb_hierarchy,
                                  os.path.join(qrefine,"plugin","yoink","dat"))
  pyoink=PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"),
                os.path.join(qrefine,"plugin","yoink","dat"),
                "cluster.xml")
  backbone_connections = fragment_utils.get_backbone_connections(
                         pdb_hierarchy)
  interaction_list,weight = pyoink.get_interactions_list()
  interaction_list += backbone_connections
  cc=clustering.betweenness_centrality_clustering(
                            interaction_list,
                            maxnum_residues_in_cluster=3)
  clusters_contain=cc.get_clusters()

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
  filename=os.path.basename(file_name)[3:7]
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  resolution = get_resolution(pdb_inp = pdb_inp)
  data_type = pdb_inp.get_experiment_type()
  if resolution <= 0.9: ##set < 0.9 
    if data_type=="X-RAY DIFFRACTION" or  data_type=="NEUTRON DIFFRACTION":
      print data_type
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      print filename
      if not_only_protein(pdb_hierarchy=pdb_hierarchy):
        print "PDB file not only protein"
        pdb_hierarchy  = remove_rna_dna(pdb_hierarchy=pdb_hierarchy)
      pdb_hierarchy = add_hydrogens_using_ReadySet(pdb_hierarchy=pdb_hierarchy)
      if have_conformers(pdb_hierarchy=pdb_hierarchy):
        pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=True)
      file_name = filename+"_update"
      pdb_hierarchy.write_pdb_file(file_name = file_name+".pdb")
      pdb_inp = iotbx.pdb.input(file_name = file_name+".pdb")
      cs = box_pdb(pdb_inp=pdb_inp,filename=file_name)
      fq = fragments(pdb_hierarchy=iotbx.pdb.input(file_name = file_name+"_box.pdb").construct_hierarchy(),
                     maxnum_residues_in_cluster=5,
                     debug=True,
                     crystal_symmetry=cs)
      print fq.clusters
      for fname in os.listdir(os.getcwd()):
        if fname.endswith("_cluster.pdb"):
          if check_missing_atom(pdb_filename = fname):
            os.remove(fname)
          else:
            print fname
            ph = completion.run(pdb_filename = fname,
                      crystal_symmetry=cs,
                      model_completion=False)
            fname = fname[:-4]+"_capping.pdb"
            charge = charges_class(pdb_filename=fname).get_total_charge()
            print charge
            os.rename(fname, fname[:-4]+'_'+str(charge)+'.pdb')
      os.mkdir(filename)
      libtbx.easy_run.fully_buffered("mv %s_*  *_cluster* %s/"%(filename,filename))

if __name__ == '__main__':
  path = "/home/yanting/pdb/pdb/"
  dpath = "/home/yanting/pdb/structure_factors/"
  of = open("".join([path,"INDEX"]),"r")
  files = ["".join([path,f]).strip() for f in of.readlines()]
  of.close()
#PDB reflection data files (list of corresponding codes)
  of = open("".join([dpath,"INDEX"]),"r")
  dfiles = [
    os.path.basename("".join([path,f]).strip())[1:5] for f in of.readlines()]
  of.close()
  for f in files:
    pdb_code = os.path.basename(f)[3:7]
    if(pdb_code in dfiles):
      try:
        run(file_name = f)
      except KeyboardInterrupt:raise
      except Exception,e:
        print "FAILED:",f
        print str(e)
        print "-"*79
