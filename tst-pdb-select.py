import os
from libtbx import easy_run
import iotbx.pdb
from iotbx.pdb import remark_2_interpretation
import mmtbx.model
from cctbx import uctbx
import libtbx
import qrefine.clustering as clustering
from qrefine.utils import yoink_utils
from qrefine.plugin.yoink.pyoink import PYoink
from qrefine.fragment import fragments
from qrefine.fragment import fragment_extracts
from qrefine.fragment import get_qm_file_name_and_pdb_hierarchy

qrefine = libtbx.env.find_in_repositories("qrefine")

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
 
 
def run(file_name):
  filename=os.path.basename(file_name)[:4]
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  resolution = get_resolution(pdb_inp = pdb_inp)
  data_type = pdb_inp.get_experiment_type()
  if resolution > 0.9: 
    print resolution
    if data_type=="X-RAY DIFFRACTION" or  data_type=="NEUTRON DIFFRACTION":
      print data_type
      easy_run.go("phenix.reduce %s > %s" %(file_name,filename+"_reduce.pdb"))
      file_name = filename+"_reduce"
      pdb_inp = iotbx.pdb.input(file_name = file_name+".pdb")
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      if not_only_protein(pdb_hierarchy=pdb_hierarchy): 
        print "PDB file not only protein"
        pdb_hierarchy  = remove_rna_dna(pdb_hierarchy=pdb_hierarchy)
      if have_conformers(pdb_hierarchy=pdb_hierarchy):
        pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=True)
      file_name = file_name+"_update"
      pdb_hierarchy.write_pdb_file(file_name = file_name+".pdb")
      pdb_inp = iotbx.pdb.input(file_name = file_name+".pdb")
      cs = box_pdb(pdb_inp=pdb_inp,filename=file_name)
      fq = fragments(pdb_hierarchy= iotbx.pdb.input(file_name = file_name+"_box.pdb").construct_hierarchy(),
                     working_folder = filename,
                     maxnum_residues_in_cluster=5,
                     debug=True,
                     crystal_symmetry=cs)
      print fq.clusters
      fq_ext = fragment_extracts(fq)
      for i in xrange(len(fq.clusters)):
        get_qm_file_name_and_pdb_hierarchy(
                          fragment_extracts=fq_ext,
                          index=i)
      libtbx.easy_run.fully_buffered("mv %s_* %s/"%(filename,filename))
      libtbx.easy_run.fully_buffered("rm *.pdb")
      
if __name__ == '__main__':
  result = run(file_name = "/home/yanting/QR/ANI/qr-work/1bdw.pdb")
