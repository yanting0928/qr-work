import os
import iotbx.pdb
from iotbx.pdb import remark_2_interpretation
import mmtbx.model
from cctbx import uctbx
import libtbx
import qrefine.clustering as clustering
from qrefine.utils import yoink_utils
from qrefine.plugin.yoink.pyoink import PYoink

qrefine = libtbx.env.find_in_repositories("qrefine")

def get_resolution(pdb_inp):
  resolution = None
  resolutions = iotbx.pdb.remark_2_interpretation.extract_resolution(
    pdb_inp.extract_remark_iii_records(2))
  if(resolutions is not None):
    resolution = resolutions[0]
  return resolution

def only_protein(pdb_hierarchy):
  get_class = iotbx.pdb.common_residue_names_get_class
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if get_class(ag.resname) == "common_rna_dna":
            return True 
  return False

def box_pdb(pdb_inp,filename):
  model = mmtbx.model.manager(model_input = pdb_inp)
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=model.get_sites_cart(),
    buffer_layer=10)
  model.set_sites_cart(box.sites_cart)
  model._crystal_symmetry = box.crystal_symmetry()
  print  model.model_as_pdb()  
  box_file=open(filename[:4]+"_box.pdb",'w')
  box_file.write(model.model_as_pdb())
  box_file.close()

def remove_rna_dna(pdb_hierarchy):
  get_class = iotbx.pdb.common_residue_names_get_class
  not_protein_resname=[]
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if get_class(ag.resname) == "common_rna_dna":
            not_protein_resname.append(ag.resname.strip())
  return not_protein_resname

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
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  resolution = get_resolution(pdb_inp = pdb_inp)
  data_type = pdb_inp.get_experiment_type()
  filename=os.path.basename(file_name)
  asc = pdb_hierarchy.atom_selection_cache()
  if resolution > 0.9: 
    print resolution
    if data_type=="X-RAY DIFFRACTION" or  data_type=="NEUTRON DIFFRACTION":
      print data_type
      if only_protein(pdb_hierarchy=pdb_hierarchy): 
        print "PDB file not only protein"
        not_protein_resname = remove_rna_dna(pdb_hierarchy=pdb_hierarchy)
        selection=" and ".join("not resname %s"%i for i in not_protein_resname)
        print selection
        sel = asc.selection(selection)
        hierarchy_new = pdb_hierarchy.select(sel)
        hierarchy_new.write_pdb_file(file_name="new_pdb.pdb")
        #print box_pdb(pdb_inp=pdb_inp,filename=filename)
  #clusters(pdb_hierarchy=pdb_hierarchy_new)
  #pdb_hierarchy_new.write_pdb_file(file_name="new_pdb.pdb")
if __name__ == '__main__':
  result = run(file_name = "/home/yanting/QR/ANI/6AI6.pdb")
