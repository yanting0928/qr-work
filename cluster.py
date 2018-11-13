import os
import iotbx.pdb
from iotbx.pdb import remark_2_interpretation
import mmtbx.model
from cctbx import uctbx
import libtbx
import qrefine.clustering as clustering
from qrefine.utils import yoink_utils,fragment_utils
from qrefine.plugin.yoink.pyoink import PYoink
from scitbx.array_family import flex


qrefine = libtbx.env.find_in_repositories("qrefine")

def one_chain(pdb_hierarchy):
  chain_contain=[]
  for chain in pdb_hierarchy.chains():
    print chain
    chain_contain.append(chain)
  if (len(chain_contain) == 1):
    return True
  return False

def pdb_hierarchy_select(atoms_size, selection):
  selection_array = flex.bool(atoms_size, False)
  for item in selection:
    if(item<=atoms_size):
      selection_array[item-1] = True
  return selection_array


def clusters(pdb_hierarchy):
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
  cc=clustering.betweenness_centrality_clustering(interaction_list,maxnum_residues_in_cluster=3)
  clusters_contain=cc.get_clusters()
  print clusters_contain
#  for i in range(len(clusters_contain)):
#    cluster_selection = pdb_hierarchy_select(
#                        pdb_hierarchy.atoms_size(), clusters_contain[i])
#    print cluster_selection
#    cluster_pdb_hierarchy = pdb_hierarchy.select(cluster_selection)
#    cluster_pdb_hierarchy.write_pdb_file(file_name=str(i)+"_cluster.pdb")
  get_class = iotbx.pdb.common_residue_names_get_class
  asc = pdb_hierarchy.atom_selection_cache()
  for i in range(len(clusters_contain)):
    selection=" or ".join(" resseq %s"%i for i in clusters_contain[i])
    sel = asc.selection(selection)
    hierarchy_new = pdb_hierarchy.select(sel)
    hierarchy_new.write_pdb_file(file_name=str(i)+"_cluster.pdb")  
  #return  cc.get_clusters()

def run(file_name):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  if one_chain(pdb_hierarchy=pdb_hierarchy):
    print "Done"
    clusters(pdb_hierarchy=pdb_hierarchy)

if __name__ == '__main__':
  result = run(file_name = "/home/yanting/QR/ANI/qr-work/1yjp.pdb")
