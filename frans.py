import os
import iotbx.pdb
import libtbx
import qrefine.clustering as clustering
from qrefine.utils import yoink_utils
from qrefine.plugin.yoink.pyoink import PYoink

if __name__ == '__main__':
  qrefine = libtbx.env.find_in_repositories("qrefine")
  pdb_inp=iotbx.pdb.input(file_name="1yjp.pdb")
  ph = pdb_inp.construct_hierarchy()
  yoink_utils.write_yoink_infiles("cluster.xml",
                                "qmmm.xml",
                                ph,
                                os.path.join(qrefine,"plugin","yoink","dat"))
  pyoink=PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"))
  pyoink.input_file="cluster.xml"
  interaction_list,weight = pyoink.get_interactions_list()
  cc=clustering.betweenness_centrality_clustering(
                          interaction_list,
                          maxnum_residues_in_cluster=3)
  clusters=cc.get_clusters()
  print clusters
