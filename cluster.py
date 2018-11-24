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
from qrefine.fragment import fragments,fragment_extracts
import qrefine.completion as completion
qrefine = libtbx.env.find_in_repositories("qrefine")
qr_yoink_path =os.path.join(qrefine, "plugin","yoink","yoink")

def clusters(pdb_hierarchy,cs):
  fq = fragments(
    pdb_hierarchy=pdb_hierarchy,
    crystal_symmetry=cs,
    charge_embedding=True,
    maxnum_residues_in_cluster=5,
    debug = True)# not run qm_calculation, just the clustering result
  print fq.clusters
  for filename in os.listdir(os.getcwd()):
    if filename.endswith("_cluster.pdb"):
      ph = completion.run(pdb_filename = os.path.join(os.getcwd(),filename),
                      crystal_symmetry=cs,
                      model_completion=False)

def run(file_name):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  clusters(pdb_hierarchy=pdb_hierarchy,cs=cs)

if __name__ == '__main__':
  result = run(file_name = "/home/yanting/QR/ANI/qr-work/1yjp.pdb")
