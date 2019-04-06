import os
from libtbx.utils import null_out
from mmtbx import model
import iotbx.pdb
from cctbx import geometry_restraints
import mmtbx
import libtbx
import numpy as np

def run(file_name):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  model = mmtbx.model.manager(
    model_input       = pdb_inp,
    build_grm         = True,
    stop_for_unknowns = False,
    log               = null_out())
  grm = model.get_restraints_manager().geometry
  sites_cart = model.get_xray_structure().sites_cart()
  delta_list = []
  pair_proxies = grm.pair_proxies(
    sites_cart=sites_cart)
  for proxy in pair_proxies.bond_proxies.simple:
    bond_d = geometry_restraints.bond(
      sites_cart   = sites_cart,
      proxy        = proxy)
    delta_list.append(bond_d.delta)
  sigma = np.std(np.array(delta_list))
  for i in delta_list:
    if abs(i)>4*sigma:
      print file_name
      libtbx.easy_run.fully_buffered("rm %s"%file_name)
  


if __name__ == '__main__':
  if 1:
    path = "/home/yanting/qr-work/work-dir/"
    for subdir, dirs, files in os.walk(path):
      for pdb_file in files:
        if pdb_file.endswith(".pdb"):
          run(os.path.join(subdir, pdb_file))
  else:
#    run("/home/yanting/qr-work/work-dir/1us0.pdb")  
    run("/home/yanting/qr-work/3X2M_31.pdb")
