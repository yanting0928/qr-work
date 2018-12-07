import iotbx.pdb
from qrefine.fragment import fragments


def run(file_name, 
        d_min=0.9, 
        maxnum_residues_in_cluster=5, 
        filter_non_protein=True):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  fq = fragments(
          pdb_hierarchy              = pdb_hierarchy,
          maxnum_residues_in_cluster = maxnum_residues_in_cluster,
          debug                      = True,
          crystal_symmetry           = pdb_inp.crystal_symmetry())
  print "clusters:", fq.clusters

if __name__ == '__main__':
  result = run(file_name = "1JXT.pdb")
