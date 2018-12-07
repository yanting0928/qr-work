import iotbx.pdb
import libtbx
from qrefine.fragment import fragments
from scitbx.array_family import flex

qrefine = libtbx.env.find_in_repositories("qrefine")

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


def run(file_name, 
        d_min=0.9, 
        maxnum_residues_in_cluster=5, 
        filter_non_protein=True):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy = keep_protein_only(pdb_hierarchy = pdb_hierarchy)
  pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=True)
  fq = fragments(
          pdb_hierarchy              = pdb_hierarchy,
          maxnum_residues_in_cluster = maxnum_residues_in_cluster,
          debug                      = True,
          crystal_symmetry           = pdb_inp.crystal_symmetry())
  print "clusters:", fq.clusters

if __name__ == '__main__':
 # result = run(file_name = "3ZQV.pdb")
  result = run(file_name = "1SSX.pdb")
