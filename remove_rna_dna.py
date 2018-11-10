import os
import iotbx.pdb

pdb_inp = iotbx.pdb.input(file_name = "6AI6.pdb")
pdb_hierarchy = pdb_inp.construct_hierarchy()
get_class = iotbx.pdb.common_residue_names_get_class
asc = pdb_hierarchy.atom_selection_cache()
not_protein_resseq=[]
for model in pdb_hierarchy.models():
  for chain in model.chains():
    for rg in chain.residue_groups():
      for ag in rg.atom_groups():
        if get_class(ag.resname) == "common_rna_dna":
          not_protein_resseq.append(rg.resseq.strip())
print not_protein_resseq
#selection=" not resseq 1  and not  resseq 2 and not  resseq 3"
selection=" and ".join("not resseq %s"%i for i in not_protein_resseq)
print selection
sel = asc.selection(selection)
print 1
hierarchy_new = pdb_hierarchy.select(sel)
print 2
hierarchy_new.write_pdb_file(file_name="new_pdb.pdb")
print 3
