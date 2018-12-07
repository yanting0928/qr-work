import iotbx.pdb

def add_hydrogens_using_ReadySet(pdb_hierarchy):
  from elbow.command_line.ready_set import run_though_all_the_options
  pdb_lines = pdb_hierarchy.as_pdb_string()
  rc = run_though_all_the_options(
    pdb_lines,
    [], # args
    hydrogens=True,
    reduce_args=['-BUILD', '-NUClear'],
    ligands=False,
    add_h_to_water=True,
    metals=False,
    output_file_name="reasyset.pdb",
    verbose=False,
    silent=True,
    )
  return rc['model_hierarchy']

def run(file_name, 
        d_min=0.9, 
        maxnum_residues_in_cluster=5, 
        filter_non_protein=True):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy = add_hydrogens_using_ReadySet(
                  pdb_hierarchy=pdb_hierarchy)

if __name__ == '__main__':
  result = run(file_name = "4PTH.pdb")
