# Error information

## Error 1

- 1JXT, 1PQ5
```
Sorry: Number of nonbonded interaction distances < 0.001: 3
Please inspect the output above (or in the log window) and correct the input model file. 
```
## Error 2
- 1EJG, 1PJX, 2F01, 3O4P
```
  File "/home/yanting/phenix-dev-3319/modules/qrefine/fragment.py", line 28, in check_atoms_integrity  
     assert len(item) in [0,2], 'error in cluster %s %s' % (key, item)   
  AssertionError: error in cluster A 2 ['"ATOM 9 CA GLU A 2 .*. C "']
```
## Error 3
- 4HP2
``` 
  File "/home/yanting/phenix-dev-3319/modules/qrefine/completion.py", line 234, in add_c_terminal_oxygens_to_atom_group
    dihedral = dihedral_angle(sites=[atom.xyz,]
  AttributeError: 'NoneType' object has no attribute 'xyz'
```  
## Error 4 
- 3ZQV, 3ZTM, 2YL7, 2YKZ, 1X6Z, 4Y9W, 2WUR, 1US0, 2PFH, 2PEV, 4XXG, 4REK, 2I16, 2I17
```
File "/home/yanting/phenix-dev-3319/modules/qrefine/completion.py", line 804, in generate_atom_group_atom_names
    assert 0    
AssertionError
```
## Error 5
- 1SSX, 2H5D, 2H5C
```
Sorry: number of groups of duplicate atom labels: 1
  total number of affected atoms:          2
  group "ATOM    .*.  H2  VAL A 120 .*.     H  "
"ATOM .*. H2 VAL A 120 .*. H "
```
## Error 6
- 4PTH
```
  RuntimeError: Bad reduce output file: "reasyset.pdb.reduce.pdb"
```
    
 
