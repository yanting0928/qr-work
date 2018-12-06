## High resolution subset of PDB

1. resolution < 0.9, X-ray or Neutron data structure
2. remove all but protein
3. remove alternative confromations
4. add H
5. place into a box
6. split into fragments

   Note: any fragment contain missing atoms amino acid will be discard

7. add capping H
8. calculate charge
9. write out PDB file with fragment and charge in its REMARK

