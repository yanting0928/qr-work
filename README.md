## High resolution subset of PDB

1. resolution < 0.8, X-ray or Neutron data structure
2. remove all but protein(remove all the H/D)
3. remove alternative confromations
4. add H using ReadySet
5. place into a box
6. split into fragments

   Note: any fragment contain missing atoms amino acid will be discard

7. add capping H
8. Verify that each fragment is reasonable(by bond rmsd, angle rmsd)
9. calculate charge
10. write out PDB file with fragment and charge in its REMARK

