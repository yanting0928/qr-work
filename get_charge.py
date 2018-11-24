from qrefine.charges import charges_class
import os

fname = "2_cluster_capping.pdb"
charge = charges_class(pdb_filename=fname).get_total_charge()
os.rename(fname, str(charge)+'_'+fname)
print charge
