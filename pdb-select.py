import os
import iotbx.pdb
from iotbx.pdb import remark_2_interpretation

def get_resolution(pdb_inp):
  resolution = None
  resolutions = iotbx.pdb.remark_2_interpretation.extract_resolution(
    pdb_inp.extract_remark_iii_records(2))
  if(resolutions is not None):
    resolution = resolutions[0]
  return resolution
  
def run(file_name):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  resolution = get_resolution(pdb_inp = pdb_inp)
  num = 0
  i = 1
  if resolution < 0.9:
    num += 1
    i += 1
    print resolution
  print num 

path = "/home/yanting/pdb/pdb/"
dpath = "/home/yanting/pdb/structure_factors/"
of = open("".join([path,"INDEX"]),"r")
files = ["".join([path,f]).strip() for f in of.readlines()]
of.close()
#PDB reflection data files (list of corresponding codes)
of = open("".join([dpath,"INDEX"]),"r")
dfiles = [
  os.path.basename("".join([path,f]).strip())[1:5] for f in of.readlines()]
of.close()
num = 0
for f in files:
  pdb_code = os.path.basename(f)[3:7]
  if(pdb_code in dfiles):
    #try:
    pdb_inp = iotbx.pdb.input(file_name = f)
    resolution = get_resolution(pdb_inp = pdb_inp)
    if resolution < 0.9:
      num += 1
      print resolution
print num
   # except KeyboardInterrupt:raise
   # except Exception,e:
   #   print "FAILED:",f
   #   print str(e)
   #   print "-"*79
