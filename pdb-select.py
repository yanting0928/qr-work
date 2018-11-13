import os
import iotbx.pdb
from iotbx.pdb import remark_2_interpretation
import mmtbx.model
from cctbx import uctbx
import libtbx
import qrefine.clustering as clustering
from qrefine.utils import yoink_utils,fragment_utils
from qrefine.plugin.yoink.pyoink import PYoink
from qrefine.fragment import fragments
from qrefine.fragment import fragment_extracts
from qrefine.fragment import get_qm_file_name_and_pdb_hierarchy


qrefine = libtbx.env.find_in_repositories("qrefine")


def get_resolution(pdb_inp):
  resolution = None
  resolutions = iotbx.pdb.remark_2_interpretation.extract_resolution(
    pdb_inp.extract_remark_iii_records(2))
  if(resolutions is not None):
    resolution = resolutions[0]
  return resolution

def not_only_protein(pdb_hierarchy):
  get_class = iotbx.pdb.common_residue_names_get_class
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if get_class(ag.resname) == "common_rna_dna":
            return True
  return False

def have_conformers(pdb_hierarchy):
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        if residue_group.have_conformers():
          return True
  return False

def box_pdb(pdb_inp,filename):
  model = mmtbx.model.manager(model_input = pdb_inp)
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=model.get_sites_cart(),
    buffer_layer=10)
  model.set_sites_cart(box.sites_cart)
  cs  = box.crystal_symmetry()
  box_file=open(filename[3:7]+"_box.pdb",'w')
  box_file.write(model.model_as_pdb())
  box_file.close()
  return cs

def remove_rna_dna(pdb_hierarchy):
  get_class = iotbx.pdb.common_residue_names_get_class
  not_protein_resname=[]
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if get_class(ag.resname) == "common_rna_dna":
            not_protein_resname.append(ag.resname.strip())
  return not_protein_resname

def clusters(pdb_hierarchy,filename):
  yoink_utils.write_yoink_infiles("cluster.xml",
                                  "qmmm.xml",
                                  pdb_hierarchy,
                                  os.path.join(qrefine,"plugin","yoink","dat"))
  pyoink=PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"),
                os.path.join(qrefine,"plugin","yoink","dat"),
                "cluster.xml")
  backbone_connections = fragment_utils.get_backbone_connections(
                         pdb_hierarchy)
  interaction_list,weight = pyoink.get_interactions_list()
  interaction_list += backbone_connections
  cc=clustering.betweenness_centrality_clustering(
                            interaction_list,
                            maxnum_residues_in_cluster=3)
  clusters_contain=cc.get_clusters()
  print clusters_contain
  get_class = iotbx.pdb.common_residue_names_get_class
  asc = pdb_hierarchy.atom_selection_cache()
  os.mkdir(filename[3:7])
  working_fold = os.path.abspath(filename[3:7])
  print os.getcwd()
  os.chdir(working_fold)
  print os.getcwd()
  for i in range(len(clusters_contain)):
    selection=" or ".join(" resseq %s"%i for i in clusters_contain[i])
    sel = asc.selection(selection)
    hierarchy_new = pdb_hierarchy.select(sel)
    hierarchy_new.write_pdb_file(file_name=str(i)+"_cluster.pdb")

def run(file_name):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  resolution = get_resolution(pdb_inp = pdb_inp)
  data_type = pdb_inp.get_experiment_type()
  filename=os.path.basename(file_name)
  asc = pdb_hierarchy.atom_selection_cache()
  if resolution < 0.9: ##set < 0.9 
    if data_type=="X-RAY DIFFRACTION" or  data_type=="NEUTRON DIFFRACTION":
      print data_type
      os.mkdir(filename[3:7])
      working_fold = os.path.abspath(filename[3:7])
      os.chdir(working_fold)
      if not_only_protein(pdb_hierarchy=pdb_hierarchy):
        print "PDB file not only protein"
        not_protein_resname = remove_rna_dna(pdb_hierarchy=pdb_hierarchy)
        selection=" and ".join("not resname %s"%i for i in not_protein_resname)
        sel = asc.selection(selection)
        pdb_hierarchy = pdb_hierarchy.select(sel)
      if have_conformers(pdb_hierarchy=pdb_hierarchy):
        pdb_hierarchy.remove_alt_confs(always_keep_one_conformer=True)
      pdb_hierarchy.write_pdb_file(file_name = filename[:3]+"_new.pdb")
      pdb_inp = iotbx.pdb.input(file_name = filename[:3]+"_new.pdb")
      cs = box_pdb(pdb_inp=pdb_inp,filename=filename)
      fq = fragments(pdb_hierarchy=pdb_hierarchy,
                     maxnum_residues_in_cluster=5,
                     debug=True,
                     crystal_symmetry=cs)
      print fq.clusters
      fq_ext = fragment_extracts(fq)
      for i in xrange(len(fq.clusters)):
        get_qm_file_name_and_pdb_hierarchy(
                          fragment_extracts=fq_ext,
                          index=i)

if __name__ == '__main__':
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
  for f in files:
    os.chdir("/home/yanting/QR/ANI/qr-work")
    pdb_code = os.path.basename(f)[3:7]
    if(pdb_code in dfiles):
      try:
        run(file_name = f)
      except KeyboardInterrupt:raise
      except Exception,e:
        print "FAILED:",f
        print str(e)
        print "-"*79
