from PyAPM import build_block_copolymer
from pkg_resources import resource_filename
import os
import shutil

if __name__ == "__main__":
	os.mkdir('test4')
	dreiding_path = resource_filename('PyAPM', 'dreiding.lt')
	shutil.copy(dreiding_path,'test4')
	os.chdir('test4')
	smi = [["COC(=O)C(C)(*)C*","CC(*)(C*)C(=O)OC1CC2CC1C1CC(CO)CC21","CC(*)(C*)C(=O)OC12CC3CC(CC(C3)C1)C2"]]
	n_moltype = 1
	n_mol = [20]
	length = [20]
	box=[75,75,75]
	component = [3]
	n_seg = [4]
	segn = [[2,8,8,2]]
	segt = [[1,2,3,1]]
	bulid_block_copolymer(smi, n_moltype, box, n_mol, length, component, n_seg, segn, segt)