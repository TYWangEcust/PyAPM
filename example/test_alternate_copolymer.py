from PyAPM import build_alternate_copolymer
from pkg_resources import resource_filename
import os
import shutil

if __name__ == "__main__":
	os.mkdir('test3')
	dreiding_path = resource_filename('PyAPM', 'dreiding.lt')
	shutil.copy(dreiding_path,'test3')
	os.chdir('test3')
	smi = [["COC(=O)C(C)(*)C*","CC(*)(C*)C(=O)OC1CC2CC1C1CC(CO)CC21"],["CC(*)(C*)C(=O)OC12CC3CC(CC(C3)C1)C2","COC(=O)C(C)(*)C*"]]
	n_moltype = 2
	n_mol = [20,30]
	length = [30,20]
	box=[120,120,120]
	build_alternate_copolymer(smi, n_moltype, box, n_mol, length)