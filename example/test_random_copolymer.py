from PyAPM import build_random_copolymer
from pkg_resources import resource_filename
import os
import shutil

if __name__ == "__main__":
	os.mkdir('test2')
	dreiding_path = resource_filename('PyAPM', 'dreiding.lt')
	shutil.copy(dreiding_path,'test2')
	os.chdir('test2')
	smi = [["COC(=O)C(C)(*)C*","CC(*)(C*)C(=O)OC1CC2CC1C1CC(CO)CC21"],["CC(*)(C*)C(=O)OC12CC3CC(CC(C3)C1)C2","COC(=O)C(C)(*)C*"]]
	n_moltype = 2
	n_mol = [20,30]
	length = [30,20]
	box=[70,70,70]
	component = [2,2]
	q = [[0.4,0.6],[0.7,0.3]]
	build_random_copolymer(smi, n_moltype, box, n_mol, length, component, q)