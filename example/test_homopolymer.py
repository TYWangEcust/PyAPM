from PyAPM import build_homopolymer
from pkg_resources import resource_filename
import os
import shutil

if __name__ == "__main__":
	os.mkdir('test1')
	dreiding_path = resource_filename('PyAPM', 'dreiding.lt')
	shutil.copy(dreiding_path,'test1')
	os.chdir('test1')
	smi = ["CC(*)(C*)C(=O)OC12CC3CC(CC(C3)C1)C2"]
	n_moltype = 1
	n_mol = [10]
	length = [25]
	box=[50,50,50]
	build_homopolymer(smi, n_moltype, n_mol, box, length)