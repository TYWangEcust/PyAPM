from PyAPM import build_sequential_copolymer
from pkg_resources import resource_filename
import os
import shutil


if __name__ == "__main__":
	os.mkdir('test5')
	dreiding_path = resource_filename('PyAPM', 'dreiding.lt')
	shutil.copy(dreiding_path,'test5')
	os.chdir('test5')
	smi = [["COC(=O)C(C)(*)C*","CC(*)(C*)C(=O)OC1CCOC1=O","CC(*)(C*)C(=O)OC12CC3CC(CC(C3)C1)C2"]]
	n_moltype = 1
	n_mol = [20]
	length = [30]
	box=[100,100,100]
	component = [3]
	n_seq = [10]
	seqt = [[1,2,3]]
	build_sequential_copolymer(smi, n_moltype, box, n_mol, length, component, n_seq, seqt)