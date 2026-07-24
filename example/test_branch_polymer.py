from PyAPM import build_branch_polymer
from pkg_resources import resource_filename
import os
import shutil


if __name__ == "__main__":
	os.mkdir('test6')
	dreiding_path = resource_filename('PyAPM', 'dreiding.lt')
	shutil.copy(dreiding_path,'test6')
	os.chdir('test6')
	smi = [["COC(=O)C(C)(*)C*","CC(*)(C*)C(=O)OC1CC2CC1C1CC(CO)CC21","CC(*)(C*)C(=O)OC12CC3CC(CC(C3)C1)C2", "RCOC(=O)C(*)(C)C*"]]
	n_moltype = 1
	box = [50,50,50]
	n_mol = [10]
	length = [10]
	component = [4]
	n_seg = [3]
	segn = [[4,2,4]]
	segt = [[1,4,2]]
	branch_point = [4]
	n_branch = [[2,3]]
	t_branch = [[3,3]]
	build_branch_polymer(smi, n_moltype, box, n_mol, length, component, n_seg, segn, segt, branch_point, n_branch, t_branch)