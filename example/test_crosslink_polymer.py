from PyAPM import build_crosslink_polymer
from pkg_resources import resource_filename
import os
import shutil

if __name__ == "__main__":
	if not os.path.exists('test9'):
		os.mkdir('test9')
		dreiding_path = resource_filename('PyAPM', 'dreiding.lt')
		shutil.copy(dreiding_path,'test9')
		automapper_path = resource_filename('PyAPM', 'AutoMapper')
		shutil.copytree(automapper_path, 'test9/AutoMapper')
	os.chdir('test9')
	smi_pre_init = ["NC1=CC=CC(S(=O)(=O)C2=CC=CC(N)=C2)=C1","C1=C(CC2=CC=C(N(CC3CO3)CC3CO3)C=C2)C=CC(N(CC2CO2)CC2CO2)=C1"]
	smi_post_init = "NC1=CC=CC(S(=O)(=O)C2=CC=CC(NCC(O)CN(CC3CO3)C3=CC=C(CC4=CC=C(N(CC5CO5)CC5CO5)C=C4)C=C3)=C2)=C1"
	bond_atom_style_init = ["N","C"]
	#副反应1
	smi_pre_1 = ["NC1=CC=CC(S(=O)(=O)C2=CC=CC(NCC(O)CN(CC3CO3)C3=CC=C(CC4=CC=C(N(CC5CO5)CC5CO5)C=C4)C=C3)=C2)=C1","C1=C(CC2=CC=C(N(CC3CO3)CC3CO3)C=C2)C=CC(N(CC2CO2)CC2CO2)=C1"]
	smi_post_1 = "NC1=CC=CC(S(=O)(=O)C2=CC(NCC(CN(CC3CO3)C3=CC=C(CC4=CC=C(N(CC5CO5)CC5CO5)C=C4)C=C3)OCC(O)CN(CC3CO3)C3=CC=C(CC4=CC=C(N(CC5CO5)CC5CO5)C=C4)C=C3)=CC=C2)=C1"
	smi_pre_2 = ["NC1=CC=CC(S(=O)(=O)C2=CC=CC(NCC(O)CN(CC3CO3)C3=CC=C(CC4=CC=C(N(CC5CO5)CC5CO5)C=C4)C=C3)=C2)=C1", "C1=C(CC2=CC=C(N(CC3CO3)CC3CO3)C=C2)C=CC(N(CC2CO2)CC2CO2)=C1"]
	smi_post_2 = "NC1=CC(S(=O)(=O)C2=CC=CC(N(CC(O)CN(CC3CO3)C3=CC=C(CC4=CC=C(N(CC5CO5)CC5CO5)C=C4)C=C3)CC(O)CN(CC3CO3)C3=CC=C(CC4=CC=C(N(CC5CO5)CC5CO5)C=C4)C=C3)=C2)=CC=C1"
	bond_atom_style_1 = ["O","C"]
	bond_atom_style_2 = ["O","C"]
	smi_pre = [smi_pre_init, smi_pre_1, smi_pre_2]
	smi_post = [smi_post_init, smi_post_1, smi_post_2]
	bond_atom_style = [bond_atom_style_init, bond_atom_style_1, bond_atom_style_2]
	n_react = 3
	box = [200,200,200]
	n_mol = [500,250]
	react_temp = 473
	build_crosslink_polymer(smi_pre, smi_post, bond_atom_style, n_react, react_temp, box, n_mol)