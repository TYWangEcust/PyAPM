from .MoleculeBuilder import MoleculeBuilder as mb
from .ChainBuilder import ChainBuilder as cb
from .AtmosphereBuilder import AtmosphereBuilder as ab
from .CrosslinkBuilder import CrosslinkBuilder as clb
from . import CGBuilder as cgb
from . import Crosslink_lib as lib
from . import ibi
import numpy as np
import os
import shutil

def build_oligomer(smi, n_moltype, n_mol, box):
	#Build Oligomers
	#smi: SMILES of oligomers, array length = n_moltype
	#n_moltype: number of molecule types, integer
	#n_mol: number of molecules of each type, array length = n_moltype, each element is an integer
	#box: system boundary lengths, array length = 3, each element is an integer
	for i in range(n_moltype):
		mol = mb(smi[i], i+1)
		mol.optchain()#Obtain single-chain structure
		mol.get_moltemplate_monomer_lt()#Obtain Moltemplate input files
	system = ab(box, n_mol, n_moltype)
	system.write_system()#Obtain input files for Moltemplate and PACKMOL
	system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
	system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
	system.run(['cleanup_moltemplate.sh'], 'system.in.settings', timeout=3600)
	#Replace LJ interaction parameters
	system.change_lj_coeff(file1 = 'system.in.init', file2 = 'system.data', file3 = 'system.in.settings')
	return 0

def build_homopolymer(smi, n_moltype, n_mol, box, length):
	#Build homopolymers
	#smi: SMILES of the repeating unit of the polymer, array length = n_moltype
	#n_moltype: number of molecule types, integer
	#n_mol: number of molecules of each type, array length = n_moltype, each element is an integer
	#box: system boundary lengths, array length = 3, each element is an integer
	#length: chain length of each polymer type, array length = n_moltype, each element is an integer
	for i in range(n_moltype):
		chain = cb(smi[i], length[i], 1)
		polymer_smi = chain.get_homopolymer_structure()#Convert monomer SMILES to polymer SMILES
		mol = mb(polymer_smi, i+1)
		mol.optchain()#Obtain single-chain structure
		mol.get_moltemplate_monomer_lt()#Obtain Moltemplate input files
	system = ab(box, n_mol, n_moltype)
	system.write_system()#Obtain input files for Moltemplate and PACKMOL
	system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
	system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
	system.run(['cleanup_moltemplate.sh'], 'system.in.settings', timeout=3600)
	#Replace LJ interaction parameters
	system.change_lj_coeff(file1 = 'system.in.init', file2 = 'system.data', file3 = 'system.in.settings')
	return 0

def build_random_copolymer(smi, n_moltype, box, n_mol, length, component, q):
	#Build random copolymers
	#smi: SMILES of polymer repeating units, 2D list with dimensions (n_moltype, component)
	#n_moltype: number of molecule types, integer
	#n_mol: number of molecules of each type, array length = n_moltype, each element is an integer
	#box: system boundary lengths, array length = 3, each element is an integer
	#length: chain length of each polymer type, array length = n_moltype, each element is an integer
	#component: number of components on each polymer chain, array length = n_moltype, each element is an integer
	#q: proportion of each monomer type, array length = n_moltype, each element is a float between 0 and 1, sum to 1
	for i in range(n_moltype):
		chain = cb(smi[i], length[i], component[i])
		co_array = chain.get_random_array(q[i])#Build copolymer template
		polymer_smi = chain.get_copolymer_structure(co_array)#Obtain polymer SMILES
		mol = mb(polymer_smi, i+1)
		mol.optchain()#Obtain single-chain structure
		mol.get_moltemplate_monomer_lt()#Obtain Moltemplate input files
	system = ab(box, n_mol, n_moltype)
	system.write_system()#Obtain input files for Moltemplate and PACKMOL
	system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
	system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
	system.run(['cleanup_moltemplate.sh'], 'system.in.settings', timeout=3600)
	#Replace LJ interaction parameters
	system.change_lj_coeff(file1 = 'system.in.init', file2 = 'system.data', file3 = 'system.in.settings')
	return 0

def build_alternate_copolymer(smi, n_moltype, box, n_mol, length):
	#Build alternating copolymers
	#smi: SMILES of polymer repeating units, 2D array with dimensions (n_moltype, 2), each chain has two components
	#n_moltype: number of molecule types, integer
	#n_mol: number of molecules of each type, array length = n_moltype, each element is an integer
	#box: system boundary lengths, array length = 3, each element is an integer
	#length: chain length of each polymer type, array length = n_moltype, each element is an integer
	for i in range(n_moltype):
		chain = cb(smi[i], length[i], 2)
		co_array = chain.get_alternate_array()#Build copolymer template
		polymer_smi = chain.get_copolymer_structure(co_array)#Obtain polymer SMILES
		mol = mb(polymer_smi, i+1)
		mol.optchain()#Obtain single-chain structure
		mol.get_moltemplate_monomer_lt()#Obtain Moltemplate input files
	system = ab(box, n_mol, n_moltype)
	system.write_system()#Obtain input files for Moltemplate and PACKMOL
	system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
	system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
	system.run(['cleanup_moltemplate.sh'], 'system.in.settings', timeout=3600)
	#Replace LJ interaction parameters
	system.change_lj_coeff(file1 = 'system.in.init', file2 = 'system.data', file3 = 'system.in.settings')
	return 0

def build_block_copolymer(smi, n_moltype, box, n_mol, length, component, n_seg, segn, segt):
	#Build block copolymers
	#smi: SMILES of polymer repeating units, 2D array with dimensions (n_moltype, component)
	#n_moltype: number of molecule types, integer
	#n_mol: number of molecules of each type, array length = n_moltype, each element is an integer
	#box: system boundary lengths, array length = 3, each element is an integer
	#length: chain length of each polymer type, array length = n_moltype, each element is an integer
	#component: number of components on each polymer chain, array length = n_moltype, each element is an integer
	#n_seg: number of blocks on each chain, array length = n_moltype
	#segn: length of each block, 2D array with dimensions (n_moltype, n_seg)
	#segt: repeating unit type of each block, 2D array with dimensions (n_moltype, n_seg), counting from 1
	for i in range(n_moltype):
		chain = cb(smi[i], length[i], component[i])
		co_array = chain.get_block_array(n_seg[i], segn[i], segt[i])#Build copolymer template
		polymer_smi = chain.get_copolymer_structure(co_array)#Obtain polymer SMILES
		mol = mb(polymer_smi, i+1)
		mol.optchain()#Obtain single-chain structure
		mol.get_moltemplate_monomer_lt()#Obtain Moltemplate input files
	system = ab(box, n_mol, n_moltype)
	system.write_system()#Obtain input files for Moltemplate and PACKMOL
	system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
	system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
	system.run(['cleanup_moltemplate.sh'], 'system.in.settings', timeout=3600)
	#Replace LJ interaction parameters
	system.change_lj_coeff(file1 = 'system.in.init', file2 = 'system.data', file3 = 'system.in.settings')
	return 0

def build_sequential_copolymer(smi, n_moltype, box, n_mol, length, component, n_seq, seqt):
	#Build sequence block copolymers
	#smi: SMILES of polymer repeating units, 2D array with dimensions (n_moltype, component)
	#n_moltype: number of molecule types, integer
	#n_mol: number of molecules of each type, array length = n_moltype, each element is an integer
	#box: system boundary lengths, array length = 3, each element is an integer
	#length: chain length of each polymer type, array length = n_moltype, each element is an integer
	#component: number of components on each polymer chain, array length = n_moltype, each element is an integer
	#n_seg: number of repeating blocks on each chain, array length = n_moltype
	#segt: repeating unit type within each block, 2D array with dimensions (n_moltype, block_length)
	for i in range(n_moltype):
		chain = cb(smi[i], length[i], component[i])
		co_array = chain.get_sequential_array(n_seq[i], seqt[i])#Build copolymer template
		polymer_smi = chain.get_copolymer_structure(co_array)#Obtain polymer SMILES
		mol = mb(polymer_smi, i+1)
		mol.optchain()#Obtain single-chain structure
		mol.get_moltemplate_monomer_lt()#Obtain Moltemplate input files
	system = ab(box, n_mol, n_moltype)
	system.write_system()#Obtain input files for Moltemplate and PACKMOL
	system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
	system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
	system.run(['cleanup_moltemplate.sh'], 'system.in.settings', timeout=3600)
	#Replace LJ interaction parameters
	system.change_lj_coeff(file1 = 'system.in.init', file2 = 'system.data', file3 = 'system.in.settings')
	return 0

def build_branch_polymer(smi, n_moltype, box, n_mol, length, component, n_seg, segn, segt, branch_point, n_branch, t_branch):
	#Build branched copolymers
	#smi: SMILES of polymer repeating units, 2D array with dimensions (n_moltype, component)
	#n_moltype: number of molecule types, integer
	#n_mol: number of molecules of each type, array length = n_moltype, each element is an integer
	#box: system boundary lengths, array length = 3, each element is an integer
	#length: backbone length of each polymer type, array length = n_moltype, each element is an integer
	#component: number of components on each polymer chain
	#segn: length of each segment on the backbone, 2D array with dimensions (n_moltype, n_seg)
	#segt: repeating unit type of each segment on the backbone, 2D array with dimensions (n_moltype, n_seg), counting from 1
	#branch_point: repeating unit type serving as branch point; currently only one branch point type per chain is supported, 
	#array length = n_moltype, each element is an integer
	#n_branch: branch chain length at each branch point; currently no further grafting on branch chains is allowed, 
	#2D array with dimensions (n_moltype, number_of_branch_points)
	#t_branch: type of branch chain at each branch point; currently each branch chain is a homopolymer chain, 
	#2D array with dimensions (n_moltype, number_of_branch_points)
	for i in range(n_moltype):
		chain = cb(smi[i], length[i], component[i])
		#Obtain polymer SMILES
		polymer_smi = chain.get_branch_structure(n_seg[i], segn[i], segt[i], branch_point[i], n_branch[i], t_branch[i])
		mol = mb(polymer_smi, i+1)
		mol.optchain()#Obtain single-chain structure
		mol.get_moltemplate_monomer_lt()#Obtain Moltemplate input files
	system = ab(box, n_mol, n_moltype)
	system.write_system()#Obtain input files for Moltemplate and PACKMOL
	system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
	system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
	system.run(['cleanup_moltemplate.sh'], 'system.in.settings', timeout=3600)
	#Replace LJ interaction parameters
	system.change_lj_coeff(file1 = 'system.in.init', file2 = 'system.data', file3 = 'system.in.settings')
	return 0 

def build_crosslink_polymer(smi_pre, smi_post, bond_atom_style, n_react, react_temp, box, n_mol):
	#Build crosslinked polymers
	#smi_pre: SMILES of components before crosslinking reaction, 
	#2D array with dimensions (n_react, 2), the first is the main reactant
	#smi_post: SMILES of dimers after one-step crosslinking reaction, 
	#array length = n_react, the first is the main reaction product
	#Type of chemical bond connecting the backbone in the crosslinking reaction, formatted as ["N","C"]
	#n_react: total number of reactions present in the system
	#n_mol: number of molecules of each type, array length = n_moltype, each element is an integer
	#box: system boundary lengths, array length = 3, each element is an integer

	#For each reaction, build single-molecule models before and after the reaction
	for i in range(n_react):
		if not os.path.exists('reaction_%d'%i):
			polymer = clb(smi_pre[i], smi_post[i], bond_atom_style[i], i)
			polymer.prereaction()
			polymer.postreaction()
			polymer.reaction()
	#Build reaction templates
	try:
		os.mkdir('reaction')
	except FileExistsError:  
		print("Directory reaction already exists!")
	shutil.copy('dreiding.lt','reaction')
	for filename in os.listdir('AutoMapper'):
		file_path = os.path.join('AutoMapper', filename)
		if os.path.isfile(file_path):
			shutil.copy(file_path, 'reaction')

	#Build amorphous system of reactants
	shutil.copy('prereact_0/system.data', 'reaction/pre_reaction.data')
	shutil.copy('postreact_0/system.data', 'reaction/post_reaction.data')
	os.chdir('reaction')
	n_reactant = len(smi_pre[0])
	for i in range(n_reactant):
		pre_reactant = mb(smi_pre[0][i], i+1)
		pre_reactant.optchain()
		pre_reactant.get_moltemplate_monomer_lt()
	system = ab(box, n_mol, n_reactant)
	system.write_system()
	system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
	system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
	system.run(['python', 'AutoMapper.py', '.', 'clean', 'system.data', 'post_reaction.data', '--coeff_file', 'system.in.settings'], 'cleanedsystem.data', timeout=3600)
	os.chdir('..')

	#Simulate crosslinking in LAMMPS
	try:
		os.mkdir('crosslink')
	except FileExistsError:  
		print("Directory crosslink already exists!")
	shutil.copy('reaction/cleanedsystem.data','crosslink')
	shutil.copy('reaction/system.in.init','crosslink')
	shutil.copy('reaction/cleanedsystem.in.settings','crosslink')
	for i in range(n_react):
		shutil.copy('reaction_%d/pre-molecule.data'%i,'crosslink/pre-molecule%d.data'%i)
		shutil.copy('reaction_%d/post-molecule.data'%i,'crosslink/post-molecule%d.data'%i)
		shutil.copy('reaction_%d/automap.data'%i,'crosslink/automap%d.data'%i)
	os.chdir('crosslink')
	lib.write_lammps_crosslink(n_react, react_temp)
	#Replace LJ interaction parameters
	system.change_lj_coeff(file1 = 'system.in.init', file2 = 'cleanedsystem.data', file3 = 'cleanedsystem.in.settings')
	command = ['mpirun', '-np', '32', 'lmp_2408', '-in', 'in.md']
	system.run(command, 'crosslink.data', timeout=3600)#Run lammps
	os.chdir('..')
	return 0

def build_cg_block_copolymer(n_mol, n_seg, segn, segt, lx, ly, lz, n_branch, branchtype, temp, n_step, path):
	bondtype, angletype, bondcoeff, anglecoeff, mass = ibi.optimization_IBI(temp, n_step, path)
	dihedraltype = [[1,1,1,1,1]]
	chain = cgb.Branched_chain(n_mol, n_seg, segn, segt, 
		lx, ly, lz, mass, bondtype, angletype, dihedraltype, 
		bondcoeff, anglecoeff, n_branch, branchtype)
	chain.get_atom()
	chain.get_bond()
	chain.get_angle()
	chain.get_dihedral()
	chain.period_box()
	chain.write_data()
	chain.write_pdb()
	return 0
