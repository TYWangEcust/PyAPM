import numpy as np
import shutil
import os
from .MoleculeBuilder import MoleculeBuilder as mb
from .AtmosphereBuilder import AtmosphereBuilder as ab
from . import Crosslink_lib as lib
import subprocess

class CrosslinkBuilder:
	def __init__(self, smi_pre, smi_post, bond_atom_style, n_react):
		self.smi_pre = smi_pre #Reactant
		self.smi_post = smi_post #Productor
		self.n_react = n_react #Current reaction index being executed
		self.bond_atom_style = bond_atom_style

	def prereaction(self):
		#Build a system containing a set of reactants
		try:
			os.mkdir('prereact_%d'%self.n_react)
		except FileExistsError:  
			print("Directory prereact already exists!")  
		shutil.copy('dreiding.lt','prereact_%d'%self.n_react)
		os.chdir('prereact_%d'%self.n_react)

		system = ab([20,20,20], [1,1], 2)
		system.write_system()

		pre_reactant_1 = mb(self.smi_pre[0], 1)
		pre_reactant_1.optchain()
		pre_reactant_1.get_moltemplate_monomer_lt()
		pre_reactant_2 = mb(self.smi_pre[1], 2)
		pre_reactant_2.optchain()
		pre_reactant_2.get_moltemplate_monomer_lt()

		process_1 = system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
		if process_1:
			print("Prereaction: packmol_sub.sh successfully completed.")

		process_2 = system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
		if process_2:  
			print("Prereaction: moltemplate_sub.sh successfully completed.")
		os.chdir('..')
		return 0 

	def postreaction(self):
		#Build a system containing a set of products
		try:
			os.mkdir('postreact_%d'%self.n_react)
		except FileExistsError:  
			print("Directory postreact already exists!")  
		shutil.copy('dreiding.lt','postreact_%d'%self.n_react)
		os.chdir('postreact_%d'%self.n_react)

		system = ab([20,20,20], [1], 1)
		system.write_system()

		post_reactant = mb(self.smi_post, 1)
		post_reactant.optchain()
		post_reactant.get_moltemplate_monomer_lt()

		process_1 = system.run(['packmol', '-i', 'mixture.inp'], 'Polymer.pdb', timeout=3600)
		if process_1:
			print("Prereaction: packmol_sub.sh successfully completed.")
		process_2 = system.run(['moltemplate.sh', '-pdb', 'Polymer.pdb', 'system.lt'], 'system.in.settings', timeout=3600)
		if process_2:  
			print("Prereaction: moltemplate_sub.sh successfully completed.")
		os.chdir('..')
		return 0

	def reaction(self):
		#Use Automapper to build a reaction template
		try:
			os.mkdir('reaction_%d'%self.n_react)
		except FileExistsError:  
			print("Directory reaction already exists!")

		shutil.copy('prereact_%d/system.data'%self.n_react, 'reaction_%d/pre_reaction.data'%self.n_react)
		shutil.copy('postreact_%d/system.data'%self.n_react, 'reaction_%d/post_reaction.data'%self.n_react)
		shutil.copy('prereact_%d/system.in.settings'%self.n_react, 'reaction_%d/system.in.settings'%self.n_react)
		shutil.copy('prereact_%d/system.in.init'%self.n_react, 'reaction_%d/system.in.init'%self.n_react)
		for filename in os.listdir('AutoMapper'):
			file_path = os.path.join('AutoMapper', filename)
			if os.path.isfile(file_path):
				shutil.copy(file_path, 'reaction_%d'%self.n_react)
		os.chdir('reaction_%d'%self.n_react)
		cmd_clean = ['python', 'AutoMapper.py', '.', 'clean', 'pre_reaction.data', 'post_reaction.data', '--coeff_file', 'system.in.settings']
		result_clean = subprocess.run(cmd_clean, capture_output=True, text=True)
		if result_clean.returncode != 0:
			print("AutoMapper clean failed:", result_clean.stderr)
			return 1
		if not os.path.exists('cleanedpost_reaction.data'):
			print("cleanedpost_reaction.data not generated")
			return 1

		#Process settings and obtain the reaction atom indices in the system before and after the reaction
		with open("cleanedpost_reaction.data", 'r') as file:
			lines = file.readlines()
			tidiedLines = lib.clean_data(lines)
			sectionIndexList = lib.find_sections(tidiedLines)
			masses = lib.get_data('Masses', tidiedLines, sectionIndexList)
			potencial_element = lib.get_style(masses)
			result_string = ' '.join(potencial_element)
		element_str =[]
		for h in range(len(masses)):
			element_str.append(masses[h][3])
		n = ((len(masses)+1)*len(masses))/2
		lib.remove_first_n_lines('cleanedsystem.in.settings', n)
		
		cmd_map = ['python', 'AutoMapper.py', '.', 'newmap', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data',
			'--save_name', 'pre-molecule.data', 'post-molecule.data',
			'--ba', self.bond_atom_style[0], self.bond_atom_style[1],
			'--ebt', result_string]
		result_map = subprocess.run(cmd_map, capture_output=True, text=True)
		if result_map.returncode != 0:
			print("AutoMapper newmap failed:", result_map.stderr)
			return 1
		if not os.path.exists('post-molecule.data'):
			print("post-molecule.data not generated")
			return 1
		os.chdir('..')
		return 0




