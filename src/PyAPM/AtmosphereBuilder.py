import subprocess
import shutil
import numpy as np
import time
import os

class AtmosphereBuilder:
	def __init__(self, box, n_chain, n_component):
		#box: box size length
		#n_chain: number of each component, vector length should equal n_component
		#n_component: number of components
		self.box = box
		self.n_chain = n_chain
		self.n_component = n_component

	def run(self, command, check_file, timeout=3600):
		#command (list or str): command to be executed, e.g., ['packmol', '-i', 'mixture.inp']
		#check_file (str): name of the file that should be generated upon successful execution
		#timeout (int): maximum waiting time (seconds)
		if isinstance(command, list) and command[0] == 'packmol':
			packmol_path = self._get_packmol_path()
			shell_cmd = f'{packmol_path} < mixture.inp'
			print(f"Running: {shell_cmd}")
			process = subprocess.Popen(shell_cmd, shell=True,
				stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
				text=True)
		else:
			print(f"Running command: {' '.join(command) if isinstance(command, list) else command}")
			process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		try:
			stdout, stderr = process.communicate(timeout=timeout)
		except subprocess.TimeoutExpired:
			process.kill()
			print(f"Command timed out after {timeout} seconds.")
			return False
		if process.returncode != 0:
			print(f"Command failed with error: {stderr}")
			return False
		# Check if the output file exists
		if os.path.exists(check_file):
			print(f"Success: {check_file} generated.")
			return True
		else:
			print(f"Error: {check_file} not found after command execution.")
			return False

	def write_system(self):
		#Write the system file in Moltemplate, read the single-chain lt and multi-chain pdb, 
		#and generate LAMMPS input files (system.data, system.in.init, system.in.settings)
		moltemplate_system = open('system.lt','w')
		for i in range(self.n_component):
			moltemplate_system.write('import "Polymer_%d.lt"\n\n'%(i+1))
		for i in range(self.n_component):
			moltemplate_system.write('molecule%d = new Polymer%d [%d].move(0,0,10)\n'%(i+1, i+1, self.n_chain[i]))
		moltemplate_system.write('\n')
		moltemplate_system.write('write_once("Data Boundary")\n {\n')
		moltemplate_system.write('   0.0    %d     xlo xhi\n   0.0    %d     ylo yhi\n   0.0    %d     zlo zhi\n}\n'%(self.box[0],self.box[1],self.box[2]))
		moltemplate_system.close()
		#Write the packmol input file to convert single-chain pdb to multi-chain pdb
		inp = open('mixture.inp','w')
		inp.write('\ntolerance 2.5\nfiletype pdb\noutput Polymer.pdb\n\n')
		for i in range(self.n_component):
			inp.write('structure Chain_%d.pdb\n  number %d\n  inside box 0. 0. 0. %d. %d. %d.\nend structure\n'%(i+1, self.n_chain[i], self.box[0],self.box[1],self.box[2]))
		inp.close()
		return 0

	def change_lj_coeff(self, file1 = 'system.in.init', file2 = 'system.data', file3 = 'system.in.settings'):
		#Replace the LJ interaction parameters in the DREIDING force field, 
		#and fix some errors in the input files written by Moltemplate
		with open(file1) as init_file:
			lines = init_file.readlines()
		with open('systemnew.in.init','w') as new_init_file:
			for line in lines:
				try:
					if line.split()[0] == 'pair_style':
						line = 'pair_style      lj/cut/coul/long 10.0\n'
				except IndexError:
					continue
				new_init_file.write('%s'%line)
		shutil.move('systemnew.in.init', file1)
		with open(file2) as date_file:
			lines = date_file.readlines()
			for i in range(len(lines)):
				if lines[i] == 'Masses\n':
					row_1 = i
				elif lines[i] == 'Atoms  # full\n' or lines[i] == 'Atoms\n':
					row_2 = i
					break
			atom = []
			atomtype = []
			for i in range(row_1+2, row_2-1):
				print(lines[i])
				atom.append(int(lines[i].split()[0]))
				atomtype.append(str(lines[i].split()[-1]))
			atom = np.array(atom)
			atomtype = np.array(atomtype)
			sorted_indices = np.argsort(atom)
			atomtype = np.take(atomtype, sorted_indices)
		with open(file3) as settings_file:
			lines = settings_file.readlines()
		with open('systemnew.in.settings','w') as newsettings_file:
			for line in lines:
				coeff = line.split()
				if coeff[0] == 'pair_coeff':
					coeff[3] = ' '
					if coeff[1] == coeff[2]:
						a = atomtype[int(coeff[1])-1]
						if a == 'H':
							coeff[4], coeff[5] = '0.01', '2.576'
						elif a == 'H_HB':
							coeff[4], coeff[5] = '0.00001', '2.576'
						elif a == 'C_3' or a == 'C_2' or a == 'C_R':
							coeff[4], coeff[5] = '0.097', '3.434'
						elif a == 'O_3' or a == 'O_2':
							coeff[4], coeff[5] = '0.218', '3.118'
						elif a == 'N_3' or a == 'N_2' or a == 'N_R' or a == 'N_2_d1_hd' or a == 'N_2_b2_d1' or a == 'N_R_d1':
							coeff[4], coeff[5] = '0.148', '3.277'
						elif a == 'F':
							coeff[4], coeff[5] = '0.308', '2.92'
						elif a == 'Cl':
							coeff[4], coeff[5] = '0.308', '3.48'
						elif a == 'Br':
							coeff[4], coeff[5] = '0.308', '3.747'
						elif a == 'S_3':
							coeff[4], coeff[5] = '0.218', '3.677'
						newsettings_file.write('%s %s %s %s %s\n'%(coeff[0], coeff[1], coeff[2], coeff[4], coeff[5]))
				else:
					newsettings_file.write('%s'%line)
		shutil.move('systemnew.in.settings', file3)
		return 0

	def _get_packmol_path(self):
		#Get the PACKMOL environment variable
		packmol_path = os.environ.get('PACKMOL_EXEC')
		if packmol_path is None:
			packmol_path = shutil.which('packmol')
		if packmol_path is None:
			raise RuntimeError(
				"packmol not found. Please set the environment variable PACKMOL_EXEC "
				"to the full path of packmol, or ensure 'packmol' is in your PATH."
			)
		return packmol_path