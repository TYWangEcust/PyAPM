import re

#CrosslinkBuilder function library
def get_data(sectionName, lines, sectionIndexList, useExcept = True):
	if useExcept: # Checks that section name is existing in LAMMPS data
		try:
			startIndex = lines.index(sectionName)
		except ValueError:
			# If doesn't exist, return empty list that can be added as normal to main list later
			data = []
			return data

	else: # Allows for later try/except blocks to catch missing section names
		startIndex = lines.index(sectionName)
	endIndex = sectionIndexList[sectionIndexList.index(startIndex) + 1]
	data = lines[startIndex+1:endIndex] # +1 means sectionName doesn't get included
	data = [val.split() for val in data]
	return data

def find_sections(lines):
	# Find index of section keywords - isalpha works as no spaces, newlines or punc in section keywords
	sectionIndexList = [lines.index(line) for line in lines if line.isalpha()]
	# Add end of file as last index
	sectionIndexList.append(len(lines))
	return sectionIndexList

def clean_data(lines):
	# Remove blank lines
	lines = [line for line in lines if line != '\n']

	# Remove comments - negative lookbehind means label comments in masses are kept e.g # C_3
	lines = [re.sub(r'(?<!\d\s\s)#(.*)', '' ,line) for line in lines]

	# Remove newline terminators
	lines = [re.sub(r'\n', '', line) for line in lines]

	# Remove empty strings in list caused by comments being removed
	lines = [line for line in lines if line != '']

	# Remove comments - negative lookbehind means label comments in masses are kept e.g # C_3
	lines = [re.sub(r'\s+$', '', line) for line in lines]

	return lines

def get_style(masses):
	potencial_element =[]
	for i in range(0,len(masses)):
		mass = masses[i][1] 
		if mass == '1.008':
			mass_str = 'H'
			potencial_element += mass_str + ' '
		elif mass == '12.011':    
			mass_str = 'C'
			potencial_element += mass_str + ' '
		elif mass == '14.007':
			mass_str = 'N'
			potencial_element += mass_str + ' '
		elif mass == '15.999':
			mass_str = 'O'
			potencial_element += mass_str + ' '
		elif mass == '32.06':
			mass_str = 'S'
			potencial_element += mass_str + ' '
	
	return potencial_element

def remove_first_n_lines(file_path, n):
	n = int(n)
	with open(file_path, 'r', encoding='utf-8') as file:
		lines = file.readlines()
	if n > len(lines):
		print(f"Warning: Attempted to remove {n} lines from a file with only {len(lines)} lines.")
		n = len(lines)
	lines = lines[n:]
	with open(file_path, 'w', encoding='utf-8') as file:
		file.writelines(lines)
	return 0 


def change_lj_coeff(masses, element_str):
	num = 1
	add_word = []
	for i in range(len(masses)-1):
		add_word.append([])
	for element in element_str:
		print(element)
		if element == 'H':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.01 2.576\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'H_HB':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.00001 2.576\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'C_3':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.097 3.434\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'C_2':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.097 3.434\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'C_R':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.097 3.434\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'N_3':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.148 3.277\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'N_2' or element == 'N_2_d1_hd' or element == 'N_2_b2_d1':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.148 3.277\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'N_R' or element == 'N_R_d1':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.148 3.277\n'.format(n1=num, n2=num))
			num = num + 1	
		elif element == 'O_3':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.218 3.118\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'O_2':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.218 3.118\n'.format(n1=num, n2=num))
			num = num + 1
		elif element == 'S_3':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.218 3.677\n'.format(n1=num, n2=num))
			num = num + 1		
		elif element == 'P_3':
			with open('cleanedsystemself.in.settings','w') as add_word:
				add_word.write('pair_coeff {n1} {n2} lj/cut/coul/long 0.218 3.815\n'.format(n1=num, n2=num))
			num = num + 1
	
	source_file_path = 'cleanedsystem.in.settings'  
	target_file_path = 'cleanedsystemself.in.settings'
	with open(source_file_path, 'r', encoding='utf-8') as source_file:  
		content = source_file.read()  
	
	with open(target_file_path, 'w', encoding='utf-8') as target_file:  
		target_file.write(content)  
	return 0

def write_lammps_crosslink(n_react, react_temp):
	with open('in.md','w') as in_file:
		in_file.write('include\t"system.in.init"\n')
		in_file.write('read_data\t"cleanedsystem.data" extra/bond/per/atom 2 extra/angle/per/atom 12 extra/dihedral/per/atom 12 extra/improper/per/atom 12 extra/special/per/atom 12\n')
		in_file.write('include\t"cleanedsystem.in.settings"\n\n')
		for n in range(n_react):
			in_file.write('molecule\tpre%d pre-molecule%d.data\nmolecule\tpost%d post-molecule%d.data\n\n'%(n,n,n,n))
		in_file.write('neighbor\t2.5 bin\nneigh_modify\tevery 1 delay 0 check yes\n\n')
		in_file.write('timestep\t1.0\n\nmin_style\tcg\nminimize\t1e-08 1e-10 1000000 10000000\nreset_timestep\t0\n\n')
		in_file.write('thermo_style\tcustom step temp epair emol etotal press vol density\nthermo\t10000\n\n')
		in_file.write('dump\t1 all custom 10000 equil.lammpstrj id mol type q x y z ix iy iz\ndump_modify\t1 sort id\n\n')
		in_file.write('fix\t1 all nvt temp 300 300 100\nrun\t100000\nunfix\t1\n\n')
		in_file.write('fix\t2 all npt temp 300 300 100 iso 3000 3000 1000\nrun\t300000\nunfix\t2\n\n')
		in_file.write('variable\tn loop 3\nlabel\there\n\n')
		in_file.write('fix\t3 all nvt temp 800 800 100\nrun\t100000\nunfix\t3\n\n')
		in_file.write('fix\t4 all nvt temp 300 300 100\nrun\t100000\nunfix\t4\n\n')
		in_file.write('fix\t5 all npt temp 300 300 100 iso 1000 1000 1000\nrun\t300000\nunfix\t5\n\n')
		in_file.write('next\tn\njump\tSELF here\n\n')
		in_file.write('fix\t6 all npt temp 300 300 100 iso 1 1 1000\nrun\t1000000\nunfix\t6\n')
		in_file.write('fix\t7 all npt temp 300 %.1f 100 iso 1 1 1000\nrun\t1000000\nunfix\t7\nundump\t1\n\n'%react_temp)
		in_file.write('write_data\tequil.data pair ij\nwrite_restart\tequil.restart\n\n')
		in_file.write('fix\tfxrct all bond/react stabilization yes statted_grp .03 ')
		for n in range(n_react):
			in_file.write('react rxn%d all 1000 1.0 4.0 pre%d post%d automap%d.data stabilize_steps 100 prob %.3f 114514 '%(n,n,n,n,1/n_react))
		in_file.write('\n')
		in_file.write('fix\t8 statted_grp_REACT npt temp %.1f %.1f 100 iso 1 1 1000\n'%(react_temp, react_temp))
		in_file.write('fix\t9 bond_react_MASTER_group temp/rescale 1 %.1f %.1f 1 1\n'%(react_temp, react_temp))
		in_file.write('thermo_style\tcustom step temp epair emol etotal press vol density f_fxrct[1]\n')
		in_file.write('dump\t2 all custom 10000 crosslink.lammpstrj id mol type q x y z ix iy iz\ndump_modify\t2 sort id\n\n')
		in_file.write('run\t500000\nwrite_data\tcrosslink.data\nunfix\t8\nunfix\t9\n')
	return 0