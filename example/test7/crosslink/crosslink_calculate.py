import numpy as np
from multiprocessing import Pool, cpu_count

def read_data(data_file):
	data_line = data_file.readlines()
	n_atoms = int(data_line[2].split()[0])
	n_atom_types = int(data_line[3].split()[0])
	n_bonds = int(data_line[4].split()[0])
	n_bond_types = int(data_line[5].split()[0])
	xsize_1 = float(data_line[13].split()[1])
	xsize_2 = float(data_line[13].split()[0])
	xsize = xsize_1 - xsize_2
	ysize_1 = float(data_line[14].split()[1]) 
	ysize_2 = float(data_line[14].split()[0])
	ysize = ysize_1 - ysize_2
	zsize_1 = float(data_line[15].split()[1])
	zsize_2 = float(data_line[15].split()[0])
	zsize = zsize_1 - zsize_2
	mass = np.zeros(n_atom_types,dtype = float)
	radiu = np.zeros(n_atom_types,dtype = float)
	atom = np.zeros((n_atoms,7), dtype = float)
	atom_number = np.zeros(n_atoms, dtype = int)
	atom_type = np.zeros(n_atoms, dtype = int)
	Xs = np.zeros(n_atoms, dtype = float)
	Ys = np.zeros(n_atoms, dtype = float)
	Zs = np.zeros(n_atoms, dtype = float)
	bond = np.zeros((n_bonds,4), dtype = int)
	for i, line in enumerate(data_line):
		if line == 'Masses\n':
			for j in range(n_atom_types):
				mass[j] = float(data_line[i+j+2].split()[1])
		elif line == 'Pair Coeffs # lj/cut\n':
			for j in range(n_atom_types):
				radiu[j] = float(data_line[i+j+2].split()[2])/2
		elif line == 'Atoms # full\n':
			for j in range(n_atoms):
				a = data_line[i+j+2].split()
				atom_number[j] = int(a[0])
				atom_type[j] = int(a[2])
				Xs[j] = float(a[4]) - xsize_2
				Ys[j] = float(a[5]) - ysize_2
				Zs[j] = float(a[6]) - zsize_2
		elif line == 'Bonds\n':
			for j in range(n_bonds):
				a = data_line[i+j+2].split()
				bond[j,0] = int(a[0])#键序号
				bond[j,1] = int(a[1])#键类型
				bond[j,2] = int(a[2])#键原子1
				bond[j,3] = int(a[3])#键原子2
	return n_atoms, n_bonds, atom_number, atom_type, Xs, Ys, Zs, mass

def order_atoms(atom_number, atom_type, Xs, Ys, Zs):  #排序
	print("***************Order!***************\n")
	sorted_indices = np.argsort(atom_number)
	atom_number = np.take(atom_number, sorted_indices)
	atom_type = np.take(atom_type, sorted_indices)
	Xs = np.take(Xs, sorted_indices)
	Ys = np.take(Ys, sorted_indices)
	Zs = np.take(Zs, sorted_indices)
	print("*****Order Successfully!*****\n")
	return atom_number, atom_type, Xs, Ys, Zs

def read_bonds(bond_file, n_bonds):
	lines = bond_file.readlines()
	n_frames = 0
	line_frame = []
	row = 0
	for line in lines:
		if line == 'ITEM: ENTRIES index c_1[1] c_1[2] c_1[3]\n':
			n_frames += 1
			line_frame.append(row + 1)
		row += 1	
	bond_number = np.zeros((n_frames, n_bonds), dtype=int)
	bond_type = np.zeros((n_frames, n_bonds), dtype=int)
	bond_atom_1 = np.zeros((n_frames, n_bonds), dtype=int)
	bond_atom_2 = np.zeros((n_frames, n_bonds), dtype=int)
	for i in range(n_frames):
		row = line_frame[i]
		for j in range(n_bonds):
			b = lines[row + j].split()
			bond_number[i,j] = int(b[0])
			bond_type[i,j] = int(b[1])
			bond_atom_1[i,j] = int(b[2])
			bond_atom_2[i,j] = int(b[3])
	return n_frames, bond_number, bond_type, bond_atom_1, bond_atom_2

def calculate_single_frame(args):
	n_atoms, n_bonds, atom_type, bond_types, bond_atoms_1, bond_atoms_2, mass, frame = args
	bond_type = bond_types[frame,:]
	bond_atom_1 = bond_atoms_1[frame,:]
	bond_atom_2 = bond_atoms_2[frame,:]
	p_exopy, p_cn_0, p_cn_1, p_cn_2, p_oh, p_coc = extent_of_reaction(n_bonds, bond_atom_1, bond_atom_2, atom_type)
	Mw, num_molecules = weight_average_molecular_weight(n_atoms, n_bonds, bond_atom_1, bond_atom_2, atom_type, mass)
	return Mw, p_exopy, p_cn_0, p_cn_1, p_cn_2, p_oh, p_coc

def extent_of_reaction(n_bonds, bond_atoms_1, bond_atoms_2, atom_type):
	n_exopy = 0
	n_cn_0 = 0
	n_cn_1 = 0
	n_cn_2 = 0
	n_oh = 0
	n_coc = 0
	for i in range(n_bonds):
		if (atom_type[bond_atoms_1[i]-1] == 5 and atom_type[bond_atoms_2[i]-1] == 3):
			for j in range(n_bonds):
				if bond_atoms_1[j] == bond_atoms_1[i]:
					if (atom_type[bond_atoms_2[j]-1] == 2 or atom_type[bond_atoms_2[j]-1] == 1):
						n_cn_1 += 1
						break
					elif atom_type[bond_atoms_2[j]-1] == 3:
						n_cn_2 += 1
						break
				elif bond_atoms_2[j] == bond_atoms_1[i]:
					if (atom_type[bond_atoms_1[j]-1] == 2 or atom_type[bond_atoms_1[j]-1] == 1):
						n_cn_1 += 1
						break
					elif atom_type[bond_atoms_1[j]-1] == 3:
						n_cn_2 += 1
						break
		elif (atom_type[bond_atoms_1[i]-1] == 3 and atom_type[bond_atoms_2[i]-1] == 5):
			for j in range(n_bonds):
				if bond_atoms_1[j] == bond_atoms_2[i]:
					if (atom_type[bond_atoms_2[j]-1] == 2 or atom_type[bond_atoms_2[j]-1] == 1):
						n_cn_1 += 1
						break
					elif atom_type[bond_atoms_2[j]-1] == 3:
						n_cn_2 += 1
						break
				elif bond_atoms_2[j] == bond_atoms_2[i]:
					if (atom_type[bond_atoms_1[j]-1] == 2 or atom_type[bond_atoms_1[j]-1] == 1):
						n_cn_1 += 1
						break
					elif atom_type[bond_atoms_1[j]-1] == 3:
						n_cn_2 += 1
						break
		elif (atom_type[bond_atoms_1[i]-1] == 5 and (atom_type[bond_atoms_2[i]-1] == 1 or atom_type[bond_atoms_2[i]-1] == 2)):
			for j in range(n_bonds):
				if (bond_atoms_1[j] == bond_atoms_1[i] and j != i):
					if (atom_type[bond_atoms_2[j]-1] == 1 or atom_type[bond_atoms_2[j]-1] == 2):
						n_cn_0 += 1
						break
				elif bond_atoms_2[j] == bond_atoms_1[i]:
					if (atom_type[bond_atoms_1[j]-1] == 1 or atom_type[bond_atoms_1[j]-1] == 2):
						n_cn_0 += 1
						break
		elif (atom_type[bond_atoms_2[i]-1] == 5 and (atom_type[bond_atoms_1[i]-1] == 1 or atom_type[bond_atoms_1[i]-1] == 2)):
			for j in range(n_bonds):
				if bond_atoms_1[j] == bond_atoms_2[i]:
					if (atom_type[bond_atoms_2[j]-1] == 1 or atom_type[bond_atoms_2[j]-1] == 2):
						n_cn_0 += 1
						break
				elif (bond_atoms_2[j] == bond_atoms_2[i] and j != i):
					if (atom_type[bond_atoms_1[j]-1] == 1 or atom_type[bond_atoms_1[j]-1] == 2):
						n_cn_0 += 1
						break
		elif (atom_type[bond_atoms_1[i]-1] == 6 and atom_type[bond_atoms_2[i]-1] == 3):
			for j in range(n_bonds):
				if (bond_atoms_1[j] == bond_atoms_1[i]) and (j != i) and (atom_type[bond_atoms_2[j]-1] == 3):
					flag = True
					for k in range(n_bonds):
						if (bond_atoms_1[k] == bond_atoms_2[i] and bond_atoms_2[k] == bond_atoms_2[j]) or (bond_atoms_2[k] == bond_atoms_2[i] and bond_atoms_1[k] == bond_atoms_2[j]):
							flag = False
							n_exopy += 1
							break
					if flag:
						n_coc += 1
						break
				elif (bond_atoms_2[j] == bond_atoms_1[i]) and (atom_type[bond_atoms_1[j]-1] == 3):
					flag = True
					for k in range(n_bonds):
						if (bond_atoms_1[k] == bond_atoms_2[i] and bond_atoms_2[k] == bond_atoms_1[j]) or (bond_atoms_2[k] == bond_atoms_2[i] and bond_atoms_1[k] == bond_atoms_1[j]):
							flag = False
							n_exopy += 1
							break
					if flag:
						n_coc += 1
						break
		elif (atom_type[bond_atoms_1[i]-1] == 3 and atom_type[bond_atoms_2[i]-1] == 6):
			for j in range(n_bonds):
				if (bond_atoms_1[j] == bond_atoms_2[i]) and (atom_type[bond_atoms_2[j]-1] == 3):
					flag = True
					for k in range(n_bonds):
						if (bond_atoms_1[k] == bond_atoms_1[i] and bond_atoms_2[k] == bond_atoms_2[j]) or (bond_atoms_2[k] == bond_atoms_1[i] and bond_atoms_1[k] == bond_atoms_2[j]):							
							flag = False
							n_exopy += 1
							break
					if flag:
						n_coc += 1
						break
				elif (bond_atoms_2[j] == bond_atoms_2[i]) and (j != i) and (atom_type[bond_atoms_1[j]-1] == 3):
					flag = True
					for k in range(n_bonds):
						if (bond_atoms_1[k] == bond_atoms_1[i] and bond_atoms_2[k] == bond_atoms_1[j]) or (bond_atoms_2[k] == bond_atoms_1[i] and bond_atoms_1[k] == bond_atoms_1[j]):
							flag = False
							n_exopy += 1
							break
					if flag:
						n_coc += 1
						break
		elif (atom_type[bond_atoms_1[i]-1] == 6 and (atom_type[bond_atoms_2[i]-1] == 1 or atom_type[bond_atoms_2[i]-1] == 2)):
			n_oh += 1
		elif (atom_type[bond_atoms_2[i]-1] == 6 and (atom_type[bond_atoms_1[i]-1] == 1 or atom_type[bond_atoms_1[i]-1] == 2)):
			n_oh += 1
	p_exopy = n_exopy / 2000
	n_cn_2 -= 1000
	p_cn_0 = n_cn_0 / 2000
	p_cn_1 = n_cn_1 / 1000
	p_cn_2 = n_cn_2 / 1000
	p_oh = n_oh / 1000
	n_coc /= 2
	p_coc = n_coc / 1000
	return p_exopy, p_cn_0, p_cn_1, p_cn_2, p_oh, p_coc

def weight_average_molecular_weight(n_atoms, n_bonds, bond_atoms_1, bond_atoms_2, atom_type, mass):
	atom_mass = np.zeros(n_atoms)
	for atom in range(n_atoms):
		atom_mass[atom] = mass[atom_type[atom]-1]
	'''
	crosslink_point_dict = np.zeros(n_atoms)
	key = 0
	for k in range(100):
		for i in range(n_bonds):
			if (crosslink_point_dict[bond_atoms_1[i]-1] == 0) and (crosslink_point_dict[bond_atoms_2[i]-1] == 0):
				#如果两个原子key值=0(初次定义分子)，则定义两者的分子key值为key_number
				key += 1
				crosslink_point_dict[bond_atoms_1[i]-1] = key
				crosslink_point_dict[bond_atoms_2[i]-1] = key
			elif (crosslink_point_dict[bond_atoms_1[i]-1] == 0) and (crosslink_point_dict[bond_atoms_2[i]-1] != 0):
				#如果原子1的key值=0, 而原子2的key不是0, 则将原子1纳入原子2的分子中
				crosslink_point_dict[bond_atoms_1[i]-1] = crosslink_point_dict[bond_atoms_2[i]-1]
			elif (crosslink_point_dict[bond_atoms_1[i]-1] != 0) and (crosslink_point_dict[bond_atoms_2[i]-1] == 0):
				#如果原子2的key值=0, 而原子1的key不是0, 则将原子2纳入原子1的分子中
				crosslink_point_dict[bond_atoms_2[i]-1] = crosslink_point_dict[bond_atoms_1[i]-1]
			elif ((crosslink_point_dict[bond_atoms_1[i]-1] != 0) and 
				(crosslink_point_dict[bond_atoms_2[i]-1] != 0) and 
				(crosslink_point_dict[bond_atoms_1[i]-1] != crosslink_point_dict[bond_atoms_2[i]-1])):
				#如果原子1和2的key值都不是0, 且二者不等, 则将原子2的分子合并到原子1的分子
				for j in range(1,key+1):
					if crosslink_point_dict[j]==crosslink_point_dict[bond_atoms_2[i]-1]:
						crosslink_point_dict[j]=crosslink_point_dict[bond_atoms_1[i]-1]
	mol = np.zeros(key)
	mol_mass = np.zeros(key)
	for i in range(key):
		for j in range(n_atoms):
			if crosslink_point_dict[j] == i+1:
				mol[i] += 1
				mol_mass[i] += atom_mass[j]
	valid_indices = [i for i in range(key) if mol[i] > 0]
	num_molecules = len(valid_indices)
	total_mass = 0.0
	total_square = 0.0
	for i in valid_indices:
		total_mass += mol_mass[i]
		total_square += mol_mass[i] * mol_mass[i]
	if num_molecules == 0:
		Mn = Mw = 0.0
	else:
		Mn = total_mass / num_molecules
		Mw = total_square / total_mass


	'''
	# 2. 使用并查集（Union-Find）识别分子
	parent = np.arange(n_atoms)  # 每个原子的父节点，初始为自己
	def find(x):
		# 路径压缩
		while parent[x] != x:
			parent[x] = parent[parent[x]]
			x = parent[x]
		return x

	def union(x, y):
		rx, ry = find(x), find(y)
		if rx != ry:
			parent[ry] = rx  # 将ry的根指向rx

	# 遍历所有键，合并原子所在的集合
	for i in range(n_bonds):
		a = bond_atoms_1[i] - 1  # 转换为0-indexed
		b = bond_atoms_2[i] - 1
		union(a, b)
	# 3. 统计每个分子的总质量和分子数目
	# 先将每个原子的根找到，并压缩路径
	for i in range(n_atoms):
		find(i)
	molecule_mass = {}
	for i in range(n_atoms):
		root = parent[i]
		molecule_mass[root] = molecule_mass.get(root, 0.0) + atom_mass[i]
	num_molecules = len(molecule_mass)
	# 4. 计算重均分子量 Mw = Σ (mi^2) / Σ mi
	total_mass = 0.0
	total_square = 0.0
	for mass_val in molecule_mass.values():
		total_mass += mass_val
		total_square += mass_val * mass_val
	if total_mass == 0:
		return 0.0, 0
	Mw = total_square / total_mass
	return Mw, num_molecules

def init_worker():
	import numpy as np
	np.random.seed()

if __name__ == '__main__':
	with open('crosslink.data') as data_file:
		n_atoms, n_bonds, atom_number, atom_type, Xs, Ys, Zs, mass = read_data(data_file)
	atom_number, atom_type, Xs, Ys, Zs = order_atoms(atom_number, atom_type, Xs, Ys, Zs)
	with open('bond.dump') as bond_file:
		n_frames, bond_number, bond_type, bond_atom_1, bond_atom_2 = read_bonds(bond_file, n_bonds)
	with Pool(processes=cpu_count(), initializer=init_worker, maxtasksperchild=10) as pool:
		args_list = [(n_atoms, n_bonds, atom_type, bond_type, bond_atom_1, bond_atom_2, mass, frame) for frame in range(n_frames)]
		results = pool.map(calculate_single_frame, args_list, chunksize=5)
	Mw = np.array([result[0] for result in results])
	p_exopy = np.array([result[1] for result in results])
	p_cn_0 = np.array([result[2] for result in results])
	p_cn_1 = np.array([result[3] for result in results])
	p_cn_2 = np.array([result[4] for result in results])
	p_oh = np.array([result[5] for result in results])
	p_coc = np.array([result[6] for result in results])
	time = 0.1*np.arange(n_frames)
	with open('Mw.txt', 'w') as file:
		for i in range(n_frames):
			file.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'%(time[i], p_exopy[i], p_cn_0[i], p_cn_1[i], p_cn_2[i], p_oh[i], p_coc[i], Mw[i]))