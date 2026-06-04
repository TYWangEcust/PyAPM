import numpy as np
import pandas as pd
import math
from . import fit

def readlammpstrj(lammpstrj):
	lammpstrj_lines = lammpstrj.readlines()
	steps = 0
	line_2 = []
	row = 0
	for line_1 in lammpstrj_lines:
		row = row + 1
		if line_1 == 'ITEM: NUMBER OF ATOMS\n':
			atoms = int(lammpstrj_lines[row].rstrip())
			steps = steps + 1
			line_2.append(row - 2)
	n = 0
	atom_number = np.zeros((atoms, steps))
	atom_type = np.zeros((atoms, steps))
	Xs = np.zeros((atoms, steps))
	Ys = np.zeros((atoms, steps))
	Zs = np.zeros((atoms, steps))
	V = []
	NumOfBins = 100 #定义分割数目
	row = 1
	while n < steps: #对于每一步
		row = line_2[n]-1
		for line_1 in lammpstrj_lines:
			row = row + 1
			if n < steps - 1:
				if ((row <= line_2[n+1]) and (row >= line_2[n])):
					if line_1 == 'ITEM: BOX BOUNDS pp pp pp\n':
						bx = lammpstrj_lines[row].split()
						by = lammpstrj_lines[row+1].split()
						bz = lammpstrj_lines[row+2].split()
						xsize = (float(bx[1])-float(bx[0]))
						ysize = (float(by[1])-float(by[0]))
						zsize = (float(bz[1])-float(bz[0]))
						if xsize <= ysize:
							rsize=xsize
						else: 
							rsize=ysize
						if zsize<=rsize: 
							rsize=zsize             

					if line_1 == 'ITEM: ATOMS id type xs ys zs\n':
						i = 0
						while i <= (atoms - 1):
							a = lammpstrj_lines[row + i].split()
							atom_number[i][n] = int(a[0])
							atom_type[i][n] = int(a[1])
							Xs[i][n] = float(a[2])
							Ys[i][n] = float(a[3])          
							Zs[i][n] = float(a[4])          
							i = i + 1
			else:
				if ((row >= line_2[n])and (row<=len(lammpstrj_lines))):
					if line_1 == 'ITEM: BOX BOUNDS pp pp pp\n':
						bx = lammpstrj_lines[row].split()
						by = lammpstrj_lines[row+1].split()
						bz = lammpstrj_lines[row+2].split()
						xsize = (float(bx[1])-float(bx[0]))
						ysize = (float(by[1])-float(by[0]))
						zsize = (float(bz[1])-float(bz[0]))
						if xsize <= ysize:
							rsize=xsize
						else: 
							rsize=ysize
						if zsize<=rsize: 
							rsize=zsize   
					if line_1 == 'ITEM: ATOMS id type xs ys zs\n':
						i = 0
						while i <= (atoms - 1): 
							a = lammpstrj_lines[row + i].split()
							atom_number[i][n] = int(a[0])
							atom_type[i][n] = int(a[1])
							Xs[i][n] = float(a[2])
							Ys[i][n] = float(a[3])          
							Zs[i][n] = float(a[4])          
							i = i + 1
		n = n + 1
	for n in range(0,steps,1):
		for i in range(0,atoms,1):
			if Xs[i][n] >= 1:
				Xs[i][n] -= 1
			elif Xs[i][n] <= 0:
				Xs[i][n] += 1
			Xs[i][n] *= xsize
			if Ys[i][n] >= 1:
				Ys[i][n] -= 1
			elif Ys[i][n] <= 0:
				Ys[i][n] += 1
			Ys[i][n] *= ysize
			if Zs[i][n] >= 1:
				Zs[i][n] -= 1
			elif Zs[i][n] <= 0:
				Zs[i][n] += 1
			Zs[i][n] *= zsize
	Volumn = xsize * ysize * zsize
	rsize = 32
	dr = 0.5 * rsize / float(NumOfBins)
	print("\n*************Calculating Radius Distribution Function************\n")
	print("Number of beads in calculation =  " + str(atoms) +", Number of Bead Type in calculation = " + str(np.max(atom_type)) +"\n")
	print("Point Number for discreting r =" +str(NumOfBins) +", Number of step for calculation = " +str(steps)+"\n")
	print("\t End of read data ........\n")
	return atom_number, atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr

def order_atoms(atom_number, atom_type, Xs, Ys, Zs, atoms, j):
	print("***************Order!***************\n")
	sorted_indices = np.argsort(atom_number[:,j])
	atom_number[:,j] = np.take(atom_number[:,j], sorted_indices)
	atom_type[:,j] = np.take(atom_type[:,j], sorted_indices)
	Xs[:,j] = np.take(Xs[:,j], sorted_indices)
	Ys[:,j] = np.take(Ys[:,j], sorted_indices)
	Zs[:,j] = np.take(Zs[:,j], sorted_indices)
	print("*****Order Successfully!*****\n")
	return atom_number, atom_type, Xs, Ys, Zs

def read_data(temp):
	data = open('25-%d-normalizing.data'%temp)
	lines = data.readlines()
	a = lines[2].split()
	Atoms = int(a[0])#总键数目
	molnumber = np.zeros(Atoms, dtype=int)
	row = 0
	for line in lines:
		if line == "Atoms # full\n":
			line_Atoms = row
			break
		row += 1
	for i in range(Atoms):
		a = lines[line_Atoms+2+i].split()
		molnumber[i] = int(a[1]) 
	mol = np.max(molnumber)
	data.close()
	return mol

def find_target(atoms, atom_type, temp):
	mol = read_data(temp)
	target = np.where(atom_type[:,-1] == 2)[0]
	len_B = target[2] - target[1]
	n_A = np.zeros((len(target),6), dtype = int)
	n_B = np.zeros((len(target),len_B-3), dtype = int)
	n_mol = np.zeros(len(target), dtype = int)
	target_per_mol = int(len(target)/mol)
	print(target_per_mol)
	for i in range(mol):
		for j in range(1, target_per_mol):
			for k in range(3):
				n_A[i * target_per_mol + j-1, k] = target[i * target_per_mol + 1] - 3*j + k
			for k in range(3):
				n_A[i * target_per_mol + j-1, k + 3] = target[i * target_per_mol + j] + k
			for k in range(len_B-3):
				n_B[i * target_per_mol + j-1, k] = target[i * target_per_mol + j] + 3 + k
			n_mol[i * target_per_mol + j-1] = i + 1
		for k in range(6):
			n_A[i * target_per_mol + target_per_mol - 1, k] = target[i * target_per_mol] - 1 + k
		n_B[i * target_per_mol + target_per_mol - 1, 0] = target[i * target_per_mol] - 1
		n_mol[i * target_per_mol + target_per_mol - 1] = i + 1
	return n_mol, n_A, n_B

def RDF(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, type_cg, bond, angle, target_1, target_2):  #计算径向分布函数gr
	
	NumOfBins = int(NumOfBins)
	rdf = np.zeros((NumOfBins+1, steps))
	r_space = dr*np.linspace(1e-14,NumOfBins,NumOfBins+1)
	rdf_ave = np.zeros(NumOfBins+1)
	print("START..........<*_*>\n")
	print("Atom 1: %d\tAtom 2: %d\n"%(target_1,target_2))
	type_1 = type_cg == target_1
	type_2 = type_cg == target_2
	if (target_1 == 1 and target_2 == 1):
		bond_atom_1, bond_atom_2, angle_atom_1, angle_atom_2 = special_atom(bond,angle)
	for j in range(steps): #对于每步
		N1 = np.sum(type_1)
		N2 = np.sum(type_2)
		print("N1=%d\tN2=%d\n"%(N1,N2))
		Xs_1 = np.array([X_cg[:,j][type_1[:]]])
		Xs_2 = np.array([X_cg[:,j][type_2[:]]])
		Ys_1 = np.array([Y_cg[:,j][type_1[:]]])
		Ys_2 = np.array([Y_cg[:,j][type_2[:]]])
		Zs_1 = np.array([Z_cg[:,j][type_1[:]]])
		Zs_2 = np.array([Z_cg[:,j][type_2[:]]])
		dx = np.abs(np.subtract.outer(Xs_1,Xs_2).flatten())
		dx = np.minimum(dx, np.abs(xsize - dx))
		dy = np.abs(np.subtract.outer(Ys_1,Ys_2).flatten())
		dy = np.minimum(dy, np.abs(ysize - dy))
		dz = np.abs(np.subtract.outer(Zs_1,Zs_2).flatten())
		dz = np.minimum(dz, np.abs(zsize - dz))
		coords = np.array([dx, dy, dz])
		r_temp = np.sqrt(np.sum(coords ** 2, axis=0))
		r_temp_bin = np.digitize(r_temp / dr, np.linspace(1e-14,NumOfBins,NumOfBins+1))
		r_temp_bin[np.where(r_temp_bin > NumOfBins)] = 0
		np.add.at(rdf, (r_temp_bin, j), 1)
		if (target_1 == 1 and target_2 == 1): 
			rdf_12=RDF_special(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, bond_atom_1, bond_atom_2, j)
			rdf_13=RDF_special(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, angle_atom_1, angle_atom_2, j)
			rdf[:,j] -= 2*rdf_12 + 2*rdf_13
	rdf_ave = np.mean(rdf, axis=1) 
	r_cut = np.where(rdf_ave < 1)[0][-1]
	rdf_ave[0:r_cut+1] = 0
	rdf_ave[r_cut+1:-1] = fit.savitzky_golay(rdf_ave[r_cut+1:-1],9, 2, deriv=0, rate=1)	
	for r in range(NumOfBins): 
		if r>=1:
			r_u = float(r_space[r])
			r_b = float(r_space[r-1])
			rdf_ave[r] /= ((N1 * N2) * (4 / 3) * math.pi * (r_u ** 3 - r_b ** 3) / Volumn)
		rdf_ave[0] = 0 
	print("Finished..........<*_*>\n")
	return rdf_ave

def RDF_special(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr,special_atom_1, special_atom_2, j):
	r_space = dr*np.linspace(1e-14,NumOfBins,NumOfBins+1)
	rdf = np.zeros(NumOfBins+1)
	rdf_ave = []
	print("START Calculate Special_bonds..........<*_*>\n")
	coords_1 = pd.DataFrame({'x': X_cg[special_atom_1,j], 'y': Y_cg[special_atom_1,j], 'z': Z_cg[special_atom_1,j]})
	coords_2 = pd.DataFrame({'x': X_cg[special_atom_2,j], 'y': Y_cg[special_atom_2,j], 'z': Z_cg[special_atom_2,j]})
	coords_1 = coords_1.reindex(coords_2.index)
	dx = np.abs(coords_1.x - coords_2.x.values)
	dy = np.abs(coords_1.y - coords_2.y.values)
	dz = np.abs(coords_1.z - coords_2.z.values)
	dx = np.minimum(dx, np.abs(xsize - dx))
	dy = np.minimum(dy, np.abs(ysize - dy))
	dz = np.minimum(dz, np.abs(zsize - dz))
	coords = np.array([dx, dy, dz])
	r_temp = np.sqrt(np.sum(coords ** 2, axis=0))
	r_temp_bin = np.digitize(r_temp / dr, np.linspace(1e-14,NumOfBins,NumOfBins+1))
	r_temp_bin[np.where(r_temp_bin > NumOfBins)] = 0
	np.add.at(rdf, (r_temp_bin), 1) 
	return rdf

def special_atom(bond,angle):
	bond_atom_1 = [int(bond[i,2] - 1) for i in range(len(bond[:,0])) if bond[i,1] == 1]
	bond_atom_2 = [int(bond[i,3] - 1) for i in range(len(bond[:,0])) if bond[i,1] == 1]
	angle_atom_1 = [int(angle[i,2] - 1) for i in range(len(angle[:,0])) if angle[i,1] == 1]
	angle_atom_2 = [int(angle[i,4] - 1) for i in range(len(angle[:,0])) if angle[i,1] == 1]
	return bond_atom_1, bond_atom_2, angle_atom_1, angle_atom_2

def Bond_distribution(X_cg, Y_cg, Z_cg, xsize, ysize, zsize, steps, bond, bondtype):
	print("Start Calculating Bond %d\n"%bondtype)
	NumOfBins = 100
	bondlength = []
	N=0
	for i in range(len(bond[:,0])):
		if (bond[i,1]==bondtype):
			for j in range(steps):
				N+=1
				dx = np.abs(X_cg[int(bond[i,2])-1,j]-X_cg[int(bond[i,3])-1,j])
				x_dist = np.minimum(dx, np.abs(xsize - dx))
				dy = np.abs(Y_cg[int(bond[i,2])-1,j]-Y_cg[int(bond[i,3])-1,j])
				y_dist = np.minimum(dy, np.abs(ysize - dy))
				dz = np.abs(Z_cg[int(bond[i,2])-1,j]-Z_cg[int(bond[i,3])-1,j])
				z_dist = np.minimum(dz, np.abs(zsize - dz))
				bondlength.append(np.sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist))
	ave_bondlength = np.sum(bondlength)/N
	l_u = 4.5
	l_b = 2
	bdf = np.zeros((NumOfBins,2))
	l_temp_bin = np.digitize(bondlength,np.linspace(l_b, l_u, NumOfBins))
	np.add.at(bdf, (l_temp_bin - 1, 1), 1)
	bdf[:,0] = np.linspace(l_b, l_u, NumOfBins)
	bdf[-1,1] = 0
	bdf[-2,1] = 0
	bdf[:,1] /= N
	bdf[:,1] /= bdf[:,0]**2
	bdf[:,1] = fit.savitzky_golay(bdf[:,1], 9, 2, deriv=0, rate=1)
	bdf[:,1] = np.clip(bdf[:,1], 0, None)
	m = np.max(bdf[:,1])
	bdf[:,1] /= m
	print("Finished..........<*_*>\n")
	return ave_bondlength, bdf

def Angel_distribution(X_cg, Y_cg, Z_cg, xsize, ysize, zsize, steps, angle, angletype):
	print("Start Calculating Angle %d\n"%angletype)
	N = 0
	angle_degree = []
	for i in range(len(angle[:,0])):
		if (angle[i,1]==angletype):
			for j in range(steps):
				N += 1
				dx_1 = X_cg[int(angle[i,3])-1,j] - X_cg[int(angle[i,2])-1,j]
				dx_1 = np.mod(dx_1 + xsize/2, xsize) - xsize/2
				dy_1 = Y_cg[int(angle[i,3])-1,j] - Y_cg[int(angle[i,2])-1,j]
				dy_1 = np.mod(dy_1 + ysize/2, ysize) - ysize/2
				dz_1 = Z_cg[int(angle[i,3])-1,j] - Z_cg[int(angle[i,2])-1,j]
				dz_1 = np.mod(dz_1 + zsize/2, zsize) - zsize/2
				dx_2 = X_cg[int(angle[i,4])-1,j] - X_cg[int(angle[i,3])-1,j]
				dx_2 = np.mod(dx_2 + xsize/2, xsize) - xsize/2
				dy_2 = Y_cg[int(angle[i,4])-1,j] - Y_cg[int(angle[i,3])-1,j]
				dy_2 = np.mod(dy_2 + ysize/2, ysize) - ysize/2
				dz_2 = Z_cg[int(angle[i,4])-1,j] - Z_cg[int(angle[i,3])-1,j]
				dz_2 = np.mod(dz_2 + zsize/2, zsize) - zsize/2
				v1 = np.array([dx_1, dy_1, dz_1])
				v2 = np.array([dx_2, dy_2, dz_2])
				dtheta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
				angle_degree.append(np.degrees(np.pi - np.arccos(dtheta)))
	ave_angle_degree = np.sum(angle_degree)/N
	theta_l = 60
	theta_u = 180
	adf = np.zeros((121,2))
	degree_temp_bin = np.digitize(angle_degree,np.linspace(theta_l,theta_u,121))
	np.add.at(adf, (degree_temp_bin - 1, 1), 1)
	adf[:,0] = np.linspace(theta_l,theta_u,121)
	adf[:,1] /= N
	adf[1:,1] /= np.sin(np.deg2rad(adf[1:,0]))
	adf[-3:,1] = 0
	adf[:,1] = fit.savitzky_golay(adf[:,1], 9, 2, deriv=0, rate=1)
	adf[:,1] = np.clip(adf[:,1], 0, None)
	m = np.max(adf[:,1])
	adf[:,1] /= m
	print("Finished..........<*_*>\n")
	return ave_angle_degree, adf

def Dihedral_distributation_BAAB(atoms, Xs, Ys, Zs, xsize, ysize, zsize, steps, n_mol, n_A, n_B):
	print("Start Calculating Dihedras BAAB\n")
	dihedral_degree = []
	N = 0
	for i in range(len(n_A)-1):
		if (n_mol[i]==n_mol[i+1]):
			for j in range(steps):
				N += 1
				dx_1 = Xs[n_A[i]-1][j] - Xs[n_B[i]-1][j]
				dx_1 = np.mod(dx_1 + xsize, xsize * 2) - xsize
				dy_1 = Ys[n_A[i]-1][j] - Ys[n_B[i]-1][j]
				dy_1 = np.mod(dy_1 + ysize, ysize * 2) - ysize
				dz_1 = Zs[n_A[i]-1][j] - Zs[n_B[i]-1][j]
				dz_1 = np.mod(dz_1 + zsize, zsize * 2) - zsize
				dx_2 = Xs[n_A[i+1]-1][j] - Xs[n_A[i]-1][j]
				dx_2 = np.mod(dx_2 + xsize, xsize * 2) - xsize
				dy_2 = Ys[n_A[i+1]-1][j] - Ys[n_A[i]-1][j]
				dy_2 = np.mod(dy_2 + ysize, ysize * 2) - ysize
				dz_2 = Zs[n_A[i+1]-1][j] - Zs[n_A[i]-1][j]
				dz_2 = np.mod(dz_2 + zsize, zsize * 2) - zsize
				dx_3 = Xs[n_B[i+1]-1][j] - Xs[n_A[i+1]-1][j]
				dx_3 = np.mod(dx_3 + xsize, xsize * 2) - xsize
				dy_3 = Ys[n_B[i+1]-1][j] - Ys[n_A[i+1]-1][j]
				dy_3 = np.mod(dy_3 + ysize, ysize * 2) - ysize
				dz_3 = Zs[n_B[i+1]-1][j] - Zs[n_A[i+1]-1][j]
				dz_3 = np.mod(dz_3 + zsize, zsize * 2) - zsize
				v1 = np.array([dx_1, dy_1, dz_1])
				v2 = np.array([dx_2, dy_2, dz_2])
				v3 = np.array([dx_3, dy_3, dz_3])
				dihedral_degree.append(np.degrees(np.arccos(np.dot(np.cross(v1,v2),np.cross(v2,v3))/(np.linalg.norm(v1)*np.linalg.norm(v2)*np.linalg.norm(v2)*np.linalg.norm(v3)))))
	ave_dihedral_degree = sum(dihedral_degree)/N
	ddf = np.zeros((181,2))
	degree_temp_bin = np.digitize(dihedral_degree,np.linspace(0,180,181))
	np.add.at(ddf, (degree_temp_bin - 1, 1), 1)
	m = 0
	ddf[:,0] = np.linspace(0,180,181)
	ddf[:,1] /= N
	m = max(ddf[:,1])
	ddf[:,1] /= m
	print("Finished..........<*_*>\n")
	return ave_dihedral_degree, ddf	

def get_data(temp):
	data = open('25-%d-normalizing.data'%temp)
	datalines = data.readlines()
	atoms = int(datalines[2].split()[0])
	xl = float(datalines[13].split()[0])
	xh = float(datalines[13].split()[1])
	xsize = xh - xl 
	yl = float(datalines[14].split()[0])
	yh = float(datalines[14].split()[1])
	ysize = yh - yl
	zl = float(datalines[15].split()[0])
	zh = float(datalines[15].split()[1])
	zsize = zh - zl  
	row = 0
	for dataline in datalines:
		if (dataline == 'Atoms # full\n'):
			line = row
		row += 1
	atom_number = np.zeros((atoms,2), dtype = int)
	atom_type = np.zeros((atoms,2), dtype = int)
	Xs = np.zeros((atoms,2), dtype = float)
	Ys = np.zeros((atoms,2), dtype = float)
	Zs = np.zeros((atoms,2), dtype = float)
	for i in range(atoms):
		a = datalines[line+2+i].split()
		atom_number[i,-1] = int(a[0])
		atom_type[i,-1] = int(a[2])
		Xs[i,-1] = float(a[4]) - xl
		Ys[i,-1] = float(a[5]) - yl
		Zs[i,-1] = float(a[6]) - zl
	sorted_indices = np.argsort(atom_number[:,-1])
	atom_number[:,-1] = np.take(atom_number[:,-1], sorted_indices)
	atom_type[:,-1] = np.take(atom_type[:,-1], sorted_indices)
	Xs[:,-1] = np.take(Xs[:,-1], sorted_indices)
	Ys[:,-1] = np.take(Ys[:,-1], sorted_indices)
	Zs[:,-1] = np.take(Zs[:,-1], sorted_indices)
	steps = 2
	n_mol, n_A, n_B = find_target(atoms, atom_type, temp)
	mass = np.zeros(atoms)
	for i in range(atoms):
		if atom_type[i,-1] == 1:
			mass[i] = 1.008
		elif atom_type[i,-1] == 2 or atom_type[i,-1] == 3:
			mass[i] = 12.011
		elif atom_type[i,-1] == 4 or atom_type[i,-1] == 5:
			mass[i] = 15.999
	n_atoms = len(n_A[:,0]) + len(n_B[:,0])
	n_mols = n_mol[-1]
	atom = np.zeros((n_atoms,7))
	atom[:,0] = np.arange(1, n_atoms+1)
	mass_A = np.sum(mass[n_A[0,:]]) + 1.008 * 5
	mass_B = np.sum(mass[n_B[0,:]]) + 1.008 * 3
	mass = [mass_A, mass_B]
	X_cg, Y_cg, Z_cg, type_cg, bond, angle = coarsegained(atom_type, Xs, Ys, Zs, atoms, steps, n_mol, n_A, n_B, xsize, ysize, zsize)
	atom[:,2] = type_cg
	atom[:,4] = X_cg[:,-1]
	atom[:,5] = Y_cg[:,-1]
	atom[:,6] = Z_cg[:,-1]
	for i in range(len(n_A[:,0])):
		atom[2*i,1] = n_mol[i]
		atom[2*i+1,1] = n_mol[i]
	n_bonds = n_atoms - n_mols
	n_angles = 3 * len(n_A[:,0]) - 4 * n_mols 
	return n_atoms, n_bonds, n_angles, xsize, ysize, zsize, atom, bond, angle, mass

def coarsegained(atom_type, Xs, Ys, Zs, atoms, steps, n_mol, n_A, n_B, xsize, ysize, zsize):
	n_A_len = len(n_A[:,0])
	n_atoms = len(n_A[:,0]) + len(n_B[:,0])
	n_mols = n_mol[-1]
	dp = int(n_atoms/n_mols)
	X_cg = np.zeros((n_atoms,steps))
	Y_cg = np.zeros((n_atoms,steps))
	Z_cg = np.zeros((n_atoms,steps))
	type_cg = np.zeros(n_atoms)
	mass = np.zeros((atoms))
	for i in range(atoms):
		if atom_type[i,-1] == 1:
			mass[i] = 1.008
		elif atom_type[i,-1] == 2 or atom_type[i,-1] == 3:
			mass[i] = 12.011
		elif atom_type[i,-1] == 4 or atom_type[i,-1] == 5:
			mass[i] = 15.999
		elif atom_type[i,-1] == 7:
			mass[i] = 18.998
	mass_A = np.sum(mass[n_A[0,:]])
	mass_B = np.sum(mass[n_B[0,:]])
	for i in range(n_A_len):
		type_cg[2*i] = 1
		type_cg[2*i+1] = 2
		for j in range(steps):
			X_cg[2*i,j] = fit.mass_center(mass[n_A[i,:]], Xs[n_A[i,:],j], xsize)
			Y_cg[2*i,j] = fit.mass_center(mass[n_A[i,:]], Ys[n_A[i,:],j], ysize)
			Z_cg[2*i,j] = fit.mass_center(mass[n_A[i,:]], Zs[n_A[i,:],j], zsize)
		for j in range(steps):
			X_cg[2*i+1,j] = fit.mass_center(mass[n_B[i,:]],Xs[n_B[i,:],j], xsize)
			Y_cg[2*i+1,j] = fit.mass_center(mass[n_B[i,:]],Ys[n_B[i,:],j], ysize)
			Z_cg[2*i+1,j] = fit.mass_center(mass[n_B[i,:]],Zs[n_B[i,:],j], zsize)
	n_bonds = n_atoms - n_mols
	bond = np.zeros((n_bonds,4))
	bond[:,0] = np.arange(1, n_bonds+1)
	j = 0
	for i in range(n_A_len):
		if (2*i+2) % dp == 0:
			j += 1
		else: 
			bond[i-j,1] = 1
			bond[i-j,2] = 2 * i + 1 
			bond[i-j,3] = 2 * i + 3
	for i in range(n_A_len):
		bond[n_A_len-n_mols+i,1] = 2
		bond[n_A_len-n_mols+i,2] = 2 * i + 1
		bond[n_A_len-n_mols+i,3] = 2 * i + 2
	n_angles = 3 * n_A_len - 4 * n_mols 
	angle = np.zeros((n_angles,5))
	angle[:,0] = np.arange(1, n_angles+1)
	j = 0
	for i in range(n_A_len):
		if (2*i+2) % dp == 0:
			j += 1
		elif (2*i+4) % dp == 0:
			j += 1
		else:
			angle[i-j,1] = 1
			angle[i-j,2] = 2 * i + 1
			angle[i-j,3] = 2 * i + 3
			angle[i-j,4] = 2 * i + 5
	j = 0
	for i in range(n_A_len):
		if (2*i+2) % dp == 0:
			j += 1
		else:
			angle[n_A_len-2*n_mols+i-j,1] = 2
			angle[n_A_len-2*n_mols+i-j,2] = 2 * i + 2
			angle[n_A_len-2*n_mols+i-j,3] = 2 * i + 1
			angle[n_A_len-2*n_mols+i-j,4] = 2 * i + 3
	j = 0
	for i in range(n_A_len):
		if (2*i) % dp == 0:
			j += 1
		else:
			angle[2*n_A_len-3*n_mols+i-j,1] = 2
			angle[2*n_A_len-3*n_mols+i-j,2] = 2 * i - 1
			angle[2*n_A_len-3*n_mols+i-j,3] = 2 * i + 1
			angle[2*n_A_len-3*n_mols+i-j,4] = 2 * i + 2
	return X_cg, Y_cg, Z_cg, type_cg, bond, angle 

def caclulate_distribution(temp):
	lammpstrj = open('25-%d-normalizing.lammpstrj'%(temp))
	atom_number, atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins,dr = readlammpstrj(lammpstrj)
	for j in range(0,steps,1):
		atom_number, atom_type, Xs, Ys, Zs = order_atoms(atom_number, atom_type, Xs, Ys, Zs, atoms, j)
	lammpstrj.close()
	n_mol, n_A, n_B = find_target(atoms, atom_type, temp)
	X_cg, Y_cg, Z_cg, type_cg, bond, angle = coarsegained(atom_type, Xs, Ys, Zs, atoms, steps, n_mol, n_A, n_B, xsize, ysize, zsize)
	bdf_file = open("Bond_distribution.txt","w")
	bdf = []
	bdf.append(0)
	ave_bondlength, bdf[0] = Bond_distribution(X_cg, Y_cg, Z_cg, xsize, ysize, zsize, steps, bond, bondtype=1)
	bdf_file.write("Bondtype:AA\nAverage length: %s\n"%(ave_bondlength))
	for i in range(len(bdf[0])):
		bdf_file.write('%f\t%f\n'%(bdf[0][i][0],bdf[0][i][1]))
	bdf_file.write("\n")
	bdf.append(0)
	ave_bondlength, bdf[1] = Bond_distribution(X_cg, Y_cg, Z_cg, xsize, ysize, zsize, steps, bond, bondtype=2)	
	bdf_file.write("Bondtype:AB\nAverage length: %s\n"%(ave_bondlength))
	for i in range(len(bdf[1])):
		bdf_file.write('%f\t%f\n'%(bdf[1][i][0],bdf[1][i][1]))
	bdf_file.write("\n")
	bdf_file.close()
	print("************Bond Calculation Finished***********")

	adf_file = open("Angle_distribution.txt","w")
	adf = []
	adf.append(0)
	ave_angle_degree, adf[0] = Angel_distribution(X_cg, Y_cg, Z_cg, xsize, ysize, zsize, steps, angle, angletype=1)
	adf_file.write("Angletype: AAA\nAverage degree: %s\n"%str(ave_angle_degree))
	for i in range(len(adf[0])):
		adf_file.write('%f\t%f\n'%(adf[0][i][0],adf[0][i][1]))
	adf_file.write("\n")
	adf.append(0)
	ave_angle_degree, adf[1] = Angel_distribution(X_cg, Y_cg, Z_cg, xsize, ysize, zsize, steps, angle, angletype=2)
	adf_file.write("Angletype: BAA\nAverage degree: %s\n"%str(ave_angle_degree))
	for i in range(len(adf[1])):
		adf_file.write('%f\t%f\n'%(adf[1][i][0],adf[1][i][1]))
	adf_file.close()
	print("************Angle Calculation Finished***********")	
	print("**********Start calculating RDF**********\n")
	rdf = []
	rdf_file = open("rdf.txt","w")
	r = dr*np.linspace(1e-14,NumOfBins,NumOfBins+1)
	rdf.append(RDF(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, type_cg, bond, angle, target_1=1, target_2=1))
	rdf_file.write("Atom 1: A\tAtom 2: A\n")
	for k in range(int(NumOfBins)):
		rdf_file.write('%f\t%f\n'%(r[k],rdf[0][k]))
	rdf_file.write("\n")
	rdf.append(RDF(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, type_cg, bond, angle, target_1=1, target_2=2))		
	rdf_file.write("Atom 1: A\tAtom 2: B\n")
	for k in range(int(NumOfBins)):
		rdf_file.write('%f\t%f\n'%(r[k],rdf[1][k]))
	rdf_file.write("\n")
	rdf.append(RDF(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, type_cg, bond, angle, target_1=2, target_2=2))	
	rdf_file.write("Atom 1: B\tAtom 2: B\n")	
	for k in range(int(NumOfBins)):
		rdf_file.write('%f\t%f\n'%(r[k],rdf[2][k]))
	rdf_file.write("\n")
	print("************RDF Calculation Finished***********")
	return r, rdf, bdf, adf




			




        
    