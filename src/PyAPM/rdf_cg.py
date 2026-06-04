import numpy as np
import pandas as pd
import os
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
	NumOfBins = 100
	row = 1
	while n < steps:
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
	print("Number of beads in calculation =  " + str(atoms) +", Number of Bead Type in calculation = " + str(atom_type) +"\n")
	print("Point Number for discreting r =" +str(NumOfBins) +", Number of step for calculation = " +str(steps)+"\n")
	print("\t End of read data ........\n")
	return atom_number, atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr

def readdata(datalines):
	row = 0
	if_Bonds = False
	if_Angles = False
	if_Dihedrals = False
	if_Impropers = False
	for line in datalines:
		if line == "Bonds\n":
			if_Bonds = True
			line_Bonds = row
		elif line == "Angles\n":
			if_Angles = True
			line_Angles = row
		elif line == "Dihedrals\n":
			if_Dihedrals = True
			line_Dihedrals = row
		elif line == "Impropers":
			if_Impropers = True
			line_Impropers = row
		row += 1
	if if_Bonds == False:
		print("No Bond!\n")
		line_Bonds = None
	if if_Angles == False:
		print("No Angle!\n")
		line_Angles = None
	if if_Dihedrals == False:
		print("No Dihedrals!\n")
		line_Dihedrals = None
	if if_Impropers == False:
		print("No Impropers!\n")
		line_Impropers = None
	a = datalines[3].split()
	atomtypes = int(a[0])
	return if_Bonds, line_Bonds, if_Angles, line_Angles, if_Dihedrals, line_Dihedrals, if_Impropers, line_Impropers, atomtypes

def readbonds(datalines, line_Bonds):
	a = datalines[4].split()
	Bonds = int(a[0])
	a = datalines[5].split()
	Bondtypes = int(a[0])
	NumOfBonds = np.zeros(Bonds)
	TypeOfBonds = np.zeros(Bonds)
	Bondatom_1 = np.zeros(Bonds)
	Bondatom_2 = np.zeros(Bonds)
	for i in range(0,Bonds,1):
		a = datalines[line_Bonds+2+i].split()
		NumOfBonds[i] = int(a[0])
		TypeOfBonds[i] = int(a[1])
		Bondatom_1[i] = int(a[2])
		Bondatom_2[i] = int(a[3])
	print("Read Bonds Successfully!\n")
	return Bonds, Bondtypes, NumOfBonds.astype(int), TypeOfBonds.astype(int), Bondatom_1.astype(int), Bondatom_2.astype(int)

def readangles(datalines, line_Angles):
	a = datalines[6].split()
	Angles = int(a[0])
	a = datalines[7].split()
	Angletypes = int(a[0])
	NumOfAngles = np.zeros(Angles)
	TypeOfAngles = np.zeros(Angles)
	Angleatom_1 = np.zeros(Angles)
	Angleatom_2 = np.zeros(Angles)
	Angleatom_3 = np.zeros(Angles)
	for i in range(0,Angles,1):
		a = datalines[line_Angles+2+i].split()
		NumOfAngles[i] = int(a[0])
		TypeOfAngles[i] = int(a[1])
		Angleatom_1[i] = int(a[2])
		Angleatom_2[i] = int(a[3])
		Angleatom_3[i] = int(a[4])
	print("Read Angles Successfully!\n")
	return Angles, Angletypes, NumOfAngles.astype(int), TypeOfAngles.astype(int), Angleatom_1.astype(int), Angleatom_2.astype(int), Angleatom_3.astype(int)

def readdihedrals(data, line_Dihedrals):
	a = datalines[8].split()
	Dihedrals = int(a[0])
	a = datalines[9].split()
	Dihedraltypes = int(a[0])
	datalines = data.readlines()
	NumofDihedrals = np.zeros(Dihedrals)
	TypeofDihedrals = np.zeros(Dihedrals)
	Dihedralatom_1 = np.zeros(Dihedrals)
	Dihedralatom_2 = np.zeros(Dihedrals)
	Dihedralatom_3 = np.zeros(Dihedrals)
	Dihedralatom_4 = np.zeros(Dihedrals)
	for i in range(0,Dihedrals,1):
		a = datalines[line_Dihedrals+2+i].split()
		NumofDihedrals[i] = int(a[0])
		TypeofDihedrals = int(a[1])
		Dihedralatom_1 = int(a[2])
		Dihedralatom_2 = int(a[3])
		Dihedralatom_3 = int(a[4])
		Dihedralatom_4 = int(a[5])
	print("Read Dihedrals Successfully!\n")
	return Dihedrals, Dihedraltypes, NumofDihedrals.astype(int), TypeofDihedrals.astype(int), Dihedralatom_1.astype(int), Dihedralatom_2.astype(int), Dihedralatom_3.astype(int), Dihedralatom_4.astype(int)

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


def RDF(atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, type_target_1, type_target_2, bond_atom_1, bond_atom_2, angle_atom_1, angle_atom_2, TypeOfBonds, TypeOfAngles):  #计算径向分布函数gr
	print("Start!................<*_*>\n")
	r_space = dr*np.linspace(1e-14,NumOfBins,NumOfBins+1)
	rdf = np.zeros((NumOfBins+1, steps))
	rdf_ave = np.zeros(NumOfBins+1)
	atom_1 = atom_type == type_target_1
	atom_2 = atom_type == type_target_2
	if (type_target_1 == 1 and type_target_2 == 1):
		bond_atom_1 = bond_atom_1[np.where(TypeOfBonds==1)]
		bond_atom_1 = np.subtract(bond_atom_1, 1)
		bond_atom_2 = bond_atom_2[np.where(TypeOfBonds==1)]
		bond_atom_2 = np.subtract(bond_atom_2, 1)
		angle_atom_1 = angle_atom_1[np.where(TypeOfAngles==1)]
		angle_atom_1 = np.subtract(angle_atom_1, 1)
		angle_atom_2 = angle_atom_2[np.where(TypeOfAngles==1)]
		angle_atom_2 = np.subtract(angle_atom_2, 1)
	for j in range(steps): #对于每步
		N1 = np.sum(atom_1[:,j])
		N2 = np.sum(atom_2[:,j])
		dx = np.abs(np.subtract.outer(Xs[:,j][atom_1[:,j]],Xs[:,j][atom_2[:,j]]).flatten())
		dy = np.abs(np.subtract.outer(Ys[:,j][atom_1[:,j]],Ys[:,j][atom_2[:,j]]).flatten())
		dz = np.abs(np.subtract.outer(Zs[:,j][atom_1[:,j]],Zs[:,j][atom_2[:,j]]).flatten())
		dx = np.minimum(dx, np.abs(xsize - dx))
		dy = np.minimum(dy, np.abs(ysize - dy))
		dz = np.minimum(dz, np.abs(zsize - dz))
		coords = np.array([dx, dy, dz])
		r_temp = np.sqrt(np.sum(coords ** 2, axis=0))
		r_temp_bin = np.digitize(r_temp / dr,np.linspace(1e-14,NumOfBins,NumOfBins+1))
		r_temp_bin[np.where(r_temp_bin > NumOfBins)] = 0
		np.add.at(rdf, (r_temp_bin, j), 1)
		if (type_target_1 == 1 and type_target_2 == 1):
			rdf_12=RDF_special(Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, bond_atom_1, bond_atom_2, j)
			rdf_13=RDF_special(Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, angle_atom_1, angle_atom_2, j)
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

def RDF_special(Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr,special_atom_1, special_atom_2, j):
	r_space = dr*np.linspace(1e-14,NumOfBins,NumOfBins+1)
	rdf = np.zeros(NumOfBins+1)
	rdf_ave = []
	print("START Calculate Special_bonds..........<*_*>\n")
	coords_1 = pd.DataFrame({'x': Xs[special_atom_1,j], 'y': Ys[special_atom_1,j], 'z': Zs[special_atom_1,j]})
	coords_2 = pd.DataFrame({'x': Xs[special_atom_2,j], 'y': Ys[special_atom_2,j], 'z': Zs[special_atom_2,j]})
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

def Bond_distribution(Bonds, bondtype_target, TypeOfBonds, Bondatom_1, Bondatom_2, Xs, Ys, Zs, xsize, ysize, zsize, steps, NumOfBins):
	print("Start Calculating Bond "+str(bondtype_target)+"\n")
	NumOfBins = 100
	bondlength = []
	N = 0
	for i in range(Bonds):
		if TypeOfBonds[i] == bondtype_target:
			for j in range(steps):
				N+=1
				dx = np.abs(Xs[Bondatom_1[i]-1,j]-Xs[Bondatom_2[i]-1,j])
				x_dist = np.minimum(dx, np.abs(xsize - dx))
				dy = np.abs(Ys[Bondatom_1[i]-1,j]-Ys[Bondatom_2[i]-1,j])
				y_dist = np.minimum(dy, np.abs(ysize - dy))
				dz = np.abs(Zs[Bondatom_1[i]-1,j]-Zs[Bondatom_2[i]-1,j])
				z_dist = np.minimum(dz, np.abs(zsize - dz))
				bondlength.append(np.sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist))
	ave_bondlength = np.sum(bondlength)/N
	l_u = 4.5
	l_b = 2
	bdf = np.zeros((NumOfBins,2))
	l_temp_bin = np.digitize(bondlength,np.linspace(l_b, l_u, NumOfBins))
	np.add.at(bdf, (l_temp_bin - 1, 1), 1)
	bdf[:,0] = np.linspace(l_b, l_u, NumOfBins)
	bdf[:,1] /= N
	bdf[:,1] /= bdf[:,0]**2
	bdf[-1,1] = 0
	bdf[-2,1] = 0
	bdf[:,1] = fit.savitzky_golay(bdf[:,1], 9, 2, deriv=0, rate=1)
	bdf[:,1] = np.clip(bdf[:,1], 0, None)
	m = max(bdf[:,1])
	bdf[:,1] /= m
	print("Finished..........<*_*>\n")
	return ave_bondlength, bdf

def Angel_distribution(Angles, angletype_target, TypeOfAngles, Angleatom_1, Angleatom_2, Angleatom_3, Xs, Ys, Zs, xsize, ysize, zsize, steps):
	print("Start Calculating Angle "+str(angletype_target)+"\n")
	N = 0
	angle_degree = []
	for i in range(Angles):
		if TypeOfAngles[i] == angletype_target:
			for j in range(steps):
				N+=1
				dx_1 = Xs[Angleatom_2[i]-1,j] - Xs[Angleatom_1[i]-1,j]
				dx_1 = np.mod(dx_1 + xsize/2, xsize) - xsize/2
				dy_1 = Ys[Angleatom_2[i]-1,j] - Ys[Angleatom_1[i]-1,j]
				dy_1 = np.mod(dy_1 + ysize/2, ysize) - ysize/2
				dz_1 = Zs[Angleatom_2[i]-1,j] - Zs[Angleatom_1[i]-1,j]
				dz_1 = np.mod(dz_1 + zsize/2, zsize) - zsize/2
				dx_2 = Xs[Angleatom_3[i]-1,j] - Xs[Angleatom_2[i]-1,j]
				dx_2 = np.mod(dx_2 + xsize/2, xsize) - xsize/2
				dy_2 = Ys[Angleatom_3[i]-1,j] - Ys[Angleatom_2[i]-1,j]
				dy_2 = np.mod(dy_2 + ysize/2, ysize) - ysize/2
				dz_2 = Zs[Angleatom_3[i]-1,j] - Zs[Angleatom_2[i]-1,j]
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
	adf[-3:,1] = 0
	adf[:,1] /= N
	adf[:,1] /= np.sin(np.deg2rad(adf[:,0]))
	adf[:,1] = fit.savitzky_golay(adf[:,1], 9, 2, deriv=0, rate=1)
	adf[:,1] = np.clip(adf[:,1], 0, None)
	m = np.max(adf[:,1])
	adf[:,1] /= m
	print("Finished..........<*_*>\n")
	return ave_angle_degree, adf

def Dihedral_distributation(Dihedrals, dihedraltype_target, TypeofDihedrals, Dihedralatom_1, Dihedralatom_2, Dihedralatom_3, Dihedralatom_4, Xs, Ys, Zs, xsize, ysize, zsize, steps):
	print("Start Calculating Dihedras "+str(dihedraltype_target)+"\n")
	N = 0
	dihedral_degree = []
	for i in range(Dihedrals):
		if TypeofDihedrals[i] == dihedraltype_target:
			for j in range(steps):
				N+=1
				dx_1 = Xs[Dihedralatom_2[i]-1,j] - Xs[Dihedralatom_1[i]-1,j]
				dx_1 = np.mod(dx_1 + xsize, xsize * 2) - xsize
				dy_1 = Ys[Dihedralatom_2[i]-1,j] - Ys[Dihedralatom_1[i]-1,j]
				dy_1 = np.mod(dy_1 + ysize, ysize * 2) - ysize
				dz_1 = Zs[Dihedralatom_2[i]-1,j] - Zs[Dihedralatom_1[i]-1,j]
				dz_1 = np.mod(dz_1 + zsize, zsize * 2) - zsize
				dx_2 = Xs[Dihedralatom_3[i]-1,j] - Xs[Dihedralatom_2[i]-1,j]
				dx_2 = np.mod(dx_2 + xsize, xsize * 2) - xsize
				dy_2 = Ys[Dihedralatom_3[i]-1,j] - Ys[Dihedralatom_2[i]-1,j]
				dy_2 = np.mod(dy_2 + ysize, ysize * 2) - ysize
				dz_2 = Zs[Dihedralatom_3[i]-1,j] - Zs[Dihedralatom_2[i]-1,j]
				dz_2 = np.mod(dz_2 + zsize, zsize * 2) - zsize
				dx_3 = Xs[Dihedralatom_4[i]-1,j] - Xs[Dihedralatom_3[i]-1,j]
				dx_3 = np.mod(dx_3 + xsize, xsize * 2) - xsize
				dy_3 = Ys[Dihedralatom_4[i]-1,j] - Ys[Dihedralatom_3[i]-1,j]
				dy_3 = np.mod(dy_3 + ysize, ysize * 2) - ysize
				dz_3 = Zs[Dihedralatom_4[i]-1,j] - Zs[Dihedralatom_3[i]-1,j]
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

def caclulate_distribution(temp):
	lammpstrj = open('%d-normalizing.lammpstrj'%(temp))
	atom_number, atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins,dr = readlammpstrj(lammpstrj)
	for j in range(0,steps,1):
		atom_number, atom_type, Xs, Ys, Zs = order_atoms(atom_number, atom_type, Xs, Ys, Zs, atoms, j)
	lammpstrj.close()
	data = open('%d-normalizing.data'%(temp))
	datalines = data.readlines()
	if_Bonds, line_Bonds, if_Angles, line_Angles, if_Dihedrals, line_Dihedrals, if_Impropers, line_Impropers, atomtypes = readdata(datalines)
	if if_Bonds == True:
		print("Bonds Reading......\n")
		Bonds, Bondtypes, NumOfBonds, TypeOfBonds, Bondatom_1, Bondatom_2 = readbonds(datalines, line_Bonds)
		bdf_file = open("%d-Bond_distribution.txt"%(temp),"w")
		bdf = []
		for j in range(1,Bondtypes+1):
			bdf.append(0)
			ave_bondlength, bdf[j-1] = Bond_distribution(Bonds, j, TypeOfBonds, Bondatom_1, Bondatom_2, Xs, Ys, Zs, xsize, ysize, zsize, steps, NumOfBins)
			bdf_file.write("Bondtype: %d\nAverage length: %f\n"%(j,ave_bondlength))
			for i in range(len(bdf[j-1])):
				bdf_file.write("%f\t%f\n"%(bdf[j-1][i][0],bdf[j-1][i][1]))
			bdf_file.write("\n")
		bdf_file.close()
		print("************Bond Calculation Finished***********")
	if if_Angles == True:
		print("Angles Reading......\n")
		Angles, Angletypes, NumOfAngles, TypeOfAngles, Angleatom_1, Angleatom_2, Angleatom_3 = readangles(datalines, line_Angles)
		adf_file = open("%d-Angle_distribution.txt"%(temp),"w")
		adf = []
		for j in range(1,Angletypes+1):
			adf.append(0)
			ave_angle_degree, adf[j-1] = Angel_distribution(Angles, j, TypeOfAngles, Angleatom_1, Angleatom_2, Angleatom_3, Xs, Ys, Zs, xsize, ysize, zsize, steps)
			adf_file.write("Angletype: %d\nAverage degree: %f\n"%(j,ave_angle_degree))
			for i in range(len(adf[j-1])):
				adf_file.write("%f\t%f\n"%(adf[j-1][i][0],adf[j-1][i][1]))
			adf_file.write("\n")
		adf_file.close()
		print("************Angle Calculation Finished***********")
	if if_Dihedrals == True:
		print("Dihedrals Reading......\n")
		Dihedrals, Dihedraltypes, NumofDihedrals, TypeofDihedrals, Dihedralatom_1, Dihedralatom_2, Dihedralatom_3, Dihedralatom_4 = readdihedrals(data,line_Dihedrals)
		ddf_file = open("%d-Dihedral_distribution.txt"%(temp),"w")
		ddf = []
		for j in range(1,Dihedraltypes+1):
			ddf.append(0)
			ave_dihedral_degree, ddf[j-1] = Dihedral_distributation(Dihedrals, j, TypeofDihedrals, Dihedralatom_1, Dihedralatom_2, Dihedralatom_3, Dihedralatom_4, Xs, Ys, Zs, xsize, ysize, zsize, steps)
			ddf_file.write("Dihedral type: %d\nAverage degree: %f\n"%(j,ave_dihedral_degree))
			for i in range(len(ddf[j-1])):
				ddf_file.write("%f\t%f\n"%(ddf[j-1][i][0],ddf[j-1][i][1]))
			ddf_file.write("\n")
		ddf_file.close()
		print("************Dihedral Calculation Finished***********")

	rdf = []
	print("**********Start calculating RDF**********\n")
	special_bonds = [0,0,1]
	r = np.zeros(NumOfBins)
	rdf_file = open("%d-rdf.txt"%(temp),"w")
	r = dr*np.linspace(1e-14,NumOfBins,NumOfBins+1)
	for i in range(1,atomtypes+1,1):
		for j in range(1,i+1,1):
			rdf.append(RDF(atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins,dr,i,j, Bondatom_1, Bondatom_2, Angleatom_1, Angleatom_3, TypeOfBonds, TypeOfAngles))
			rdf_file.write("Atom 1: %d\tAtom 2: %d\n"%(i,j))
			for k in range(int(NumOfBins)):
				rdf_file.write("%f\t%f\n"%(r[k],rdf[i+j-2][k]))
			rdf_file.write("\n")
	print("************RDF Calculation Finished***********")
	print("************Finish Tempterature=%dK***********"%(temp))
	return r, rdf, bdf, adf



			




        
    