import numpy as np
import pandas as pd
import math

def readlammpstrj(lammpstrj):
	lammpstrj_lines = lammpstrj.readlines()
	steps = 0
	line_2 = []
	row = 0
	for line_1 in lammpstrj_lines:
		row = row + 1
		if line_1 == 'ITEM: NUMBER OF ATOMS\n':
			atoms = int(lammpstrj_lines[row].rstrip())   #记录体系中的原子数目
			steps = steps + 1                                                #记录文件中有多少步(抽样时间点)
			line_2.append(row - 2)                                #记录每步在文件中的行数 
	n = 0
	atom_number = np.zeros((atoms, steps))
	atom_type = np.zeros((atoms, steps))
	Xs = np.zeros((atoms, steps))
	Ys = np.zeros((atoms, steps))
	Zs = np.zeros((atoms, steps))
	xsize = np.zeros(steps)
	ysize = np.zeros(steps)
	zsize = np.zeros(steps)
	V = []
	NumOfBins = 100 #定义分割数目
	row = 1
	while n < steps: #对于每一步
		row = line_2[n]-1
		for line_1 in lammpstrj_lines:
			row = row + 1
			if n < steps - 1:
				if ((row <= line_2[n+1]) and (row >= line_2[n])):  #在每一步的行中搜索
					if line_1 == 'ITEM: BOX BOUNDS pp pp pp\n':
						bx = lammpstrj_lines[row].split()    #记录盒子尺寸
						by = lammpstrj_lines[row+1].split()
						bz = lammpstrj_lines[row+2].split()
						xsize[n] = (float(bx[1])-float(bx[0]))
						ysize[n] = (float(by[1])-float(by[0]))
						zsize[n] = (float(bz[1])-float(bz[0]))        

					if line_1 == 'ITEM: ATOMS id mol type q x y z ix iy iz\n':  #定位坐标初始行
						i = 0
						while i <= (atoms - 1):
							a = lammpstrj_lines[row + i].split()
							atom_number[i][n] = int(a[0])  	#记录原子编号
							atom_type[i][n] = int(a[2])     #记录原子种类
							Xs[i][n] = float(a[4])          #记录原子位置
							Ys[i][n] = float(a[5])          
							Zs[i][n] = float(a[6])          
							i = i + 1
			else:
				if ((row >= line_2[n])and (row<=len(lammpstrj_lines))):  #在每一步的行中搜索
					if line_1 == 'ITEM: BOX BOUNDS pp pp pp\n':
						bx = lammpstrj_lines[row].split()    #记录盒子尺寸
						by = lammpstrj_lines[row+1].split()
						bz = lammpstrj_lines[row+2].split()
						xsize[n] = (float(bx[1])-float(bx[0]))
						ysize[n] = (float(by[1])-float(by[0]))
						zsize[n] = (float(bz[1])-float(bz[0]))

					if line_1 == 'ITEM: ATOMS id mol type q x y z ix iy iz\n':  #定位坐标初始行
						i = 0
						while i <= (atoms - 1): 
							a = lammpstrj_lines[row + i].split()
							atom_number[i][n] = int(a[0])  	#记录原子编号
							atom_type[i][n] = int(a[2])     #记录原子种类
							Xs[i][n] = float(a[4])          #记录原子位置
							Ys[i][n] = float(a[5])          
							Zs[i][n] = float(a[6])          
							i = i + 1
		n = n + 1
	for n in range(0,steps,1):
		for i in range(0,atoms,1):
			if Xs[i][n] >= xsize[n]:
				Xs[i][n] -= xsize[n]
			elif Xs[i][n] <= 0:
				Xs[i][n] += xsize[n]
			if Ys[i][n] >= ysize[n]:
				Ys[i][n] -= ysize[n]
			elif Ys[i][n] <= 0:
				Ys[i][n] += ysize[n]
			if Zs[i][n] >= zsize[n]:
				Zs[i][n] -= zsize[n]
			elif Zs[i][n] <= 0:
				Zs[i][n] += zsize[n]
	Volumn = xsize * ysize * zsize
	rsize = 32
	dr = 0.5 * rsize / float(NumOfBins) #定义dr
	print("\n*************Calculating Radius Distribution Function************\n")
	print("Number of beads in calculation =  " + str(atoms) +", Number of Bead Type in calculation = " + str(np.max(atom_type)) +"\n")
	print("Point Number for discreting r =" +str(NumOfBins) +", Number of step for calculation = " +str(steps)+"\n")
	print("\t End of read data ........\n")
	return atom_number, atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr

def order_atoms(atom_number, atom_type, Xs, Ys, Zs, atoms, j):  #冒泡排序
	print("***************Order!***************\n")
	sorted_indices = np.argsort(atom_number[:,j])
	atom_number[:,j] = np.take(atom_number[:,j], sorted_indices)
	atom_type[:,j] = np.take(atom_type[:,j], sorted_indices)
	Xs[:,j] = np.take(Xs[:,j], sorted_indices)
	Ys[:,j] = np.take(Ys[:,j], sorted_indices)
	Zs[:,j] = np.take(Zs[:,j], sorted_indices)
	print("*****Order Successfully!*****\n")
	return atom_number, atom_type, Xs, Ys, Zs

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
	if not (isinstance(window_size, int) and isinstance(order, int)):#判断是否为int类型
		raise ValueError("window_size and order must be of type int")
	if window_size % 2 != 1 or window_size < 1:
		raise TypeError('window_size must be a positive odd number')
	if window_size < order + 2:
		raise TypeError('window_size is too small for the polynomials order')
	order_range = range(order+1)
	half_window = (window_size - 1) // 2
	b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window +1)])
	m = np.linalg.pinv(b).A[deriv] * rate**deriv * math.factorial(deriv)
	firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
	lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
	y = np.concatenate((firstvals, y, lastvals))
	return np.convolve(m[::-1], y, mode='valid')

def RDF(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, atom_types, target_1, target_2):  #计算径向分布函数gr
	
	NumOfBins = int(NumOfBins)
	rdf = np.zeros((NumOfBins+1, steps))
	r_space = dr*np.linspace(1e-14,NumOfBins,NumOfBins+1)
	print("START..........<*_*>\n")
	print("Atom 1: %d\tAtom 2: %d\n"%(target_1,target_2))	
	for j in range(steps): #对于每步
		atom_type = atom_types[:,j]
		type_1 = atom_type == target_1
		type_2 = atom_type == target_2
		N1 = np.sum(type_1)
		N2 = np.sum(type_2)
		print("N1=%d\tN2=%d\n"%(N1,N2))
		Xs_1 = np.array([X_cg[:,j][type_1]])
		Xs_2 = np.array([X_cg[:,j][type_2]])
		Ys_1 = np.array([Y_cg[:,j][type_1]])
		Ys_2 = np.array([Y_cg[:,j][type_2]])
		Zs_1 = np.array([Z_cg[:,j][type_1]])
		Zs_2 = np.array([Z_cg[:,j][type_2]])
		dx = np.abs(np.subtract.outer(Xs_1,Xs_2).flatten())
		dx = np.minimum(dx, np.abs(xsize[j] - dx))  # 距离不能超过盒子尺寸一半
		dy = np.abs(np.subtract.outer(Ys_1,Ys_2).flatten())
		dy = np.minimum(dy, np.abs(ysize[j] - dy))
		dz = np.abs(np.subtract.outer(Zs_1,Zs_2).flatten())
		dz = np.minimum(dz, np.abs(zsize[j] - dz))
		coords = np.array([dx, dy, dz])
		r_temp = np.sqrt(np.sum(coords ** 2, axis=0))
		r_temp_bin = np.digitize(r_temp / dr, np.linspace(1e-14, NumOfBins, NumOfBins+1))#分配至对数空间的索引
		r_temp_bin[np.where(r_temp_bin > NumOfBins)] = 0
		np.add.at(rdf, (r_temp_bin, j), 1)
		'''
		indices = np.where(rdf[:,j] < 1)[0]
		if len(indices) > 0:
			r_cut = indices[-1]
		else:
			r_cut = 0
		rdf[0:r_cut+1, j] = 0
		rdf[r_cut+1:-1, j] = savitzky_golay(rdf[r_cut+1:-1, j], 9, 2, deriv=0, rate=1)
		'''
	for r in range(NumOfBins): 
		if r>=1:
			r_u = float(r_space[r])
			r_b = float(r_space[r-1])
			rdf[r,:] /= ((N1 * N2) * (4 / 3) * math.pi * (r_u ** 3 - r_b ** 3) / Volumn[j])
		rdf[0,:] = 0 
	print("Finished..........<*_*>\n")
	return rdf

def RDF_special(X_cg, Y_cg, Z_cg, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr,special_atom_1, special_atom_2, j):
	r_space = dr*np.linspace(1e-14,NumOfBins,NumOfBins+1)
	rdf = np.zeros(NumOfBins+1)
	rdf_ave = []
	print("START Calculate Special_bonds..........<*_*>\n")
	coords_1 = pd.DataFrame({'x': X_cg[special_atom_1,j], 'y': Y_cg[special_atom_1,j], 'z': Z_cg[special_atom_1,j]})# 求原子对之间距离
	coords_2 = pd.DataFrame({'x': X_cg[special_atom_2,j], 'y': Y_cg[special_atom_2,j], 'z': Z_cg[special_atom_2,j]})# 利用pandas数据框实现行间的差值计算
	coords_1 = coords_1.reindex(coords_2.index)
	dx = np.abs(coords_1.x - coords_2.x.values)# 用广播操作代替循环
	dy = np.abs(coords_1.y - coords_2.y.values)
	dz = np.abs(coords_1.z - coords_2.z.values)
	dx = np.minimum(dx, np.abs(xsize - dx))  # 距离不能超过盒子尺寸一半
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

def find_peaks(r, rdf, min_dist=5):
	from scipy.signal import argrelextrema
	maxima_idx = argrelextrema(rdf, np.greater)[0]
	if len(maxima_idx) < 2:
		return []
	# 按峰值大小排序
	peaks = sorted([(r[i], rdf[i]) for i in maxima_idx], key=lambda x: x[1], reverse=True)
	# 筛选出前两个且距离间隔足够大的峰（避免同一峰附近的假峰）
	valid_peaks = []
	for p in peaks:
		if not valid_peaks or abs(p[0] - valid_peaks[-1][0]) > r[min_dist]:
			valid_peaks.append(p)
		if len(valid_peaks) == 2:
			break
	return valid_peaks


def caclulate_distribution():
	with open('crosslink.lammpstrj') as lammpstrj:
		atom_number, atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins,dr = readlammpstrj(lammpstrj)
	for j in range(0,steps,1):
		atom_number, atom_type, Xs, Ys, Zs = order_atoms(atom_number, atom_type, Xs, Ys, Zs, atoms, j)
	rdf_N_O = RDF(Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, atom_type, target_1=5, target_2=6)
	rdf_N_C = RDF(Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, atom_type, target_1=5, target_2=3)
	rdf_O_C = RDF(Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins, dr, atom_type, target_1=6, target_2=3)
	r = dr * np.linspace(1e-14, NumOfBins, NumOfBins+1)
	with open('rdf_peaks.txt', 'w') as peaks_file:
		for j in range(steps):
			peaks_N_O = find_peaks(r, rdf_N_O[:,j], min_dist=5)
			peaks_N_C = find_peaks(r, rdf_N_C[:,j], min_dist=5)
			peaks_O_C = find_peaks(r, rdf_O_C[:,j], min_dist=5)
			with open('rdf_N_O_%d.txt'%j, 'w') as file:
				for i in range(len(r)):
					file.write('%f\t%f\n'%(r[i], rdf_N_O[i,j]))
			with open('rdf_N_C_%d.txt'%j, 'w') as file:
				for i in range(len(r)):
					file.write('%f\t%f\n'%(r[i], rdf_N_C[i,j]))
			with open('rdf_O_C_%d.txt'%j, 'w') as file:
				for i in range(len(r)):
					file.write('%f\t%f\n'%(r[i], rdf_O_C[i,j]))
			peaks_file.write('%d\t%f\t%f\t%f\t%f\t%f\t%f\n'%(j, peaks_N_O[0][1], peaks_N_O[1][1], peaks_N_C[0][1], peaks_N_C[1][1], peaks_O_C[0][1], peaks_O_C[1][1]))

if __name__ == '__main__':
	caclulate_distribution()