from multiprocessing import Pool, cpu_count
#from numba import njit, prange
from scipy.spatial import KDTree
import time
import numpy as np
import pandas as pd
import os
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

def free_volumn_parallel(args):
	frame, n_atoms, atom_type, radiu, Xs, Ys, Zs, Lx, Ly, Lz, probe_radius, n_samples = args
	"""
	自由体积计算核心函数
	参数说明：
	probe_radius - 探针半径（单位：A）
	n_samples    - 每个网格的采样点数
	"""
	lx = Lx[frame]
	ly = Ly[frame]
	lz = Lz[frame]
	n_grid_x = max(1, int(lx // 10))#创建当前帧网格
	n_grid_y = max(1, int(ly // 10))
	n_grid_z = max(1, int(lz // 10))
	#print('%d %d %d\n'%(n_grid_x,n_grid_y,n_grid_z))
	grid_x = np.linspace(0, lx, n_grid_x+1)
	grid_y = np.linspace(0, ly, n_grid_y+1)
	grid_z = np.linspace(0, lz, n_grid_z+1)
	#将当前帧的原子分配至网格
	grid_atoms = assign_grid(Xs, Ys, Zs, n_atoms, frame, 
		grid_x, grid_y, grid_z, 
		lx, ly, lz, 
		n_grid_x, n_grid_y, n_grid_z)
	free_count = 0
	# 遍历所有网格
	for i in range(n_grid_x):
		for j in range(n_grid_y):
			for k in range(n_grid_z):
				# 生成当前网格的采样点
				x_samples = np.random.uniform(grid_x[i], grid_x[i+1], n_samples)
				y_samples = np.random.uniform(grid_y[j], grid_y[j+1], n_samples)
				z_samples = np.random.uniform(grid_z[k], grid_z[k+1], n_samples)
				# 收集相邻26个网格的原子索引
				neighbor_indices = []
				for di in [-1, 0, 1]:
					for dj in [-1, 0, 1]:
						for dk in [-1, 0, 1]:
							ii = (i + di) % n_grid_x
							jj = (j + dj) % n_grid_y
							kk = (k + dk) % n_grid_z
							neighbor_indices.extend(grid_atoms[ii][jj][kk])
				if (not neighbor_indices) or (len(neighbor_indices) == 0):
					free_count += n_samples
					continue
				# 向量化计算
				x_atoms = Xs[neighbor_indices, frame].astype(np.float64)
				y_atoms = Ys[neighbor_indices, frame].astype(np.float64)
				z_atoms = Zs[neighbor_indices, frame].astype(np.float64)
				types = atom_type[neighbor_indices, frame].astype(np.int)
				radii = radiu[types-1].astype(np.float64)  # 原子类型到半径映射
				collision_flags = compute_collisions(x_atoms, y_atoms, z_atoms, radii,
					x_samples, y_samples, z_samples, lx, ly, lz, probe_radius)
				free_count += np.sum(~collision_flags)
	# 计算自由体积比例
	total_samples = n_grid_x * n_grid_y * n_grid_z * n_samples
	free_frac = free_count / total_samples
	return free_frac

def assign_grid(Xs, Ys, Zs, n_monomers, frame, grid_x, grid_y, grid_z, l_x, l_y, l_z, n_grid_x, n_grid_y, n_grid_z):
	dx = grid_x[1] - grid_x[0]
	dy = grid_y[1] - grid_y[0]
	dz = grid_z[1] - grid_z[0]
	# 提取当前帧的坐标并处理周期性边界条件
	x = np.mod(Xs[:,frame], l_x)  # 使用np.mod向量化处理
	y = np.mod(Ys[:,frame], l_y)
	z = np.mod(Zs[:,frame], l_z)
	# 向量化计算网格索引 (i,j,k)
	i = np.floor_divide(x, dx).astype(int)
	j = np.floor_divide(y, dy).astype(int)
	k = np.floor_divide(z, dz).astype(int)
	# 使用np.clip限制索引范围
	i = np.clip(i, 0, n_grid_x-1)
	j = np.clip(j, 0, n_grid_y-1)
	k = np.clip(k, 0, n_grid_z-1)
	# 初始化三维网格容器（使用更高效的dtype=object）
	grid_atoms = np.empty((n_grid_x, n_grid_y, n_grid_z), dtype=object)
	for ii in range(n_grid_x):
		for jj in range(n_grid_y):
			for kk in range(n_grid_z):
				grid_atoms[ii][jj][kk] = []  # 每个网格初始化为空列表
	# 批量获取所有原子的网格坐标并填充到对应网格
	ijk = np.stack([i, j, k], axis=1)
	unique_ijk, inverse = np.unique(ijk, axis=0, return_inverse=True)
	# 使用分组技巧快速填充原子索引
	for idx, (ii, jj, kk) in enumerate(unique_ijk):
		mask = (inverse == idx)  # 属于当前网格的原子掩码
		atom_indices = np.where(mask)[0].tolist()  # 直接获取索引列表
		grid_atoms[ii][jj][kk].extend(atom_indices)
	return grid_atoms

#@njit(parallel=True)
def compute_collisions(x_atoms, y_atoms, z_atoms, radii, x_samples, y_samples, z_samples, lx, ly, lz, probe_radius):
	n_atoms = len(x_atoms)
	n_samples = len(x_samples)
	collision_flags = np.zeros(n_samples, dtype=np.bool_)
	sum_radii = radii + probe_radius
	sum_radii_sq = sum_radii**2
	# 并行计算距离矩阵
	for j in range(n_samples):
		x = x_samples[j]
		y = y_samples[j]
		z = z_samples[j]
		for i in range(n_atoms):
			dx = x_atoms[i] - x
			dx -= lx * np.round(dx / lx)
			dy = y_atoms[i] - y
			dy -= ly * np.round(dy / ly)
			dz = z_atoms[i] - z
			dz -= lz * np.round(dz / lz)
			dist_sq = dx**2 + dy**2 + dz**2
			if dist_sq < sum_radii_sq[i]:
				collision_flags[j] = True
				break  # 检测到碰撞即跳出
	return collision_flags



if __name__ == '__main__':
	with open('crosslink.lammpstrj') as lammpstrj:
		atom_number, atom_type, Xs, Ys, Zs, atoms, steps, Volumn, xsize, ysize, zsize, NumOfBins,dr = readlammpstrj(lammpstrj)
	for j in range(0,steps,1):
		atom_number, atom_type, Xs, Ys, Zs = order_atoms(atom_number, atom_type, Xs, Ys, Zs, atoms, j)
	diameters = np.array([2.576, 2.576, 3.434, 3.434, 3.277, 3.118, 3.118, 3.677])
	radius = diameters/2
	probe_radius = 1.5
	n_samples = 1000
	with Pool(cpu_count()) as p:
		args_list = [(frame, atoms, atom_type, radius, Xs, Ys, Zs, xsize, ysize, zsize, probe_radius, n_samples) for frame in range(steps)]
		results = p.map(free_volumn_parallel, args_list)
	with open('ffv.txt','w') as file:
		for i in range(steps):
			file.write('%f\n'%results[i])
