from PyAPM import build_branch_copolymer
import numpy as np
import os
import shutil

if __name__ == "__main__":
	n_mol = 80
	n_seg = 1
	segn = [50]
	segt = [1]
	lx, ly, lz = 100, 100, 100
	n_branch = np.full((1,50),1,dtype=int)
	branchtype = np.full((1,50),2,dtype=int)
	temp = 323
	n_step = 1
	path = 'origin' 
	build_cg_block_copolymer(n_mol, n_seg, segn, segt, lx, ly, lz, n_branch, branchtype, temp, n_step, path)