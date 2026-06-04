import numpy as np
import math
import shutil
import os

class Chain:
	def __init__(self, n_mol, n_seg, segn, segt, lx, ly, lz, mass, bondtype, angletype, dihedraltype, bondcoeff, anglecoeff):
		self.n_mol = n_mol                                             #Total number of molecules
		self.n_seg = n_seg                                             #Number of blocks per molecule
		self.segn = np.array(segn)                                     #Number of atoms per block
		self.segt = np.array(segt)                                     #Atom type name (number) for each block
		self.lx = lx                                                   #Box boundary lengths
		self.ly = ly
		self.lz = lz
		self.mass = np.array(mass)                                     #Mass of each atom type
		self.bondtype = np.array(bondtype)                             #Bond types: [[1 bond type, bond atom 1, bond atom 2], [2 bond type, bond atom 1, bond atom 2], ...]
		self.angletype = np.array(angletype)                           #Angle types: [[1 angle type, angle atom 1, angle atom 2, angle atom 3], [2 angle type, angle atom 1, angle atom 2, angle atom 3], ...]
		self.dihedraltype = np.array(dihedraltype)                     #Dihedral types: [[1 dihedral type, dihedral atom 1, dihedral atom 2, dihedral atom 3, dihedral atom 4], [2 dihedral type, dihedral atom 1, dihedral atom 2, dihedral atom 3, dihedral atom 4], ...]
		self.bondcoeff = bondcoeff			                           #Bond parameters: [[1 bond amplitude, 1 bond variance, 1 bond mean], [...] , ...]
		self.anglecoeff = anglecoeff			                       #Angle parameters
		self.n_atom = self.n_mol * np.sum(self.segn)                   #Total number of atoms
		self.n_atomtype = np.max(self.segt)                            #Total number of atom types
		self.n_bond = self.n_mol * (np.sum(self.segn) - 1)             #Total number of bonds
		self.n_angle = self.n_mol * (np.sum(self.segn) - 2)            #Total number of angles
		self.n_dihedral = self.n_mol * (np.sum(self.segn) - 3)         #Total number of dihedrals
		self.atom = np.zeros((self.n_atom, 7))                         #Atom information matrix: columns are atom index, molecule index, atom type, atom charge, atom coordinates x, y, z
		self.bond = np.zeros((self.n_bond, 4), dtype = int)            #Bond information matrix: columns are bond index, bond type, bond atom 1, bond atom 2
		self.angle = np.zeros((self.n_angle, 5), dtype = int)          #Angle information matrix: columns are angle index, angle type, angle atom 1, angle atom 2, angle atom 3
		self.dihedral = np.zeros((self.n_dihedral, 6), dtype = int)    #Dihedral information matrix: columns are dihedral index, dihedral type, dihedral atom 1, dihedral atom 2, dihedral atom 3, dihedral atom 4
		self.file = open('system.data','w')                            #Output file name
	
	def	get_atom(self):
		na = 0
		nm = 1
		self.atom[:,0] = np.arange(1, self.n_atom + 1)
		for i in range(self.n_mol):
			for j in range(self.n_seg):
				for k in range(self.segn[j]):
					self.atom[na,2] = self.segt[j]
					self.atom[na,1] = nm
					self.atom[na,3] = 0.0   #Atomic charge
					if k == 0 and j == 0:   #If it is the first atom on each chain
						self.atom[na,4] = np.random.uniform(0, self.lx, size = 1)
						self.atom[na,5] = np.random.uniform(0, self.ly, size = 1)
						self.atom[na,6] = np.random.uniform(0, self.lz, size = 1)
						ox1 = self.atom[na,4]
						oy1 = self.atom[na,5]
						oz1 = self.atom[na,6]
					elif (k == 1 and j == 0) or (self.segn[0] == 1 and k == 0 and j == 1):  #If it is the second atom on each chain
						for l in range(int(np.max(self.bondtype[:,0]))):
							if ((self.atom[na,2] == self.bondtype[l,1] and self.atom[na-1,2] == self.bondtype[l,2]) or (self.atom[na,2] == self.bondtype[l,2] and self.atom[na-1,2] == self.bondtype[l,1])):
								peak = int(self.bondcoeff[self.bondtype[l,0]-1][0])
								bondlength = np.zeros(peak)
								amplitude = np.zeros(peak)
								for m in range(peak):
									mean = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+3])
									std_dev = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+2])
									bondlength[m] = np.random.normal(mean, std_dev, size = 1)
									amplitude[m] = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+1])
								amplitude = amplitude/np.sum(amplitude)
								chosen_index = np.random.choice(len(bondlength), p=amplitude)
								chosen_bondlength = bondlength[chosen_index] 
						theta = np.random.uniform(0, 2*np.pi, size = 1)#Horizontal azimuth angle
						phi = np.random.uniform(0, np.pi, size = 1)#Vertical polar angle
						self.atom[na,4] = ox1 + chosen_bondlength * np.sin(phi) * np.cos(theta)   #Randomly generate the position of the second atom on a spherical surface
						self.atom[na,5] = oy1 + chosen_bondlength * np.sin(phi) * np.sin(theta)
						self.atom[na,6] = oz1 + chosen_bondlength * np.cos(phi)
						ox2 = ox1    #Record the positions of the first two atoms
						oy2 = oy1
						oz2 = oz1
						ox1 = self.atom[na,4]
						oy1 = self.atom[na,5]
						oz1 = self.atom[na,6]
					else:
						for l in range(int(np.max(self.bondtype[:,0]))): #Select appropriate bond and bond angle parameters
							if ((self.atom[na,2] == self.bondtype[l,1] and self.atom[na-1,2] == self.bondtype[l,2]) or (self.atom[na,2] == self.bondtype[l,2] and self.atom[na-1,2] == self.bondtype[l,1])):
								peak = int(self.bondcoeff[self.bondtype[l,0]-1][0])
								bondlength = np.zeros(peak)
								amplitude = np.zeros(peak)
								for m in range(peak):
									mean = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+3])
									std_dev = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+2])
									bondlength[m] = np.random.normal(mean, std_dev, size = 1)
									amplitude[m] = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+1])
								amplitude = amplitude/np.sum(amplitude)
								chosen_index = np.random.choice(len(bondlength), p=amplitude)
								chosen_bondlength = bondlength[chosen_index]
						for l in range(int(np.max(self.angletype[:,0]))):
							if (self.atom[na-1,2] == self.angletype[l,2] and ((self.atom[na-2,2] == self.angletype[l,1] and self.atom[na,2] == self.angletype[l,3]) or (self.atom[na-2,2] == self.angletype[l,3] and self.atom[na,2] == self.angletype[l,1]))):
								peak = int(self.anglecoeff[self.angletype[l,0]-1][0])
								angledegree = np.zeros(peak)
								amplitude = np.zeros(peak)
								for m in range(peak):
									mean = float(self.anglecoeff[self.angletype[l,0]-1][3*m+3])
									std_dev = float(self.anglecoeff[self.angletype[l,0]-1][3*m+2])
									angledegree[m] =  np.random.normal(mean, std_dev, size = 1)
									amplitude[m] = float(self.anglecoeff[self.angletype[l,0]-1][3*m+1])
								amplitude = amplitude/np.sum(amplitude)
								chosen_index = np.random.choice(len(angledegree), p=amplitude)
								theta = np.radians(angledegree[chosen_index])
						vab_unit = unit_vector(ox2,oy2,oz2,ox1,oy1,oz1) #Calculate the unit vector between the first two points
						vbd = vab_unit * chosen_bondlength * (-np.cos(theta))
						xd = ox1 + vbd[0]
						yd = oy1 + vbd[1]
						zd = oz1 + vbd[2]
						v1, v2, v3 = generate_orthonormal_basis(vab_unit)  #Take ab as the z-axis to obtain a set of orthonormal basis vectors
						phi = np.random.uniform(0, 2 * np.pi, size = 1)   #Generate a polar angle on the xy-plane of the new coordinate system
						vdc_new = chosen_bondlength * np.sin(theta) * np.array([np.cos(phi), np.sin(phi), [0]]) #Obtain the coordinates of the DC vector in the new coordinate system
						vdc = transformation(v1, v2, v3, vdc_new) #Obtain the coordinates of the DC vector in the original coordinate system
						self.atom[na,4] = xd + vdc[0]
						self.atom[na,5] = yd + vdc[1]
						self.atom[na,6] = zd + vdc[2]
						ox2 = ox1    #Record the positions of the first two atoms
						oy2 = oy1
						oz2 = oz1
						ox1 = self.atom[na,4]
						oy1 = self.atom[na,5]
						oz1 = self.atom[na,6]
					na += 1
			nm += 1

	def get_bond(self):
		self.bond[:,0] = np.arange(1, self.n_bond + 1)
		na = 0
		nb = 1
		for i in range(self.n_mol):
			for j in range(self.n_seg):
				for k in range(self.segn[j]):
					na += 1
					if (k**2 + j**2 != 0): 
						self.bond[nb-1,2] = na - 1
						self.bond[nb-1,3] = na
						for l in range(int(np.max(self.bondtype[:,0]))):
							if ((self.atom[self.bond[nb-1,2]-1,2] == self.bondtype[l,1] and self.atom[self.bond[nb-1,3]-1,2] == self.bondtype[l,2]) or (self.atom[self.bond[nb-1,2]-1,2] == self.bondtype[l,2] and self.atom[self.bond[nb-1,3]-1,2] == self.bondtype[l,1])):
								self.bond[nb-1,1] = self.bondtype[l,0]
						nb += 1

	def get_angle(self):
		self.angle[:,0] = np.arange(1, self.n_angle + 1)
		na = 0
		ng = 1
		atom_per_chain = np.sum(self.segn[:])
		for i in range(self.n_mol):
			for j in range(0,self.n_seg):
				for k in range(self.segn[j]):
					na += 1
					if (na%atom_per_chain!=0 and (na-1)%atom_per_chain!=0): 
						self.angle[ng-1,2] = na - 1
						self.angle[ng-1,3] = na
						self.angle[ng-1,4] = na + 1
						for l in range(int(np.max(self.angletype[:,0]))):
							if (self.atom[self.angle[ng-1,3]-1,2] == self.angletype[l,2] and ((self.atom[self.angle[ng-1,2]-1,2] == self.angletype[l,1] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[l,3]) or (self.atom[self.angle[ng-1,2]-1,2] == self.angletype[l,3] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[l,1]))):
								self.angle[ng-1,1] = self.angletype[l,0]
						ng += 1
	
	def get_dihedral(self):
		self.dihedral[:,0] = np.arange(1, self.n_dihedral + 1)
		na = 0
		nd = 1
		atom_per_chain = np.sum(self.segn[:])
		for i in range(self.n_mol):
			for j in range(0,self.n_seg):
				for k in range(self.segn[j]):
					na += 1
					if (na%atom_per_chain!=0 and (na-1)%atom_per_chain!=0 and (na+1)%atom_per_chain!=0):
						self.dihedral[nd-1,2] = na - 1
						self.dihedral[nd-1,3] = na
						self.dihedral[nd-1,4] = na + 1
						self.dihedral[nd-1,5] = na + 2
						for l in range(int(np.max(self.dihedraltype[:,0]))):
							if((self.atom[self.dihedral[nd-1,2]-1,2] == self.dihedraltype[l,1] and 
															self.atom[self.dihedral[nd-1,3]-1,2] == self.dihedraltype[l,2] and 
															self.atom[self.dihedral[nd-1,4]-1,2] == self.dihedraltype[l,3] and 
															self.atom[self.dihedral[nd-1,5]-1,2] == self.dihedraltype[l,4]
															) or 
															(self.atom[self.dihedral[nd-1,2]-1,2] == self.dihedraltype[l,4] and 
															self.atom[self.dihedral[nd-1,3]-1,2] == self.dihedraltype[l,3] and 
															self.atom[self.dihedral[nd-1,4]-1,2] == self.dihedraltype[l,2] and 
															self.atom[self.dihedral[nd-1,5]-1,2] == self.dihedraltype[l,1]
															)):
								self.dihedral[nd-1,1] = self.dihedraltype[l,0]
						nd += 1
	
	def period_box(self):
		for i in range(self.n_atom):
			period_x = math.floor(self.atom[i,4]/self.lx)
			self.atom[i,4] -= self.lx * period_x
			period_y = math.floor(self.atom[i,5]/self.ly)
			self.atom[i,5] -= self.ly * period_y
			period_z = math.floor(self.atom[i,6]/self.lz)
			self.atom[i,6] -= self.lz * period_z

	def write_data(self):
		self.file.write('LAMMPS Description\n\n')
		self.file.write('%d atoms\n%d bonds\n%d angles\n%d dihedrals\n0 impropers\n'%(self.n_atom,self.n_bond,self.n_angle,self.n_dihedral))
		self.file.write('%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n0 improper types\n'%(self.n_atomtype,np.max(self.bondtype[:,0]),np.max(self.angletype[:,0]),np.max(self.dihedraltype[:,0])))
		self.file.write('0 %f xlo xhi\n0 %f ylo yhi\n0 %f zlo zhi\n\n'%(self.lx,self.ly,self.lz))
		self.file.write('Masses\n\n')
		for i in range(1,self.n_atomtype+1):
			self.file.write('%d %f\n'%(i,self.mass[i-1]))
		self.file.write('\nAtoms  # full\n\n')
		for i in range(self.n_atom):
			self.file.write('%d %d %d %f %f %f %f\n'%(self.atom[i,0],self.atom[i,1],self.atom[i,2],self.atom[i,3],self.atom[i,4],self.atom[i,5],self.atom[i,6]))
		if self.n_bond != 0:
			self.file.write('\nBonds\n\n')
			for i in range(self.n_bond):
				self.file.write('%d %d %d %d\n'%(self.bond[i,0],self.bond[i,1],self.bond[i,2],self.bond[i,3]))
		if self.n_angle != 0:
			self.file.write('\nAngles\n\n')
			for i in range(self.n_angle):
				self.file.write('%d %d %d %d %d\n'%(self.angle[i,0],self.angle[i,1],self.angle[i,2],self.angle[i,3],self.angle[i,4]))
		if self.n_dihedral != 0:
			self.file.write('\nDihedrals\n\n')
			for i in range(self.n_dihedral):
				self.file.write('%d %d %d %d %d %d\n'%(self.dihedral[i,0],self.dihedral[i,1],self.dihedral[i,2],self.dihedral[i,3],self.dihedral[i,4],self.dihedral[i,5]))
		self.file.close()

	def write_pdb(self):
		with open('system.pdb', 'w') as pdb_file:
			pdb_file.write("REMARK   LAMMPS data file converted to PDB format\n")
			for i in range(self.n_atom):
				if self.atom[i,2] == 1:
					atomtype = 'C'
				elif self.atom[i,2] == 2:
					atomtype = 'O'
				atom_name = atomtype.rjust(4) # Atom name, right-aligned to 4 characters
				element = atomtype           # Element symbol, usually two characters in PDB, no alignment needed
				residue_name = 'UNL'         # Residue name, exactly 3 characters, left-aligned (already default)
				chain_identifier = 'A'       # Chain identifier, one character
				pdb_file.write('HETATM%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n' % (
					self.atom[i, 0],  # Atom serial number, right-aligned to 5 characters
					atom_name,        # Atom name, right-aligned to 4 characters (but usually no alignment needed)
					residue_name,     # Residue name, left-aligned to 3 characters
					chain_identifier, # Chain identifier, 1 character wide
					self.atom[i, 1],  # Residue sequence number, right-aligned to 4 characters
					self.atom[i, 4],  # X coordinate, 8 characters wide, 3 digits after decimal point
					self.atom[i, 5],  # Y coordinate, same as above
					self.atom[i, 6],  # Z coordinate, same as above
					1.00,             # Occupancy, 6 characters wide, 2 digits after decimal point (fixed to 1.00 here)
					0.00,             # Temperature factor, 6 characters wide, 2 digits after decimal point (fixed to 0.00 here)
					element           # Element symbol, left-aligned to 2 characters (but already default)
				))
			pdb_file.write('END\n')

class Sequential_chain(Chain):	
	def	get_atom(self):
		na = 0
		nm = 1
		self.atom[:,0] = np.arange(1, self.n_atom + 1)
		for i in range(self.n_mol):
			for j in range(self.n_seg):
				for k in range(self.segn[j]):
					self.atom[na,2] = self.segt[k]
					self.atom[na,1] = nm
					self.atom[na,3] = 0.0   
					if k == 0 and j == 0:   
						self.atom[na,4] = np.random.uniform(0, self.lx, size = 1)
						self.atom[na,5] = np.random.uniform(0, self.ly, size = 1)
						self.atom[na,6] = np.random.uniform(0, self.lz, size = 1)
						ox1 = self.atom[na,4]
						oy1 = self.atom[na,5]
						oz1 = self.atom[na,6]
					elif (k == 1 and j == 0) or (self.segn[0] == 1 and k == 0 and j == 1):  
						flag = True
						for l in range(int(np.max(self.bondtype[:,0]))):
							if ((self.atom[na,2] == self.bondtype[l,1] and self.atom[na-1,2] == self.bondtype[l,2]) or (self.atom[na,2] == self.bondtype[l,2] and self.atom[na-1,2] == self.bondtype[l,1])):
								peak = int(self.bondcoeff[self.bondtype[l,0]-1][0])
								bondlength = np.zeros(peak)
								amplitude = np.zeros(peak)
								for i in range(peak):
									mean = float(self.bondcoeff[self.bondtype[l,0]-1][3*i+3])
									std_dev = float(self.bondcoeff[self.bondtype[l,0]-1][3*i+2])
									bondlength[i] = np.random.normal(mean, std_dev, size = 1)
									amplitude[i] = float(self.bondcoeff[self.bondtype[l,0]-1][3*i+1])
								amplitude = amplitude/np.sum(amplitude)
								chosen_index = np.random.choice(len(bondlength), p=amplitude)
								chosen_bondlength = bondlength[chosen_index] 
						while flag:
							theta = np.random.uniform(0, 2*np.pi, size = 1)
							phi = np.random.uniform(0, np.pi, size = 1)
							self.atom[na,4] = ox1 + chosen_bondlength * np.sin(phi) * np.cos(theta)   
							self.atom[na,5] = oy1 + chosen_bondlength * np.sin(phi) * np.sin(theta)
							self.atom[na,6] = oz1 + chosen_bondlength * np.cos(phi)
							if (self.atom[na,4]>=0 and self.atom[na,4]<self.lx and self.atom[na,5]>=0 and self.atom[na,5]<self.ly and self.atom[na,6]>=0 and self.atom[na,6]<self.lz):
								ox2 = ox1    
								oy2 = oy1
								oz2 = oz1
								ox1 = self.atom[na,4]
								oy1 = self.atom[na,5]
								oz1 = self.atom[na,6]
								flag = False
					else:
						flag = True
						for l in range(int(np.max(self.bondtype[:,0]))):
							if ((self.atom[na,2] == self.bondtype[l,1] and self.atom[na-1,2] == self.bondtype[l,2]) or (self.atom[na,2] == self.bondtype[l,2] and self.atom[na-1,2] == self.bondtype[l,1])):
								peak = int(self.bondcoeff[self.bondtype[l,0]-1][0])
								bondlength = np.zeros(peak)
								amplitude = np.zeros(peak)
								for m in range(peak):
									mean = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+3])
									std_dev = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+2])
									bondlength[m] = np.random.normal(mean, std_dev, size = 1)
									amplitude[m] = float(self.bondcoeff[self.bondtype[l,0]-1][3*m+1])
								amplitude = amplitude/np.sum(amplitude)
								chosen_index = np.random.choice(len(bondlength), p=amplitude)
								chosen_bondlength = bondlength[chosen_index]
						for l in range(int(np.max(self.angletype[:,0]))):
							if (self.atom[na-1,2] == self.angletype[l,2] and ((self.atom[na-2,2] == self.angletype[l,1] and self.atom[na,2] == self.angletype[l,3]) or (self.atom[na-2,2] == self.angletype[l,3] and self.atom[na,2] == self.angletype[l,1]))):
								peak = int(self.anglecoeff[self.angletype[l,0]-1][0])
								angledegree = np.zeros(peak)
								amplitude = np.zeros(peak)
								for m in range(peak):
									mean = float(self.anglecoeff[self.angletype[l,0]-1][3*m+3])
									std_dev = float(self.anglecoeff[self.angletype[l,0]-1][3*m+2])
									angledegree[m] =  np.random.normal(mean, std_dev, size = 1)
									amplitude[m] = float(self.anglecoeff[self.angletype[l,0]-1][3*m+1])
								amplitude = amplitude/np.sum(amplitude)
								chosen_index = np.random.choice(len(angledegree), p=amplitude)
								theta = np.radians(angledegree[chosen_index])
						vab_unit = unit_vector(ox2,oy2,oz2,ox1,oy1,oz1) 
						vbd = vab_unit * chosen_bondlength * (-np.cos(theta))
						xd = ox1 + vbd[0]
						yd = oy1 + vbd[1]
						zd = oz1 + vbd[2]
						step = 0
						while flag:
							step += 1
							if step < 50:
								v1, v2, v3 = generate_orthonormal_basis(vab_unit)  
								phi = np.random.uniform(0, 2 * np.pi, size = 1)  
								vdc_new = chosen_bondlength * np.sin(theta) * np.array([np.cos(phi), np.sin(phi), [0]]) 
								vdc = transformation(v1, v2, v3, vdc_new)
								self.atom[na,4] = xd + vdc[0]
								self.atom[na,5] = yd + vdc[1]
								self.atom[na,6] = zd + vdc[2]
								if (self.atom[na,4]>=0 and self.atom[na,4]<self.lx and self.atom[na,5]>=0 and self.atom[na,5]<self.ly and self.atom[na,6]>=0 and self.atom[na,6]<self.lz):
									ox2 = ox1    
									oy2 = oy1
									oz2 = oz1
									ox1 = self.atom[na,4]
									oy1 = self.atom[na,5]
									oz1 = self.atom[na,6]
									flag = False
							else:
								theta = np.random.uniform(0, 2*np.pi, size = 1)
								phi = np.random.uniform(0, np.pi, size = 1)
								self.atom[na,4] = ox1 + chosen_bondlength * np.sin(phi) * np.cos(theta)
								self.atom[na,5] = oy1 + chosen_bondlength * np.sin(phi) * np.sin(theta)
								self.atom[na,6] = oz1 + chosen_bondlength * np.cos(phi)
								if (self.atom[na,4]>=0 and self.atom[na,4]<self.lx and self.atom[na,5]>=0 and self.atom[na,5]<self.ly and self.atom[na,6]>=0 and self.atom[na,6]<self.lz):
									ox2 = ox1    
									oy2 = oy1
									oz2 = oz1
									ox1 = self.atom[na,4]
									oy1 = self.atom[na,5]
									oz1 = self.atom[na,6]
									flag = False
						
					na += 1
			nm += 1

class Branched_chain(Chain):
	def __init__(self, n_mol, n_seg, segn, segt, lx, ly, lz, mass, bondtype, angletype, dihedraltype, bondcoeff, anglecoeff, n_branch, branchtype):
		Chain.__init__(self, n_mol, n_seg, segn, segt, lx, ly, lz, mass, bondtype, angletype, dihedraltype, bondcoeff, anglecoeff)
		self.n_branch = np.array(n_branch)   #n_branch is a 2D array, [[branch length for each atom in block 1], [branch length for each atom in block 2], ...]
		self.branchtype = np.array(branchtype)  #branchtype is a 2D array, [[branch atom type for each atom in block 1], [branch atom type for each atom in block 2], ...]
		self.n_atom_back = self.n_atom
		self.n_bond_back = self.n_bond
		self.n_angle_back = self.n_angle
		self.n_dihedral_back = self.n_dihedral
		self.n_atom += np.sum(n_branch) * n_mol
		self.n_atomtype = max(np.max(self.segt), np.max(branchtype))
		self.n_bond += np.sum(n_branch) * n_mol
		A = np.subtract(n_branch, 1)
		A [A<0] = 0
		B = 0
		for j in range(self.n_seg):
			for k in range(self.segn[j]):
				if self.n_branch[j,k] != 0:
					if k == 0 and j == 0:
						B += 1
					elif (k == self.segn[j] - 1 and j == self.n_seg - 1):
						B += 1
					else:
						B += 2
		self.n_angle += (np.sum(A) + B) * n_mol
		self.n_dihedral += (np.sum(self.segn) - 1) * n_mol
		self.atom = np.zeros((self.n_atom, 7))
		self.bond = np.zeros((self.n_bond, 4), dtype = int)
		self.angle = np.zeros((self.n_angle, 5), dtype = int)
		self.dihedral = np.zeros((self.n_dihedral,6), dtype = int)
		print(self.n_atom)
		print(self.n_bond)
		print(self.n_angle)

	def get_atom(self):
		Chain.get_atom(self)
		na = self.n_atom_back
		nm = 1
		atom_per_backchain = np.sum(self.segn)
		for i in range(self.n_mol):
			for j in range(self.n_seg):
				for k in range(self.segn[j]):
					for l in range(self.n_branch[j,k]):
						self.atom[na,2] = self.branchtype[j,k]
						self.atom[na,1] = nm
						self.atom[na,3] = 0.0   #电荷
						if l == 0:
							if k == 0 and j == 0:
								preatom = i * np.sum(self.segn) + np.sum(self.segn[0:j]) + k
								ox1 = self.atom[preatom, 4]
								oy1 = self.atom[preatom, 5]
								oz1 = self.atom[preatom, 6]
								ox2 = self.atom[preatom + 1, 4]
								oy2 = self.atom[preatom + 1, 5]
								oz2 = self.atom[preatom + 1, 6]
								for n in range(int(np.max(self.bondtype[:,0]))):
									if ((self.atom[na,2] == self.bondtype[n,1] and self.atom[preatom,2] == self.bondtype[n,2]) or (self.atom[na,2] == self.bondtype[n,2] and self.atom[preatom,2] == self.bondtype[n,1])):	
										peak = int(self.bondcoeff[self.bondtype[n,0]-1][0])
										bondlength = np.zeros(peak)
										amplitude = np.zeros(peak)
										for m in range(peak):
											mean = float(self.bondcoeff[self.bondtype[n,0]-1][3*m+3])
											std_dev = float(self.bondcoeff[self.bondtype[n,0]-1][3*m+2])
											bondlength[m] = np.random.normal(mean, std_dev, size = 1)
											amplitude[m] = float(self.bondcoeff[self.bondtype[n,0]-1][3*m+1])
										amplitude = amplitude/np.sum(amplitude)
										chosen_index = np.random.choice(len(bondlength), p=amplitude)
										chosen_bondlength = bondlength[chosen_index]
								for n in range(int(np.max(self.angletype[:,0]))):
									if (self.atom[preatom,2] == self.angletype[n,2] and ((self.atom[preatom + 1,2] == self.angletype[n,1] and self.atom[na,2] == self.angletype[n,3]) or (self.atom[preatom + 1,2] == self.angletype[n,3] and self.atom[na,2] == self.angletype[n,1]))):
										peak = int(self.anglecoeff[self.angletype[n,0]-1][0])
										angledegree = np.zeros(peak)
										amplitude = np.zeros(peak)
										for m in range(peak):
											mean = float(self.anglecoeff[self.angletype[n,0]-1][3*m+3])
											std_dev = float(self.anglecoeff[self.angletype[n,0]-1][3*m+2])
											angledegree[m] =  np.random.normal(mean, std_dev, size = 1)
											amplitude[m] = float(self.anglecoeff[self.angletype[n,0]-1][3*m+1])
										amplitude = amplitude/np.sum(amplitude)
										chosen_index = np.random.choice(len(angledegree), p=amplitude)
										theta = np.radians(angledegree[chosen_index])
							else:
								preatom = i * np.sum(self.segn) + np.sum(self.segn[0:j]) + k
								ox1 = self.atom[preatom, 4]
								oy1 = self.atom[preatom, 5]
								oz1 = self.atom[preatom, 6]
								ox2 = self.atom[preatom - 1, 4]
								oy2 = self.atom[preatom - 1, 5]
								oz2 = self.atom[preatom - 1, 6]
								for n in range(int(np.max(self.bondtype[:,0]))):
									if ((self.atom[na,2] == self.bondtype[n,1] and self.atom[preatom,2] == self.bondtype[n,2]) or (self.atom[na,2] == self.bondtype[n,2] and self.atom[preatom,2] == self.bondtype[n,1])):	
										peak = int(self.bondcoeff[self.bondtype[n,0]-1][0])
										bondlength = np.zeros(peak)
										amplitude = np.zeros(peak)
										for m in range(peak):
											mean = float(self.bondcoeff[self.bondtype[n,0]-1][3*m+3])
											std_dev = float(self.bondcoeff[self.bondtype[n,0]-1][3*m+2])
											bondlength[m] = np.random.normal(mean, std_dev, size = 1)
											amplitude[m] = float(self.bondcoeff[self.bondtype[n,0]-1][3*m+1])
										amplitude = amplitude/np.sum(amplitude)
										chosen_index = np.random.choice(len(bondlength), p=amplitude)
										chosen_bondlength = bondlength[chosen_index]
								for n in range(int(np.max(self.angletype[:,0]))):
									if (self.atom[preatom,2] == self.angletype[n,2] and ((self.atom[preatom - 1,2] == self.angletype[n,1] and self.atom[na,2] == self.angletype[n,3]) or (self.atom[preatom - 1,2] == self.angletype[n,3] and self.atom[na,2] == self.angletype[n,1]))):
										peak = int(self.anglecoeff[self.angletype[n,0]-1][0])
										angledegree = np.zeros(peak)
										amplitude = np.zeros(peak)
										for m in range(peak):
											mean = float(self.anglecoeff[self.angletype[n,0]-1][3*m+3])
											std_dev = float(self.anglecoeff[self.angletype[n,0]-1][3*m+2])
											angledegree[m] =  np.random.normal(mean, std_dev, size = 1)
											amplitude[m] = float(self.anglecoeff[self.angletype[n,0]-1][3*m+1])
										amplitude = amplitude/np.sum(amplitude)
										chosen_index = np.random.choice(len(angledegree), p=amplitude)
										theta = np.radians(angledegree[chosen_index]) 
						flag = True
						vab_unit = unit_vector(ox2,oy2,oz2,ox1,oy1,oz1)
						vbd = vab_unit * chosen_bondlength * (-np.cos(theta))
						xd = ox1 + vbd[0]
						yd = oy1 + vbd[1]
						zd = oz1 + vbd[2]
						step = 0
						while flag:
							step += 1
							if step < 50:
								v1, v2, v3 = generate_orthonormal_basis(vab_unit)
								phi = np.random.uniform(0, 2 * np.pi, size = 1)
								vdc_new = chosen_bondlength * np.sin(theta) * np.array([np.cos(phi), np.sin(phi), [0]])
								vdc = transformation(v1, v2, v3, vdc_new)
								self.atom[na,4] = xd + vdc[0]
								self.atom[na,5] = yd + vdc[1]
								self.atom[na,6] = zd + vdc[2]
								if (self.atom[na,4]>=0 and self.atom[na,4]<self.lx and self.atom[na,5]>=0 and self.atom[na,5]<self.ly and self.atom[na,6]>=0 and self.atom[na,6]<self.lz):
									ox2 = ox1
									oy2 = oy1
									oz2 = oz1
									ox1 = self.atom[na,4]
									oy1 = self.atom[na,5]
									oz1 = self.atom[na,6]
									flag = False
							else:
								theta = np.random.uniform(0, 2*np.pi, size = 1)
								phi = np.random.uniform(0, np.pi, size = 1)
								self.atom[na,4] = ox1 + chosen_bondlength * np.sin(phi) * np.cos(theta)
								self.atom[na,5] = oy1 + chosen_bondlength * np.sin(phi) * np.sin(theta)
								self.atom[na,6] = oz1 + chosen_bondlength * np.cos(phi)
								if (self.atom[na,4]>=0 and self.atom[na,4]<self.lx and self.atom[na,5]>=0 and self.atom[na,5]<self.ly and self.atom[na,6]>=0 and self.atom[na,6]<self.lz):
									ox2 = ox1
									oy2 = oy1
									oz2 = oz1
									ox1 = self.atom[na,4]
									oy1 = self.atom[na,5]
									oz1 = self.atom[na,6]
									flag = False		
						na += 1
			nm += 1

	def get_bond(self):
		Chain.get_bond(self)
		na = self.n_atom_back + 1
		nb = self.n_bond_back + 1
		for i in range(self.n_mol):
			for j in range(self.n_seg):
				for k in range(self.segn[j]):
					for l in range(self.n_branch[j,k]):
						if l == 0:
							preatom = i * np.sum(self.segn) + np.sum(self.segn[0:j]) + k + 1
							self.bond[nb-1,2] = preatom
							self.bond[nb-1,3] = na
						else:
							self.bond[nb-1,2] = na - 1
							self.bond[nb-1,3] = na
						for m in range(int(np.max(self.bondtype[:,0]))):
							if ((self.atom[self.bond[nb-1,2]-1,2] == self.bondtype[m,1] and self.atom[self.bond[nb-1,3]-1,2] == self.bondtype[m,2]) or (self.atom[self.bond[nb-1,2]-1,2] == self.bondtype[m,2] and self.atom[self.bond[nb-1,3]-1,2] == self.bondtype[m,1])):
								self.bond[nb-1,1] = self.bondtype[m,0]
						na += 1
						nb += 1

	def get_angle(self):
		Chain.get_angle(self)
		na = self.n_atom_back + 1
		ng = self.n_angle_back + 1
		atom_per_chain = np.sum(self.segn)
		for i in range(self.n_mol):
			for j in range(0,self.n_seg):
				for k in range(self.segn[j]):
					for l in range(self.n_branch[j,k]):
						if l == 0:
							if k == 0 and j == 0:
								preatom = i * atom_per_chain + np.sum(self.segn[0:j]) + k + 1
								self.angle[ng-1,2] = preatom + 1
								self.angle[ng-1,3] = preatom
								self.angle[ng-1,4] = na
								for m in range(int(np.max(self.angletype[:,0]))):
									if (self.atom[self.angle[ng-1,3]-1,2] == self.angletype[m,2] and ((self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,1] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,3]) or (self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,3] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,1]))):
										self.angle[ng-1,1] = self.angletype[m,0]
								na += 1
								ng += 1
							elif (k == self.segn[j] - 1 and j == self.n_seg - 1):
								preatom = i * atom_per_chain + np.sum(self.segn[0:j]) + k + 1
								self.angle[ng-1,2] = preatom - 1
								self.angle[ng-1,3] = preatom
								self.angle[ng-1,4] = na
								for m in range(int(np.max(self.angletype[:,0]))):
									if (self.atom[self.angle[ng-1,3]-1,2] == self.angletype[m,2] and ((self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,1] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,3]) or (self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,3] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,1]))):
										self.angle[ng-1,1] = self.angletype[m,0]
								na += 1
								ng += 1
							else:
								preatom = i * atom_per_chain + np.sum(self.segn[0:j]) + k + 1
								self.angle[ng-1,2] = preatom - 1
								self.angle[ng-1,3] = preatom
								self.angle[ng-1,4] = na
								for m in range(int(np.max(self.angletype[:,0]))):
									if (self.atom[self.angle[ng-1,3]-1,2] == self.angletype[m,2] and ((self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,1] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,3]) or (self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,3] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,1]))):
										self.angle[ng-1,1] = self.angletype[m,0]
								self.angle[ng,2] = preatom + 1
								self.angle[ng,3] = preatom
								self.angle[ng,4] = na
								for m in range(int(np.max(self.angletype[:,0]))):
									if (self.atom[self.angle[ng,3]-1,2] == self.angletype[m,2] and ((self.atom[self.angle[ng,2]-1,2] == self.angletype[m,1] and self.atom[self.angle[ng,4]-1,2] == self.angletype[m,3]) or (self.atom[self.angle[ng,2]-1,2] == self.angletype[m,3] and self.atom[self.angle[ng,4]-1,2] == self.angletype[m,1]))):
										self.angle[ng,1] = self.angletype[m,0]
								na += 1
								ng += 2
						elif l == 1:
							preatom = i * atom_per_chain + np.sum(self.segn[0:j]) + k + 1
							self.angle[ng-1,2] = na
							self.angle[ng-1,3] = na - 1
							self.angle[ng-1,4] = preatom
							for m in range(int(np.max(self.angletype[:,0]))):
								if (self.atom[self.angle[ng-1,3]-1,2] == self.angletype[m,2] and ((self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,1] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,3]) or (self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,3] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,1]))):
									self.angle[ng-1,1] = self.angletype[m,0]
							na += 1
							ng += 1
						else:
							self.angle[ng-1,2] = na
							self.angle[ng-1,3] = na - 1
							self.angle[ng-1,4] = na - 2
							for m in range(int(np.max(self.angletype[:,0]))):
								if (self.atom[self.angle[ng-1,3]-1,2] == self.angletype[m,2] and ((self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,1] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,3]) or (self.atom[self.angle[ng-1,2]-1,2] == self.angletype[m,3] and self.atom[self.angle[ng-1,4]-1,2] == self.angletype[m,1]))):
									self.angle[ng-1,1] = self.angletype[m,0]
							na += 1
							ng += 1
	def get_dihedral(self):
		Chain.get_dihedral(self)
		na = self.n_atom_back + 1
		nd = self.n_dihedral_back + 1
		atom_per_chain = np.sum(self.segn)
		for i in range(self.n_mol):
			for j in range(0,self.n_seg):
				for k in range(self.segn[j]):
					preatom = i * atom_per_chain + np.sum(self.segn[0:j]) + k + 1
					if preatom % atom_per_chain != 0:
						self.dihedral[nd-1,2] = na
						self.dihedral[nd-1,3] = preatom
						self.dihedral[nd-1,4] = preatom + 1
						self.dihedral[nd-1,5] = na + self.n_branch[j,k]
						for l in range(int(np.max(self.dihedraltype[:,0]))):
							if((self.atom[self.dihedral[nd-1,2]-1,2] == self.dihedraltype[l,1] and 
															self.atom[self.dihedral[nd-1,3]-1,2] == self.dihedraltype[l,2] and 
															self.atom[self.dihedral[nd-1,4]-1,2] == self.dihedraltype[l,3] and 
															self.atom[self.dihedral[nd-1,5]-1,2] == self.dihedraltype[l,4]
															) or 
															(self.atom[self.dihedral[nd-1,2]-1,2] == self.dihedraltype[l,4] and 
															self.atom[self.dihedral[nd-1,3]-1,2] == self.dihedraltype[l,3] and 
															self.atom[self.dihedral[nd-1,4]-1,2] == self.dihedraltype[l,2] and 
															self.atom[self.dihedral[nd-1,5]-1,2] == self.dihedraltype[l,1]
															)):
								self.dihedral[nd-1,1] = self.dihedraltype[l,0]
						na += 1
						nd += 1	

def unit_vector(x1,y1,z1,x2,y2,z2): #Calculate the unit vector between two points
	v = np.array([x2 - x1, y2 - y1, z2 - z1])
	norm = np.linalg.norm(v)
	if norm != 0: #If it is not a zero vector
		v = v / norm
	else:
		v = np.zeros_like(v)
	return v

def generate_orthonormal_basis(v3):                                            #Generate a set of orthonormal basis vectors from v3 via Gram-Schmidt orthonormalization
	assert np.isclose(np.linalg.norm(v3), 1.0), "Input must be a unit vector"  #Ensure the input is a unit vector
	v2_raw = np.random.rand(3)
	v2_raw = v2_raw - np.dot(v2_raw, v3) * v3                                  
	while np.linalg.norm(v2_raw)<0.001:                                        #Generate a vector v2_raw that is not collinear with v3
		v2_raw = np.random.rand(3)
		v2_raw = v2_raw - np.dot(v2_raw, v3) * v3                              #Gram-Schmidt orthonormalization
	v2 = v2_raw / np.linalg.norm(v2_raw)                                       #Normalize v2_raw to obtain v2
	v1 = np.cross(v2, v3)                                                      #Use cross product to generate the third vector v1
	return v1, v2, v3

def transformation(v1, v2, v3, v): #Convert a vector in the new coordinate system to a vector in the Cartesian coordinate system
	transformation_matrix = np.column_stack((v1, v2, v3)) #The transformation matrix is the coordinates of the basis vectors of the new coordinate system in the Cartesian coordinate system
	v_cartesian = np.dot(transformation_matrix, v) #Multiply the transformation matrix by the column vector to obtain the column vector in the Cartesian coordinate system
	return v_cartesian

