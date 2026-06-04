from openbabel import pybel
import openbabel as ob
import numpy as np
import MDAnalysis as mda

class MoleculeBuilder:
	#Obtain atomic information of the molecule: smi is the SMILES of the molecule, index is the molecule serial number

	def __init__(self, smi, n_moltype):
		self.mol = pybel.readstring("smi", smi)
		self.n_moltype = n_moltype

	def optchain(self):
		#Obtain the .pdb file of a single molecule
		self.mol.addh()
		self.mol.make3D(forcefield='mmff94', steps=5000)   
		self.mol.write("pdb", "Chain_%d.pdb"%self.n_moltype, overwrite=True)# Write the optimized molecule to a pdb file
		return 0

	def get_moltemplate_monomer_lt(self):
		#Obtain atom types and charges of the molecule, match with Dreiding force field, and generate the .lt file for Moltemplate
		self.mol.addh()  #Add H atoms
		self.mol.make3D(forcefield="uff")  #Generate 3D coordinates and optimize using the UFF force field
		N = self.mol.OBMol.NumAtoms()  #Get the number of atoms
		atoms = np.zeros((N, 6)).tolist()
		partialCharges = ob.OBChargeModel.FindType("Gasteiger")  # Use the Gasteiger charge model  
		partialCharges.ComputeCharges(self.mol.OBMol) 
		for idx, atom in enumerate(self.mol.atoms):
			atoms[idx][0] = atom.type  # Obtain atom types
			print(atoms[idx][0])
			if atoms[idx][0] == "C3":
				atoms[idx][1] = 'C_3'
			elif atoms[idx][0] == "C2":
				atoms[idx][1] = 'C_2'
			elif atoms[idx][0] == "C1":
				atoms[idx][1] = 'C_1'
			elif atoms[idx][0] == "Nar":
				atoms[idx][1] = 'N_R_d1'
			elif atoms[idx][0] == "Car":
				atoms[idx][1] = 'C_R'
			elif atoms[idx][0] == "O3":
				atoms[idx][1] = 'O_3'
			elif atoms[idx][0] == "O2":
				atoms[idx][1] = 'O_2'
			elif atoms[idx][0] == "N3":
				atoms[idx][1] = 'N_3'
			elif atoms[idx][0] == "Npl":
				atoms[idx][1] = 'N_3'
			elif atoms[idx][0] == "Nam":
				atoms[idx][1] = 'N_2_d1_hd'
			elif atoms[idx][0] == "N2":
				atoms[idx][1] = 'N_2_b2_d1'
			elif atoms[idx][0] == "N1":
				atoms[idx][1] = 'N_1'
			elif atoms[idx][0] == "H":
				atoms[idx][1] = 'H'
			elif atoms[idx][0] == "F":
				atoms[idx][1] = 'F'
			elif atoms[idx][0] == "HO":
				atoms[idx][1] = 'H_HB'
			elif atoms[idx][0] == "Sac":
				atoms[idx][1] = 'S_3'
			elif atoms[idx][0] == "S3":
				atoms[idx][1] = 'S_3'
			elif atoms[idx][0] == "So2":
				atoms[idx][1] = 'S_3'
			elif atoms[idx][0] == "S2":
				atoms[idx][1] = 'S_2'
			elif atoms[idx][0] == "P":
				atoms[idx][1] = 'P'	
			pos = atom.coords  # Compute atomic coordinate
			atoms[idx][3] = float(pos[0])
			atoms[idx][4] = float(pos[1])
			atoms[idx][5] = float(pos[2])
			atoms[idx][2] = atom.OBAtom.GetPartialCharge() 
		moltemplate_monomer = open("Polymer_%d.lt"%self.n_moltype,"w")
		moltemplate_monomer.write('''import "dreiding.lt" \n\nPolymer%d inherits DREIDING \n{\n\nwrite("Data Atoms") \n{\n'''%self.n_moltype)
		for i in range(N):
			moltemplate_monomer.write('    $atom:%d  $mol:0  @atom:%s\t%f\t%f\t%f\t%f\n'%(i,atoms[i][1],atoms[i][2],atoms[i][3] ,atoms[i][4] ,atoms[i][5]))
		moltemplate_monomer.write('''  }\nwrite("Data Bond List")\n   {\n''')
		
		#Obtain bonded atom indices via MDAnalysis
		u = mda.Universe('Chain_%d.pdb'%self.n_moltype)
		i = 0
		for bond in u.bonds:
			atom1 = bond[0].index
			atom2 = bond[1].index
			moltemplate_monomer.write('$bond:%d\t$atom:%d\t$atom:%d\n'%(i,atom1,atom2))
			i+=1
		moltemplate_monomer.write('''  }\n\n\n}\n''')
		moltemplate_monomer.close()
		return 0