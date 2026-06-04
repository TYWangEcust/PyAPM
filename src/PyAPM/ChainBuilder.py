from rdkit.Chem import AllChem
import numpy as np


class ChainBuilder:
	def __init__(self, smi, length, component):
		self.smi = smi     #Input smi: SMILES of the repeating unit. If component == 1, smi is a string; if component > 1, smi is an array
		self.n = length         #n: degree of polymerization, number of repetitions of the repeating unit
		self.component = component

	def get_homopolymer_structure(self):
		self.smi = self.smi.replace("[*]","X").replace("*","X") #Unify the end symbol to X
		self.smi = self.smi.replace("X","[Ar]",1).replace("X","[He]",1) #Then represent the end symbols as Ar and He
		try:
			mol=AllChem.MolFromSmiles(self.smi)
		except:
			print("Could not load molecule\n")
			return None
		
		#Prepare two types of reactions
		reactants=(mol,mol)
		rxn_str1 ='[Ar:1][c,C,n,N,O,s,S,Si,P:2].[He:3][c,C,n,N,O,s,S,Si,P:4]>>[c,C,n,N,O,s,S,Si,P:2]-[c,C,n,N,O,s,S,Si,P:4].[He:3][Ar:1]'
		rxn_str2 ='[Ar:1]=[c,C,n,N,O,s,S,Si,P:2].[He:3]=[c,C,n,N,O,s,S,Si,P:4]>>[c,C,n,N,O,s,S,Si,P:2]=[c,C,n,N,O,s,S,Si,P:4].[He:3]=[Ar:1]'
		exception_code=0
		#Test dimer, confirm reaction rxn
		try:
			rxn=AllChem.ReactionFromSmarts(rxn_str1)
			products=rxn.RunReactants(reactants)
			mol_test=products[0][0] #Dimer mol object
			rxn_str=rxn_str1
		except:
			print("Failed to load reaction from rxn_str1\n")
			try:
				rxn=AllChem.ReactionFromSmarts(rxn_str2)
				products=rxn.RunReactants(reactants)
				mol_test=products[0][0]
				rxn_str=rxn_str2
				exception_code=1
			except:
				print("Failed to load reaction from rxn_str2\n")
				print("Could not load reaction\n")
				return None
		rxn=AllChem.ReactionFromSmarts(rxn_str)

		#Perform reaction
		mol_temp=mol
		for i in range(0,self.n-1):
			reactants_temp=(mol,mol_temp)
			products=rxn.RunReactants(reactants_temp)
			mol_temp=products[0][0]

		self.smi=AllChem.MolToSmiles(mol_temp)    #Obtain SMILES of the n-mer
		if exception_code:
			end_logo="[C]"
		else:
			end_logo="[H]"
		# Convert the product SMILES back to a string ending with *
		self.smi = self.smi.replace("[He]",end_logo).replace("[Ar]",end_logo)
		#Canonicalize SMILES
		chain = AllChem.MolFromSmiles(self.smi)
		chain = AllChem.RemoveHs(chain)
		self.smi = AllChem.MolToSmiles(chain)
		print(self.smi)
		return self.smi

	def get_copolymer_structure(self, co_array):  
		# Input self.smi array, containing self.component components
		# self.n: degree of polymerization, number of repetitions of the repeating unit
		# co_array: vector of repeating unit types for the copolymer
		for i in range(self.component):
			self.smi[i] = self.smi[i].replace("[*]","X").replace("*","X")
			self.smi[i] = self.smi[i].replace("X","[Ar]",1).replace("X","[He]",1)

		# Initialize reactants and reactions
		reactants=[AllChem.MolFromSmiles(self.smi[co_array[0]-1]), AllChem.MolFromSmiles(self.smi[co_array[1]-1])] 
		rxn_str1 ='[Ar:1][c,C,n,N,O,s,S,Si,P:2].[He:3][c,C,n,N,O,s,S,Si,P:4]>>[c,C,n,N,O,s,S,Si,P:2]-[c,C,n,N,O,s,S,Si,P:4].[He:3][Ar:1]'
		rxn_str2 ='[Ar:1]=[c,C,n,N,O,s,S,Si,P:2].[He:3]=[c,C,n,N,O,s,S,Si,P:4]>>[c,C,n,N,O,s,S,Si,P:2]=[c,C,n,N,O,s,S,Si,P:4].[He:3]=[Ar:1]'
		#Test dimer, confirm reaction rxn
		exception_code=0
		try:
			rxn = AllChem.ReactionFromSmarts(rxn_str1)
			products=rxn.RunReactants(reactants)
			mol_test=products[0][0] #Dimer mol object
			rxn_str = rxn_str1
		except:
			print("Failed to load reaction from rxn_str1\n")
			try:
				rxn = AllChem.ReactionFromSmarts(rxn_str2)
				products=rxn.RunReactants(reactants)
				mol_test=products[0][0]
				rxn_str = rxn_str2
				exception_code = 1
			except:
				print("Failed to load reaction from rxn_str2\n")
				print("Could not load reaction")
				return None  

		rxn=AllChem.ReactionFromSmarts(rxn_str)#Select reaction
		
		# Start building the copolymer
		mol_temp = reactants[0]
		for i in range(1, self.n):
			next_mol = AllChem.MolFromSmiles(self.smi[co_array[i]-1])
			reactants_temp = (mol_temp, next_mol)
			products = rxn.RunReactants(reactants_temp)
			if products:
				mol_temp = products[0][0]
			else:
				print("Reaction failed at iteration", i)
				return None
				
		# Process the final molecule
		self.smi = AllChem.MolToSmiles(mol_temp)
		if exception_code:
			end_logo = "[C]"
		else:
			end_logo = "[H]"
		self.smi = self.smi.replace("[He]", end_logo).replace("[Ar]", end_logo)
		chain = AllChem.MolFromSmiles(self.smi)
		chain = AllChem.RemoveHs(chain)
		self.smi = AllChem.MolToSmiles(chain)
		print(self.smi)
		return self.smi

	def get_random_array(self, q):
		# Generate random copolymer
		# q is the proportion of each component, length = self.component, sum should be 1
		q = np.array(q)
		if not np.isclose(q.sum(), 1):
			raise ValueError("The sum of proportions `q` should be close to 1.")
		cumulative_q = np.cumsum(q)#Calculate the cumulative sum of q
		random_floats = np.random.rand(self.n)#Generate an array of random floats of length n
		random_array = np.searchsorted(cumulative_q, random_floats, side='right')
		#Use np.searchsorted to find the component each random float belongs to
		#Note: side='right' should be used to ensure correct interval division
		random_array += 1 #Since components are numbered starting from 1, add 1 to the resul
		return random_array 

	def get_alternate_array(self):
		#Generate alternating copolymer, only for two components
		indices = np.arange(self.n)
		alternating_array = np.where(indices % 2 == 0, 1, 2)
		return alternating_array

	def get_block_array(self, n_seg, segn, segt):
		#Generate block copolymer
		#n_seg: number of blocks per molecule
		#segn: number of atoms per block
		#segt: atom type name (number) for each block
		block_array = np.zeros(self.n, dtype=int)
		k = 0
		for i in range(n_seg):
			for j in range(segn[i]):
				block_array[k] = segt[i]
				k += 1
		return block_array

	def get_sequential_array(self, n_seq, seqt):
		#Generate sequence block copolymer
		#n_seq: number of sequences per molecule
		#seqt: composition of each sequence
		sequential_array = np.zeros(self.n, dtype=int)
		for i in range(n_seq):
			for j in range(len(seqt)):
				sequential_array[i*len(seqt)+j] = seqt[j]
		return sequential_array
	
	def get_branch_structure(self, n_seg, segn, segt, branch_point, n_branch, t_branch):
		#segn: length of each segment on the main chain
		#segt: monomer type of each segment on the main chain
		#branch_point: monomer type at the branching point
		#n_branch: length of branch chain attached to each branching monomer on the main chain; if 0, no branch at that position; array length equals the number of branch points; if an element is an array, the branch has multiple segments
		#t_branch: monomer type of branch chain attached to each branching monomer on the main chain; array length equals the number of branch points; if an element is an array, the branch has multiple segments
		#Convert SMILES
		for i in range(self.component):
			if i == branch_point - 1:
				self.smi[i] = self.smi[i].replace("[*]","X").replace("*","X")
				print(self.smi[i])
				self.smi[i] = self.smi[i].replace("X","[Ar]",1).replace("X","[He]",1).replace("R","[Ne]",1)
				print(self.smi[i])
			else:
				self.smi[i] = self.smi[i].replace("[*]","X").replace("*","X")
				self.smi[i] = self.smi[i].replace("X","[Ar]",1).replace("X","[He]",1)
		
		#Confirm main chain reaction
		reactants=[AllChem.MolFromSmiles(self.smi[0]), AllChem.MolFromSmiles(self.smi[1])] 
		rxn_str1 ='[Ar:1][c,C,n,N,O,s,S,Si,P:2].[He:3][c,C,n,N,O,s,S,Si,P:4]>>[c,C,n,N,O,s,S,Si,P:2]-[c,C,n,N,O,s,S,Si,P:4].[He:3][Ar:1]'
		rxn_str2 ='[Ar:1]=[c,C,n,N,O,s,S,Si,P:2].[He:3]=[c,C,n,N,O,s,S,Si,P:4]>>[c,C,n,N,O,s,S,Si,P:2]=[c,C,n,N,O,s,S,Si,P:4].[He:3]=[Ar:1]'
		exception_code=0
		try:
			rxn = AllChem.ReactionFromSmarts(rxn_str1)
			products=rxn.RunReactants(reactants)
			mol_test=products[0][0]
			rxn_str = rxn_str1
		except:
			print("Failed to load reaction from rxn_str1\n")
			try:
				rxn = AllChem.ReactionFromSmarts(rxn_str2)
				products=rxn.RunReactants(reactants)
				mol_test=products[0][0]
				rxn_str = rxn_str2
				exception_code = 1
			except:
				print("Failed to load reaction from rxn_str2\n")
				print("Could not load reaction")
				return None
		rxn = AllChem.ReactionFromSmarts(rxn_str)

		#Confirm branching reaction
		reactants=[AllChem.MolFromSmiles(self.smi[branch_point - 1]), AllChem.MolFromSmiles(self.smi[0])] 
		rxn_str1 ='[Ne:1][c,C,n,N,O,s,S,Si,P:2].[He:3][c,C,n,N,O,s,S,Si,P:4]>>[c,C,n,N,O,s,S,Si,P:2]-[c,C,n,N,O,s,S,Si,P:4].[He:3][Ne:1]'
		rxn_str2 ='[Ne:1]=[c,C,n,N,O,s,S,Si,P:2].[He:3]=[c,C,n,N,O,s,S,Si,P:4]>>[c,C,n,N,O,s,S,Si,P:2]=[c,C,n,N,O,s,S,Si,P:4].[He:3]=[Ne:1]'
		exception_code=0
		try:
			rxn_branch = AllChem.ReactionFromSmarts(rxn_str1)
			products = rxn_branch.RunReactants(reactants)
			mol_test = products[0][0]
			rxn_str_branch = rxn_str1
		except:
			print("Failed to load reaction from rxn_str1\n")
			try:
				rxn_branch = AllChem.ReactionFromSmarts(rxn_str2)
				products = rxn_branch.RunReactants(reactants)
				mol_test = products[0][0]
				rxn_str_branch = rxn_str2
				exception_code = 1
			except:
				print("Failed to load reaction from rxn_str2\n")
				print("Could not load reaction")
				return None
		rxn_branch = AllChem.ReactionFromSmarts(rxn_str_branch)

		#Perform reaction
		n_branchpoint = 0
		mol = AllChem.MolFromSmiles(self.smi[segt[0]-1])
		for i in range(n_seg):
			for j in range(segn[i]):
				if i**2 + j**2 != 0:
					#For branch points, build the branch chain first, then attach it entirely to the branch point
					if segt[i] == branch_point:
						mol_branch = AllChem.MolFromSmiles(self.smi[segt[i]-1])
						for k in range(n_branch[n_branchpoint]):
							next_branch_mol = AllChem.MolFromSmiles(self.smi[t_branch[n_branchpoint]-1])
							reactants_branch = (mol_branch, next_branch_mol)
							if k == 0:
								products = rxn_branch.RunReactants(reactants_branch)
							else:
								products = rxn.RunReactants(reactants_branch)
							if products:
								mol_branch = products[0][0]
							else:
								print("Branch Reaction failed at iteration", i)
								return None
						reactants_branch = (mol,mol_branch)
						products = rxn.RunReactants(reactants_branch)
						if products:
							mol = products[0][0]
						else:
							print("Reaction failed at iteration", i)
							return None
						n_branchpoint += 1
					#Connect main chain
					else:
						next_mol = AllChem.MolFromSmiles(self.smi[segt[i]-1])
						reactants_temp = (mol, next_mol)
						products = rxn.RunReactants(reactants_temp)
						if products:
							mol = products[0][0]
						else:
							print("Reaction failed at iteration", i)
							return None

		self.smi = AllChem.MolToSmiles(mol)
		if exception_code:
			end_logo = "[C]"
		else:
			end_logo = "[H]"
		self.smi = self.smi.replace("[He]", end_logo).replace("[Ar]", end_logo)
		chain = AllChem.MolFromSmiles(self.smi)
		chain = AllChem.RemoveHs(chain)
		self.smi = AllChem.MolToSmiles(chain)
		print(self.smi)
		return self.smi