import numpy as np
import subprocess
import shutil
import os
from . import fit
from . import rdf_aa
from . import rdf_cg
import time

def target(temp):
	r, rdf, bdf, adf= rdf_aa.caclulate_distribution(temp)
	bdf_AA = np.array(bdf[0])
	bdf_AB = np.array(bdf[1])
	adf_AAA = np.array(adf[0])
	adf_BAA = np.array(adf[1])
	return r, rdf, bdf_AA, bdf_AB, adf_AAA, adf_BAA

def write_data(n_atoms, n_bonds, n_angles, xsize, ysize, zsize, atom, bond, angle, mass):
	data=open('system.data','w')
	data.write('LAMMPS Description\n\n')
	data.write('%d  atoms\n'%n_atoms)
	data.write('%d  bonds\n'%n_bonds)
	data.write('%d  angles\n'%n_angles)
	data.write('0  dihedrals\n0  impropers\n\n')
	data.write('2  atom types\n2  bond types\n2  angle types\n0  dihedral types\n0  improper types\n\n')
	data.write('0 %f xlo xhi\n0 %f ylo yhi\n0 %f zlo zhi\n\n'%(xsize,ysize,zsize))
	data.write('Masses\n\n1 %f  # A\n2 %f  # B\n'%(mass[0], mass[1]))
	data.write('Atoms  # full\n\n')
	for i in range(n_atoms):
		values = list(atom[i,j] for j in range(0,7))
		data.write('%d %d %d %3.1f %f %f %f\n'%(values[0], values[1], values[2], values[3], values[4], values[5], values[6]))
	data.write('\nBonds\n\n')
	for i in range(n_bonds):
		values = list(bond[i,j] for j in range(0,4))
		data.write('%d %d %d %d\n'%(values[0], values[1], values[2], values[3]))
	data.write('\nAngles\n\n')
	for i in range(n_angles):
		values = list(angle[i,j] for j in range(0,5))
		data.write('%d %d %d %d %d\n'%(values[0], values[1], values[2], values[3], values[4]))
	data.close()
	return 0

def write_table(r, pot_non_AA, pot_non_AB, pot_non_BB, press_factor):
	pot_non_AA = pot_non_AA - press_factor*(1-r/r[-1])
	pot_non_AB = pot_non_AB - press_factor*(1-r/r[-1])
	pot_non_BB = pot_non_BB - press_factor*(1-r/r[-1])
	pot_non_AA = fit.tail_correction(r, pot_non_AA,r[-10])
	pot_non_AB = fit.tail_correction(r, pot_non_AB,r[-10])
	pot_non_BB = fit.tail_correction(r, pot_non_BB,r[-10])
	r_AA, pot_non_AA, force_non_AA = fit.diff(pot_non_AA,r)
	r_AB, pot_non_AB, force_non_AB = fit.diff(pot_non_AB,r)
	r_BB, pot_non_BB, force_non_BB = fit.diff(pot_non_BB,r)
	table_AA = open('ljAA.table','w')
	table_AA.write('LJ_AA\n')
	table_AA.write('N %d\n\n'%len(r))
	for i in range(len(r_AA)):
		table_AA.write('%d %.4e %f %f\n'%(i, r_AA[i], pot_non_AA[i], force_non_AA[i]))
	if len(r_AA)<len(r):
		for i in range(len(r_AA), len(r)):
			table_AA.write('%d %.4e %f %f\n'%(i, r[i], 0, 0))
	table_AA.close()
	table_AB = open('ljAB.table','w')
	table_AB.write('LJ_AB\n')
	table_AB.write('N %d\n\n'%len(r))
	for i in range(len(r_AB)):
		table_AB.write('%d %.4e %f %f\n'%(i, r_AB[i], pot_non_AB[i], force_non_AB[i]))
	if len(r_AB)<len(r):
		for i in range(len(r_AB), len(r)):
			table_AB.write('%d %.4e %f %f\n'%(i, r[i], 0, 0))
	table_AB.close()
	table_BB = open('ljBB.table','w')
	table_BB.write('LJ_BB\n')
	table_BB.write('N %d\n\n'%len(r))
	for i in range(len(r_BB)):
		table_BB.write('%d %.4e %f %f\n'%(i, r_BB[i], pot_non_BB[i], force_non_BB[i]))
	if len(r_BB)<len(r):
		for i in range(len(r_BB), len(r)):
			table_BB.write('%d %.4e %f %f\n'%(i, r[i], 0, 0))
	table_BB.close()
	return 0

def write_lammps(temp, state):
	in_file = open('in.md','w')
	in_file.write('atom_style\tfull\nunits\treal\nboundary\tp p p\ntimestep\t1.0\n\n')
	in_file.write('neighbor\t3.0 bin\nneigh_modify\tdelay 0 every 1 check yes one 5000 page 50000\n\n')
	in_file.write('special_bonds\tlj  0  0  1\n')
	in_file.write('read_data\tsystem.data\n')
	in_file.write('include\tsystem.in.settings\n')
	in_file.write('timestep\t1.0\n')
	in_file.write('thermo_style\tcustom step temp epair emol etotal press vol density\nthermo\t10000\n\n')
	in_file.write('velocity\tall create %d 4928459 dist gaussian\n'%temp)
	in_file.write('fix\t1 all nve/limit 0.001\ndump\t1 all atom 100000 relax.lammpstrj\nrun\t400000\nunfix\t1\nundump\t1\n\n')
	in_file.write('velocity\tall create %d 4928459 dist gaussian\n'%temp)
	if state == 'nvt':
		in_file.write('compute\tmypress all pressure thermo_temp\n')
		in_file.write('fix\t2 all ave/time 1000 200 1000000 c_mypress file press.profile\n')
		in_file.write('fix\t3 all nvt temp %d %d $(100*dt)\n'%(temp,temp))
		in_file.write('dump\t2 all atom 100000 %d-normalizing.lammpstrj\n'%temp)
		in_file.write('run\t1000000\nunfix\t2\nunfix\t3\nundump\t2\nwrite_data\t%d-normalizing.data\n'%temp)
	elif state == 'npt':
		in_file.write('variable\td1 equal density\n')
		in_file.write('fix\t2 all ave/time 1000 200 1000000 v_d1 file density.profile\n')
		in_file.write('fix\t3 all npt temp %d %d $(100*dt) iso 1.0 1.0 $(1000*dt)\n'%(temp,temp))
		in_file.write('dump\t2 all atom 100000 %d-normalizing.lammpstrj\n'%temp)
		in_file.write('run\t1000000\nunfix\t2\nunfix\t3\nundump\t2\nwrite_data\t%d-normalizing.data\n'%temp)
	in_file.close()
	return 0

def write_settings(n_peak_bond_AA, pram_bond_AA, n_peak_angle_AAA, pram_angle_AAA, n_peak_bond_AB, pram_bond_AB, n_peak_angle_BAA, pram_angle_BAA):
	settings = open('system.in.settings','w')
	settings.write('pair_style\ttable linear 101\n')
	settings.write('pair_coeff\t1 1 ljAA.table LJ_AA\n')
	settings.write('pair_coeff\t1 2 ljAB.table LJ_AB\n')
	settings.write('pair_coeff\t2 2 ljBB.table LJ_BB\n')
	settings.write('bond_style\tgaussian\n')
	settings.write('bond_coeff\t1 %d %d '%(temp, n_peak_bond_AA))
	for i in range(n_peak_bond_AA):
		settings.write('%f %f %f '%(pram_bond_AA[3*i], pram_bond_AA[3*i+1], pram_bond_AA[3*i+2]))
	settings.write('\n')
	settings.write('bond_coeff\t2 %d %d '%(temp, n_peak_bond_AB))
	for i in range(n_peak_bond_AB):
		settings.write('%f %f %f '%(pram_bond_AB[3*i], pram_bond_AB[3*i+1], pram_bond_AB[3*i+2]))
	settings.write('\n')
	settings.write('angle_style\tgaussian\n')
	settings.write('angle_coeff\t1 %d %d '%(temp, n_peak_angle_AAA))
	for i in range(n_peak_angle_AAA):
		settings.write('%f %f %f '%(pram_angle_AAA[3*i], pram_angle_AAA[3*i+1], pram_angle_AAA[3*i+2]))
	settings.write('\n')
	settings.write('angle_coeff\t2 %d %d '%(temp, n_peak_angle_BAA))
	for i in range(n_peak_angle_BAA):
		settings.write('%f %f %f '%(pram_angle_BAA[3*i], pram_angle_BAA[3*i+1], pram_angle_BAA[3*i+2]))
	settings.write('\n')
	settings.close()
	return 0

def optimization_IBI(temp, n_step, path):
	print('Optimizing!\n')
	os.chdir('origin')
	r, rdf_target, bdf_AA_target, bdf_AB_target, adf_AAA_target, adf_BAA_target = target(temp)
	n_peak_bond_AA, pram_bond_AA = fit.fitting_gaussian_bond(bdf_AA_target[:,0],bdf_AA_target[:,1])
	n_peak_bond_AB, pram_bond_AB = fit.fitting_gaussian_bond(bdf_AB_target[:,0] ,bdf_AB_target[:,1])
	n_peak_angle_AAA, pram_angle_AAA = fit.fitting_gaussian_angle(adf_AAA_target[:,0],adf_AAA_target[:,1])
	n_peak_angle_BAA, pram_angle_BAA = fit.fitting_gaussian_angle(adf_BAA_target[:,0],adf_BAA_target[:,1])
	n_atoms, n_bonds, n_angles, xsize, ysize, zsize, atom, bond, angle, mass = rdf_aa.get_data(temp)
	pot_non_AA = fit.Boltzmann_inversion(rdf_target[0], temp)
	pot_non_AB = fit.Boltzmann_inversion(rdf_target[1], temp)
	pot_non_BB = fit.Boltzmann_inversion(rdf_target[2], temp)
	is_inf = np.isinf(pot_non_AA)
	inf_index = np.where(is_inf)[0][-1]
	pot_non_AA = fit.linear_head_correction(r, pot_non_AA, inf_index)
	is_inf = np.isinf(pot_non_AB)
	inf_index = np.where(is_inf)[0][-1]
	pot_non_AB = fit.linear_head_correction(r, pot_non_AB, inf_index)
	is_inf = np.isinf(pot_non_BB)
	inf_index = np.where(is_inf)[0][-1]
	pot_non_BB = fit.linear_head_correction(r, pot_non_BB, inf_index)
	os.chdir('..')
	if not os.path.exists('opt'):
		os.mkdir('opt')
	os.chdir('opt')
	if not os.path.exists('Iteration'):
		os.mkdir('Iteration')
	os.chdir('Iteration')
	write_data(n_atoms, n_bonds, n_angles, xsize, ysize, zsize, atom, bond, angle, mass)
	os.chdir('..')
	i = 0
	fit_func = 0
	optlog = open('optlog.txt','w')
	optlog.write('Iter\tFitting value\n')
	optlog.close()
	press_factor = 0.14
	density_aa = 1.04976
	if n_step > 1:
		while (fit_func < 0.9995) and (i < n_step):
			print('*********************Iteration %d!*************************\n'%i)
			i += 1
			shutil.copytree('Iteration','Iteration_%d'%i)
			os.chdir('Iteration_%d'%i)
			if i < int(n_step/2):
				write_lammps(temp, state = 'nvt')
				write_settings(n_peak_bond_AA, pram_bond_AA, n_peak_angle_AAA, pram_angle_AAA, n_peak_bond_AB, pram_bond_AB, n_peak_angle_BAA, pram_angle_BAA)
				rdf, pot_non_AA, pot_non_AB, pot_non_BB, bdf_AA, bdf_AB, adf_AAA, adf_BAA = iteration(r, rdf_target, pot_non_AA, pot_non_AB, pot_non_BB, temp, press_factor, alpha = 0.05)
				if i % 5 == 0 and i <= int(n_step/3):
					n_peak_bond_AA, pram_bond_AA = iteration_bond(bdf_AA_target, bdf_AA, n_peak_bond_AA, pram_bond_AA)
					n_peak_bond_AB, pram_bond_AB = iteration_bond(bdf_AB_target, bdf_AB, n_peak_bond_AB, pram_bond_AB)
					n_peak_angle_AAA, pram_angle_AAA = iteration_angle(adf_AAA_target, adf_AAA, n_peak_angle_AAA, pram_angle_AAA)
					n_peak_angle_BAA, pram_angle_BAA = iteration_angle(adf_BAA_target, adf_BAA, n_peak_angle_BAA, pram_angle_BAA)
			else:
				write_lammps(temp, state = 'npt')
				write_settings(n_peak_bond_AA, pram_bond_AA, n_peak_angle_AAA, pram_angle_AAA, n_peak_bond_AB, pram_bond_AB, n_peak_angle_BAA, pram_angle_BAA)
				rdf, pot_non_AA_new, pot_non_AB_new, pot_non_BB_new, bdf_AA, bdf_AB, adf_AAA, adf_BAA = iteration(r, rdf_target, pot_non_AA, pot_non_AB, pot_non_BB, temp, press_factor, alpha = 0.05)
				with open('density.profile') as density_file:
					line = density_file.readlines()[-1]
					density = float(line.split()[-1])
				press_factor = density_correction(density, density_aa, press_factor)
			fit_func = (fit.match(rdf[0], rdf_target[0])+fit.match(rdf[1], rdf_target[1])+fit.match(rdf[2], rdf_target[2]))/3
			os.chdir('..')
			print('*********************Finish Iteration %d!*************************\n'%i)
			optlog = open('optlog.txt','a')
			optlog.write('%d\t%f\t%f\n'%(i,press_factor,fit_func))
			optlog.close()
			print('Fitting value: %f\n'%fit_func)
			print("<(` ^')>")
	print('************Optmize Successfully!****************\n')
	bondtype = [[1,1,1],[2,1,2]]
	angletype = [[1,1,1,1],[2,2,1,1]]
	bondcoeff = [pram_bond_AA.insert(0,n_peak_bond_AA), pram_bond_AB.insert(0,n_peak_bond_AB)]
	anglecoeff = [pram_angle_AAA.insert(0,n_peak_angle_AAA), pram_angle_BAA.insert(0,n_peak_angle_BAA)]
	return bondtype, angletype, bondcoeff, anglecoeff, mass

def iteration(r, rdf_target, pot_non_AA, pot_non_AB, pot_non_BB, temp, press_factor, alpha):
	write_table(r, pot_non_AA, pot_non_AB, pot_non_BB, press_factor)
	returncode = run_lammps(temp)
	if returncode != 0:
		print('Lammps Error!\n')
	r, rdf, bdf, adf = rdf_cg.caclulate_distribution(temp)	
	pot_AA_next = fit.iteration_Boltzmann(pot_non_AA, rdf[0], rdf_target[0], temp, alpha)
	pot_AA_next = fit.head_correction(r, pot_AA_next, pot_non_AA, form='linear')
	pot_AB_next = fit.iteration_Boltzmann(pot_non_AB, rdf[1], rdf_target[1], temp, alpha)
	pot_AB_next = fit.head_correction(r, pot_AB_next, pot_non_AB, form='linear')
	pot_BB_next = fit.iteration_Boltzmann(pot_non_BB, rdf[2], rdf_target[2], temp, alpha)
	pot_BB_next = fit.head_correction(r, pot_BB_next, pot_non_BB, form='linear')
	bdf_AA = np.array(bdf[0])
	bdf_AB = np.array(bdf[1])
	adf_AAA = np.array(adf[0])
	adf_BAA = np.array(adf[1])
	return rdf, pot_AA_next, pot_AB_next, pot_BB_next, bdf_AA, bdf_AB, adf_AAA, adf_BAA

def iteration_bond(bdf_target, bdf, n_peak_bond, pram_bond):
	bdf_next = np.zeros(len(bdf[:,0]))
	constant = np.sqrt(np.pi / 2)
	for i in range(len(bdf[:,0])):
		if bdf[i,1] != 0:
			factor = bdf_target[i,1] / bdf[i,1]
		elif bdf[i,1] == 0:
			if bdf_target[i,1] == 0:
				factor = 0
			else:
				factor = 1
		bdf_func = 0
		for j in range(n_peak_bond):
			bdf_func += pram_bond[3*j]/(pram_bond[3*j+1] * constant) * np.exp(-2 * ((bdf[i,0] - pram_bond[3*j+2]) / pram_bond[3*j+1]) ** 2)
		bdf_next[i] = bdf_func * (factor ** 0.2)
	n_peak_bond, pram_bond = fit.fitting_gaussian_bond(bdf[:,0],bdf_next)
	return n_peak_bond, pram_bond

def iteration_angle(adf_target, adf, n_peak_angle, pram_angle):
	adf_next = np.zeros(len(adf[:,0]))
	constant = np.sqrt(np.pi / 2)
	for i in range(len(adf[:,0])):
		if adf[i,1] != 0:
			factor = adf_target[i,1] / adf[i,1]
		elif adf[i,1] == 0:
			if adf_target[i,1] == 0:
				factor = 0
			else:
				factor = 1
		adf_func = 0
		for j in range(n_peak_angle):
			adf_func += pram_angle[3*j]/(pram_angle[3*j+1] * constant) * np.exp(-2 * ((np.deg2rad(adf[i,0]) - np.deg2rad(pram_angle[3*j+2])) / pram_angle[3*j+1]) ** 2)
		adf_next[i] = adf_func * (factor ** 0.1)
	n_peak_angle, pram_angle = fit.fitting_gaussian_angle(adf[:,0],adf_next)
	return n_peak_angle, pram_angle

def press_iteration(r, pot_non_AA, pot_non_BB, temp, press_factor):
	write_table(r, pot_non_AA, pot_non_BB, press_factor)
	returncode = run_lammps(temp)
	if returncode != 0:
		print('Lammps Error!\n')
	r, rdf, bdf, adf = rdf_cg.caclulate_distribution(temp)
	press_file = open('press.profile')
	line = press_file.readlines()[-1]
	press = float(line.split()[-1])
	press_file.close()
	return pot_non_AA, pot_non_BB, press, rdf

def density_correction(density, density_aa, factor):
	factor_new = factor * ((density_aa/density))
	if factor_new > (factor+0.1):
		factor_new = factor+0.1
	elif factor_new < (factor-0.1):
		factor_new = factor-0.1
	else:
		factor_new = factor_new
	return factor_new

def run_lammps(temp):
	lmp_exec = shutil.which('lmp_mpi') or os.environ.get('LAMMPS_EXEC', 'lmp_mpi')
	command = ['mpirun', '-np', '8', lmp_exec, '-in', 'in.md']
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	while not os.path.exists('%d-normalizing.data'%temp):
		time.sleep(60)
		pass
	time.sleep(60)	
	return process.wait()


