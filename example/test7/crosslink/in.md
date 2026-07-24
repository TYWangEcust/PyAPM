include	"system.in.init"
read_data	"cleanedsystem.data" extra/bond/per/atom 2 extra/angle/per/atom 12 extra/dihedral/per/atom 12 extra/improper/per/atom 12 extra/special/per/atom 12
include	"cleanedsystem.in.settings"

molecule	pre0 pre-molecule0.data
molecule	post0 post-molecule0.data

molecule	pre1 pre-molecule1.data
molecule	post1 post-molecule1.data

molecule	pre2 pre-molecule2.data
molecule	post2 post-molecule2.data

neighbor	2.5 bin
neigh_modify	every 1 delay 0 check yes

timestep	1.0

min_style	cg
minimize	1e-08 1e-10 1000000 10000000
reset_timestep	0

thermo_style	custom step temp epair emol etotal press vol density
thermo	10000

dump	1 all custom 10000 equil.lammpstrj id mol type q x y z ix iy iz
dump_modify	1 sort id

fix	1 all nvt temp 300 300 100
run	100000
unfix	1

fix	2 all npt temp 300 300 100 iso 3000 3000 1000
run	300000
unfix	2

variable	n loop 3
label	here

fix	3 all nvt temp 800 800 100
run	100000
unfix	3

fix	4 all nvt temp 300 300 100
run	100000
unfix	4

fix	5 all npt temp 300 300 100 iso 1000 1000 1000
run	300000
unfix	5

next	n
jump	SELF here

fix	6 all npt temp 300 300 100 iso 1 1 1000
run	1000000
unfix	6
fix	7 all npt temp 300 473.0 100 iso 1 1 1000
run	1000000
unfix	7
undump	1

write_data	equil.data pair ij
write_restart	equil.restart

fix	fxrct all bond/react stabilization yes statted_grp .03 react rxn0 all 1000 1.0 4.0 pre0 post0 automap0.data stabilize_steps 100 prob 0.333 114514 react rxn1 all 1000 1.0 4.0 pre1 post1 automap1.data stabilize_steps 100 prob 0.333 114514 react rxn2 all 1000 1.0 4.0 pre2 post2 automap2.data stabilize_steps 100 prob 0.333 114514 
fix	8 statted_grp_REACT npt temp 473.0 473.0 100 iso 1 1 1000
fix	9 bond_react_MASTER_group temp/rescale 1 473.0 473.0 1 1
thermo_style	custom step temp epair emol etotal press vol density f_fxrct[1]
dump	2 all custom 10000 crosslink.lammpstrj id mol type q x y z ix iy iz
dump_modify	2 sort id

run	500000
write_data	crosslink.data
unfix	8
unfix	9
