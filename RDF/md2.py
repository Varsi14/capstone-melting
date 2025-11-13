import numpy as np
import matplotlib.pylab as plt
from pyiron_atomistics import Project
import ase.units as units
from ase.io import write
import math
import pandas as pd
import os

custom_potential = pd.DataFrame({
    'Name': ['grace_fs_NiAlH'],
    'Filename': [[]],               # No external LAMMPS potential file (we point to a YAML in pair_coeff)
    'Model': ['Custom'],
    'Species': [['Ni', 'Al', 'H']],
    'Config': [[
        'atom_style atomic\n',
        'pair_style grace/fs\n',
        'pair_coeff * * /cluster/home/theorystud08/.cache/grace/GRACE-FS-OMA.yaml Ni Al H\n'
    ]]
})

a=int(input("Enter cell size:"))
ele=input(f"Enter element({custom_potential['Species'][0]}):")
temp=round(int(input("Enter melting temperature:")),-2)
os.makedirs(f"{ele}_{a}x{a}x{a}_1",exist_ok=True)
os.chdir(f"{ele}_{a}x{a}x{a}_1")
pr = Project(path=f"{ele}_{a}x{a}x{a}")

T=np.linspace(temp-600,temp+600,3)
basis= pr.create.structure.ase.bulk(ele, cubic=True)
supercell = basis.repeat([a, a, a])
write(f"{ele}.xyz", supercell.to_ase())


for i in T:
	name = f"T_{int(i)}K"
	job = pr.create_job(job_type=pr.job_type.Lammps, job_name=name)
	job.structure = supercell
	job.potential = custom_potential
	job.calc_md(temperature=i, pressure=0, n_ionic_steps=1e4, n_print=10, time_step=1)
	job.run(run_mode="queue", delete_existing_job=True) 
	pr.wait_for_job(job)
	structure1 = job.get_structure()
	write(f"{ele}_{int(i)}K.xyz", structure1.to_ase())
	plt.plot(job["output/generic/energy_pot"])
	plt.xlabel("step")
	plt.ylabel("Potential energy [eV]");
	plt.savefig(f"{ele}_{int(i)}K_pot.png")
	plt.close()

	struct = job.get_structure()
	num_neighbors = struct.get_numbers_of_neighbors_in_sphere(cutoff_radius=12).max()
	neighbors = struct.get_neighbors(num_neighbors=12)
	ele_indices = struct.select_index(ele)
	neigh_indices = np.hstack(np.array(neighbors.indices)[ele_indices])
	neigh_distances = np.hstack(np.array(neighbors.distances)[ele_indices])
	ele_neigh_indices = np.in1d(neigh_indices, ele_indices)
	ele_neigh_distances = neigh_distances[ele_neigh_indices]
	traj = job["output/generic/positions"]
	nsteps = len(traj)
	stepincrement = int(nsteps / 10)
	snapshots = range(stepincrement, nsteps - stepincrement, stepincrement)
	bins = np.linspace(0, 6, 200)
	for j in snapshots:
		struct.positions = traj[j]
		neighbors = struct.get_neighbors(num_neighbors)
		neigh_indices = np.hstack(np.array(neighbors.indices)[ele_indices])
		neigh_distances = np.hstack(np.array(neighbors.distances)[ele_indices])
		ele_neigh_indices = np.in1d(neigh_indices, ele_indices)
		ele_neigh_distances = np.concatenate((ele_neigh_distances, neigh_distances[ele_neigh_indices]))
	ele_gr = np.histogram(ele_neigh_distances, bins=bins)
	dr = bins[1] - bins[0]
	normfac = 4 * np.pi * len(struct) / struct.get_volume() * dr * len(snapshots) * len(ele_indices)
	
	plt.bar(ele_gr[1][0:-1], ele_gr[0] / (normfac * (0.5 * (ele_gr[1][1:] + ele_gr[1][:-1])) ** 2), dr)
	plt.title(f"{int(i)}K")
	plt.xlim(min(ele_gr[1][0:-1]),max(ele_gr[1][0:-1]))
	plt.xlabel(r"d$_{Al-Al}$ [$\AA$]")
	plt.ylabel("$g_{Al-Al}(r)$")
	plt.savefig(f"{ele}_{int(i)}K_d.png")
	plt.close()

