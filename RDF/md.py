import numpy as np
import matplotlib.pylab as plt
from pyiron_atomistics import Project
import ase.units as units
from ase.io import write
import math
import pandas as pd

import pandas

pr = Project(path="first_steps")
T=np.linspace(300,1500,6)
jobT=[]


basis= pr.create.structure.ase.bulk("Al", cubic=True)
supercell_3x3x3 = basis.repeat([3, 3, 3])
write("Al.xyz", supercell_3x3x3.to_ase())


for i in T:
	name = f"T_{int(i)}K"
	job = pr.create_job(job_type=pr.job_type.Lammps, job_name=name)
	job.structure = supercell_3x3x3
	job.potential = job.list_potentials()[0]
	job.calc_md(temperature=i, pressure=0, n_ionic_steps=1e4, n_print=10, time_step=1)
	job.run()
	structure1 = job.get_structure(iteration_step=-1)
	write(f"Al_{int(i)}K.xyz", structure1.to_ase())
	plt.plot(job["output/generic/energy_pot"])
	plt.xlabel("step")
	plt.ylabel("Potential energy [eV]");
	plt.savefig(f"{int(i)}K_pot.png")
	plt.close()

	struct = job.get_structure(iteration_step=-1)
	num_neighbors = struct.get_numbers_of_neighbors_in_sphere(cutoff_radius=12).max()
	neighbors = struct.get_neighbors(num_neighbors=12)
	Al_indices = struct.select_index("Al")
	neigh_indices = np.hstack(np.array(neighbors.indices)[Al_indices])
	neigh_distances = np.hstack(np.array(neighbors.distances)[Al_indices])
	Al_neigh_indices = np.in1d(neigh_indices, Al_indices)
	Al_neigh_distances = neigh_distances[Al_neigh_indices]
	traj = job["output/generic/positions"]
	nsteps = len(traj)
	stepincrement = int(nsteps / 10)
	snapshots = range(stepincrement, nsteps - stepincrement, stepincrement)
	bins = np.linspace(0, 6, 200)
	for j in snapshots:
		struct.positions = traj[j]
		neighbors = struct.get_neighbors(num_neighbors)
		neigh_indices = np.hstack(np.array(neighbors.indices)[Al_indices])
		neigh_distances = np.hstack(np.array(neighbors.distances)[Al_indices])
		Al_neigh_indices = np.in1d(neigh_indices, Al_indices)
		Al_neigh_distances = np.concatenate((Al_neigh_distances, neigh_distances[Al_neigh_indices]))
	Al_gr = np.histogram(Al_neigh_distances, bins=bins)
	dr = bins[1] - bins[0]
	normfac = 4 * np.pi * len(struct) / struct.get_volume() * dr * len(snapshots) * len(Al_indices)
	
	plt.bar(Al_gr[1][0:-1], Al_gr[0] / (normfac * (0.5 * (Al_gr[1][1:] + Al_gr[1][:-1])) ** 2), dr)
	plt.title(f"{int(i)}K")
	plt.xlim(min(Al_gr[1][0:-1]),max(Al_gr[1][0:-1]))
	plt.xlabel(r"d$_{Al-Al}$ [$\AA$]")
	plt.ylabel("$g_{Al-Al}(r)$")
	plt.savefig(f"{int(i)}K_d.png")
	plt.close()

