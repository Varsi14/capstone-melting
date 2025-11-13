import numpy as np
import matplotlib.pylab as plt
from pyiron_atomistics import Project
import ase.units as units
from ase.io import write
import math
import pandas as pd
import os

name= os.path.basename(os.path.dirname(os.path.abspath(__file__)))
ele= name.split("_")[0]
a= name.split("x")[0]
os.chdir("hdf5")
T=[]
fe=[]
for i in os.listdir():
	if os.path.isdir(i):
		T.append(int(i.split("_")[1].split("K")[0]))
os.chdir("..")
pr=Project(path="hdf5")

for i in T:
	job=pr.load(f"T_{int(i)}K")
	struct = job.get_structure()
	write(f"{ele}_{int(i)}K.xyz", struct.to_ase())
	plt.plot(job["output/generic/energy_pot"])
	plt.xlabel("step")
	plt.ylabel("Potential energy [eV]");
	plt.savefig(f"{ele}_{int(i)}K_pot.png")
	plt.close()
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
	fe.append(np.mean(job["output/generic/energy_pot"][1]))
plt.plot(T,fe, marker="o", linestyle="none")
plt.xlabel("Temperature")
plt.ylabel("free Energy [eV]");
plt.savefig("Free_Energy.png")
plt.close()

	
