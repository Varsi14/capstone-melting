import numpy as np
import matplotlib.pylab as plt
from pyiron_atomistics import Project
import ase.units as units
from ase.io import write
import math
import pandas as pd
import os
from pyiron_atomistics.atomistics.job.atomistic import Trajectory
from ase import Atoms
from matplotlib import cm
import glob


name= os.path.basename(os.path.dirname(os.path.abspath(__file__)))
ele= name.split("_")[0]
a= name.split("_")[1].split("x")[0]
b= name.split("_")[1].split("x")[1]
c= name.split("_")[1].split("x")[2]

os.chdir("hdf5")
T=[]
fe=[]
x=[]
y=[]
z=[]
for i in os.listdir():
	if os.path.isdir(i):
		T.append(int(i.split("_")[1].split("K")[0]))
os.chdir("..")
pr=Project(path="hdf5")
T.sort()
#T=T[0:4]
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
	plt.xlabel(f"r d$_{{{ele}-{ele}}}$ [$\\AA$]")
	plt.ylabel(f"$g_{{{ele}-{ele}}}(r)$")
	plt.savefig(f"{ele}_{int(i)}K_d.png")
	plt.close()
	traj = job["output/generic/positions"]
	symbols = job["output/generic/atoms"]
	ase_traj = []
	for j in range(len(traj)):
		atoms = Atoms(symbols=symbols,positions=traj[j],pbc=False)
		ase_traj.append(atoms)
	write(f"{ele}_{int(i)}K_gif.xyz",ase_traj)
	x.append(ele_gr[1][0:-1])
	y.append(ele_gr[0] / (normfac * (0.5 * (ele_gr[1][1:] + ele_gr[1][:-1])) ** 2))
	z.append(i)
	print(i)

for i in range(len(T)):
	plt.plot(x[i],y[i],color=(0,0,0.1+i/len(T)*0.9,0.8))
plt.axhline(y=1, color='red', linestyle='--', linewidth=0.5)
plt.xlabel(f"r d$_{{{ele}-{ele}}}$ [$\\AA$]")
plt.ylabel(f"$g_{{{ele}-{ele}}}(r)$")
plt.savefig(f"{ele}_2D.png")


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for i in range(len(T)):
	ax.plot(x[i],y[i],zs=z[i],color=(0,0,0.1+i/len(T)*0.9,0.8))
ax.set_xlabel(f"r d$_{{{ele}-{ele}}}$ [$\\AA$]")
ax.set_ylabel(f"$g_{{{ele}-{ele}}}(r)$")
ax.set_zlabel("Temperature [K]")
ax.view_init(vertical_axis='y',azim=45)
plt.savefig(f"{ele}_3D.png")
plt.close()

x,z= np.meshgrid(x[0], z)

plt.figure()
cp = plt.contourf(x,z,y, cmap=cm.viridis, levels=100)
plt.xlabel(f"r d$_{{{ele}-{ele}}}$ [$\\AA$]")
plt.ylabel("Temperature [K]")
plt.colorbar(cp, label=f"$g_{{{ele}-{ele}}}(r)$")
plt.savefig(f"{ele}_contour.png")
plt.close()