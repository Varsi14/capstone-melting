import numpy as np
import matplotlib.pylab as plt
from pyiron_atomistics import Project
import ase.units as units
from ase.io import write
import math
import pandas as pd
import os
import shutil

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
mor=input("Temperature, manual (m) or range (r):")
if mor=="m":
	m=input("Enter temperature separated by ',':")
	T=np.array(m.split(","),dtype=int)
else:
	up=int(input("Enter upper bound:"))
	low=int(input("Enter lower bound:"))
	ndata=int(input("Enter number of plots:"))
	T=np.linspace(low,up,ndata)

os.makedirs(f"{ele}_{a}x{a}x{a}",exist_ok=True)
shutil.copy("md3plot.py",f"{ele}_{a}x{a}x{a}/")
os.chdir(f"{ele}_{a}x{a}x{a}")
pr = Project(path="hdf5")


basis= pr.create.structure.ase.bulk(ele)
supercell = basis.repeat([a, a, a])
write(f"{ele}.xyz", supercell.to_ase())

for i in T:
	name = f"T_{int(i)}K"
	job = pr.create_job(job_type=pr.job_type.Lammps, job_name=name)
	job.structure = supercell
	job.potential = custom_potential
	job.calc_md(temperature=i, n_ionic_steps=2*1e6, n_print=10, time_step=1)
	job.run(run_mode="queue", delete_existing_job=True)

	

