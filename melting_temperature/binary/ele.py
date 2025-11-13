import os
import mendeleev as men
import yaml
import shutil
import subprocess
import time
import numpy as np
import sys
import glob
from pyiron_atomistics import Project
from ase.io import write

elen1 = sys.argv[1]
elen2 = sys.argv[2]
elen = elen1 + elen2
tm=int(sys.argv[3])
a=float(sys.argv[4])
tpm=200
Temp=int(tm-tpm/2)
id=[0,0]
result=[0,0]
s=True
s2=False
ph=["solid","liquid"]
yfl="input.yaml"
sfl="submit.sh"
pfl="analyse.py"

os.makedirs(f"{elen}",exist_ok=True)
os.makedirs(f"{elen}/{Temp}",exist_ok=True)
shutil.copy(yfl,f"{elen}/{Temp}")
shutil.copy(sfl,f"{elen}/{Temp}")
shutil.copy(pfl,f"{elen}/{Temp}")
os.chdir(f"{elen}/{Temp}")

bas = Project(".").create.structure.ase.crystal((elen1 , elen2), basis=[(0,0,0), (0.5, 0.5, 0.5)],spacegroup=221, cellpar=[a, a, a, 90, 90, 90])
supercell = bas.repeat([5, 5, 5])
write("bas5.data", supercell, format="lammps-data")
supercell = bas.repeat([11, 11, 11])
write("bas11.data", supercell, format="lammps-data")

file=open(yfl,"r")
l=file.readlines()
file.close()

for i in range(len(l)):
	if "element:" in l[i]:
		l[i+1]=f"{l[i+1].rstrip()} {elen1}\n"
		l[i+2]=f"{l[i+2].rstrip()} {elen2}\n"
	elif "lattice:" in l[i]:
		l[i]=f"{l[i].rstrip()} bas5.data \n"
	elif "mass:" in l[i]:
		l[i+1]=f"{l[i+1].rstrip()} {men.element(elen1).atomic_weight}\n"
		l[i+2]=f"{l[i+2].rstrip()} {men.element(elen2).atomic_weight}\n"
	elif "pair_coeff:" in l[i]:
		l[i]=f"{l[i].rstrip()} {elen1} {elen2}'\n"
	elif "temperature:" in l[i]:
		l[i+1]=f"{l[i+1].rstrip()} {Temp}\n"
		l[i+2]=f"{l[i+2].rstrip()} {Temp+tpm}\n"

file=open(yfl,"w")
l=file.writelines(l)
file.close()

result[0] = subprocess.run(["sbatch", sfl, "0"],capture_output=True, text=True)
id[0] = result[0].stdout.strip().split()[-1]
result[1] = subprocess.run(["sbatch", sfl, "1"],capture_output=True, text=True)
id[1] = result[1].stdout.strip().split()[-1]

while s:
	for i in range(len(id)):
		if ("spring" in subprocess.run(["cat", f"ts-bas5.data-{ph[i]}-{Temp}-0/time.out"], capture_output=True, text=True).stdout) or ("spring" in subprocess.run(["cat", f"ts-bas11.data-{ph[i]}-{Temp}-0/time.out"], capture_output=True, text=True).stdout):
			subnbprocess.run(["scancel", id[i]])
		if id[i] not in subprocess.run(["squeue", "-j", id[i]], capture_output=True, text=True).stdout:
			print("nice")
			time.sleep(60)
			if ("temper" not in subprocess.run(["ls", f"ts-bas5.data-{ph[i]}-{Temp}-0/"],capture_output=True,text=True).stdout):
				s2=True
				print(subprocess.run(["ls", f"ts-bas5.data-{ph[i]}-{Temp}-0/"],capture_output=True,text=True))
			else:
				id = np.delete(id, i)
				result = np.delete(result, i)
				ph = np.delete(ph, i)
				if len(id) == 0:
					subprocess.run(["python",pfl])
					time.sleep(20)
					if (int(float(open("temp.dat","r").readlines()[0].strip())) == int(Temp)) or (int(float(open("temp.dat","r").readlines()[0].strip())) == int(Temp+tpm)):
						s2=True
					else:
						Temp1=int(float(open("temp.dat","r").readlines()[0].strip()))-tpm
						Temp2=Temp1+2*tpm
						os.makedirs(f"../11cell",exist_ok=True)
						os.chdir(f"../11cell")
						shutil.copy(f"../{Temp}/bas5.data","bas5.data")
						shutil.copy(f"../{Temp}/bas11.data","bas11.data")
						shutil.copy(f"../../{yfl}",yfl)
						shutil.copy(f"../../{sfl}",sfl)
						shutil.copy(f"../../{pfl}",pfl)
						file=open(yfl,"r")
						l=file.readlines()
						file.close()
						for k in range(len(l)):
							if "element:" in l[k]:
								l[k+1]=f"{l[k+1].rstrip()} {elen1}\n"
								l[k+2]=f"{l[k+2].rstrip()} {elen2}\n"
							elif "lattice:" in l[k]:
								l[k]=f"{l[k].rstrip()} bas5.data \n"
							elif "mass:" in l[k]:
								l[k+1]=f"{l[k+1].rstrip()} {men.element(elen1).atomic_weight}\n"
								l[k+2]=f"{l[k+2].rstrip()} {men.element(elen2).atomic_weight}\n"
							elif "pair_coeff:" in l[k]:
								l[k]=f"{l[k].rstrip()} {elen1} {elen2}'\n"
							elif "temperature:" in l[k]:
								l[k+1]=f"{l[k+1].rstrip()} {Temp}\n"
								l[k+2]=f"{l[k+2].rstrip()} {Temp+tpm}\n"
						file=open(yfl,"w")
						l=file.writelines(l)
						file.close()
						id=[0,0]
						result=[0,0]
						ph=["solid","liquid"]
						result[0] = subprocess.run(["sbatch", sfl, "0"],capture_output=True, text=True)
						id[0] = result[0].stdout.strip().split()[-1]
						result[1] = subprocess.run(["sbatch", sfl, "1"],capture_output=True, text=True)
						id[1] = result[1].stdout.strip().split()[-1]
						time.sleep(60)
						break
					
			if s2:
				if ("System melted" in subprocess.run(["cat", "time.out"], capture_output=True, text=True).stdout) or (os.path.exists("temp.dat") and int(float(open("temp.dat","r").readlines()[0].strip())) == int(Temp)) or ("spring" in subprocess.run(["cat", "time.out"], capture_output=True, text=True).stdout):
					for j in range(len(id)):
						subprocess.run(["scancel", id[j]])
					Tempr=Temp
					Temp=int(Temp-tpm/2)
					if Temp <= 0:
						sys.exit("T below 0")
					os.makedirs(f"../{Temp}",exist_ok=True)
					os.chdir(f"../{Temp}")
					shutil.copy(f"../{Tempr}/bas5.data","bas5.data")
					shutil.copy(f"../{Tempr}/bas11.data","bas11.data")
					shutil.copy(f"../../{yfl}",yfl)
					shutil.copy(f"../../{sfl}",sfl)
					shutil.copy(f"../../{pfl}",pfl)
					file=open(yfl,"r")
					l=file.readlines()
					file.close()
					for k in range(len(l)):
						if "element:" in l[k]:
							l[k+1]=f"{l[k+1].rstrip()} {elen1}\n"
							l[k+2]=f"{l[k+2].rstrip()} {elen2}\n"
						elif "lattice:" in l[k]:
							l[k]=f"{l[k].rstrip()} bas5.data \n"
						elif "mass:" in l[k]:
							l[k+1]=f"{l[k+1].rstrip()} {men.element(elen1).atomic_weight}\n"
							l[k+2]=f"{l[k+2].rstrip()} {men.element(elen2).atomic_weight}\n"
						elif "pair_coeff:" in l[k]:
							l[k]=f"{l[k].rstrip()} {elen1} {elen2}'\n"
						elif "temperature:" in l[k]:
							l[k+1]=f"{l[k+1].rstrip()} {Temp}\n"
							l[k+2]=f"{l[k+2].rstrip()} {Temp+tpm}\n"
					file=open(yfl,"w")
					l=file.writelines(l)
					file.close()
					id=[0,0]
					result=[0,0]
					ph=["solid","liquid"]
					s2=False
					result[0] = subprocess.run(["sbatch", sfl, "0"],capture_output=True, text=True)
					id[0] = result[0].stdout.strip().split()[-1]
					result[1] = subprocess.run(["sbatch", sfl, "1"],capture_output=True, text=True)
					id[1] = result[1].stdout.strip().split()[-1]
					break

				elif ("System solidified" in subprocess.run(["cat", "time.out"], capture_output=True, text=True).stdout) or (os.path.exists("temp.dat") and int(float(open("temp.dat","r").readlines()[0].strip())) == int(Temp+tpm)):
					for j in range(len(id)):
						subprocess.run(["scancel", id[j]])
					Tempr=Temp
					Temp=int(Temp+int(tpm/3)+1)
					os.makedirs(f"../{Temp}",exist_ok=True)
					os.chdir(f"../{Temp}")
					shutil.copy(f"../{Tempr}/bas5.data","bas5.data")
					shutil.copy(f"../{Tempr}/bas11.data","bas11.data")
					shutil.copy(f"../../{yfl}",yfl)
					shutil.copy(f"../../{sfl}",sfl)
					shutil.copy(f"../../{pfl}",pfl)
					file=open(yfl,"r")
					l=file.readlines()
					file.close()
					for k in range(len(l)):
						if "element:" in l[k]:
							l[k+1]=f"{l[k+1].rstrip()} {elen1}\n"
							l[k+2]=f"{l[k+2].rstrip()} {elen2}\n"
						elif "lattice:" in l[k]:
							l[k]=f"{l[k].rstrip()} bas5.data \n"
						elif "mass:" in l[k]:
							l[k+1]=f"{l[k+1].rstrip()} {men.element(elen1).atomic_weight}\n"
							l[k+2]=f"{l[k+2].rstrip()} {men.element(elen2).atomic_weight}\n"
						elif "pair_coeff:" in l[k]:
							l[k]=f"{l[k].rstrip()} {elen1} {elen2}'\n"
						elif "temperature:" in l[k]:
							l[k+1]=f"{l[k+1].rstrip()} {Temp}\n"
							l[k+2]=f"{l[k+2].rstrip()} {Temp+tpm}\n"
					file=open(yfl,"w")
					l=file.writelines(l)
					file.close()
					id=[0,0]
					result=[0,0]
					ph=["solid","liquid"]
					s2=False
					result[0] = subprocess.run(["sbatch", sfl, "0"],capture_output=True, text=True)
					id[0] = result[0].stdout.strip().split()[-1]
					result[1] = subprocess.run(["sbatch", sfl, "1"],capture_output=True, text=True)
					id[1] = result[1].stdout.strip().split()[-1]
					break

			break
		if s:
			time.sleep(20)
