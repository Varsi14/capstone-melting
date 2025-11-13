import os
import mendeleev as men
import yaml
import shutil
import subprocess
import time
import numpy as np
import sys
import glob

elen = sys.argv[1]


ele=men.element(elen)
struct=ele.lattice_structure

if struct == "HEX":
	struct="HCP"
if struct not in ["BCC","FCC","HCP","SC","DIA"]:
	sys.exit(struct)
struct=struct.lower()
const=ele.lattice_constant
tm=int(round(int(ele.melting_point),-2))
if tm<2000:
	tpm=int(tm/10)
else:
	tpm=200
Temp=int(tm-tpm/2)
mass=ele.atomic_weight
id=[0,0]
result=[0,0]
a=5
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

file=open(yfl,"r")
l=file.readlines()
file.close()

for i in range(len(l)):
	if "element:" in l[i]:
		l[i]=f"{l[i].rstrip()} {elen}\n"
	elif "lattice:" in l[i]:
		l[i]=f"{l[i].rstrip()} {struct}\n"
	elif "lattice_constant:" in l[i]:
		l[i]=f"{l[i].rstrip()} {const}\n"
	elif "mass:" in l[i]:
		l[i]=f"{l[i].rstrip()} {mass}\n"
	elif "pair_coeff:" in l[i]:
		l[i]=f"{l[i].rstrip()} {elen}'\n"
	elif "repeat:" in l[i]:
		l[i+1]=f"{l[i+1].rstrip()} {a}\n"
		l[i+2]=f"{l[i+2].rstrip()} {a}\n"
		l[i+3]=f"{l[i+3].rstrip()} {a}\n"
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
		if id[i] not in subprocess.run(["squeue", "-j", id[i]], capture_output=True, text=True).stdout:
			print("nice")
			time.sleep(60)
			if "temper" not in subprocess.run(["ls", f"ts-{struct}-{ph[i]}-{Temp}-0/"],capture_output=True,text=True).stdout:
				s2=True
				print(subprocess.run(["ls", f"ts-{struct}-{ph[i]}-{Temp}-0/"],capture_output=True,text=True))
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
						a=11
						os.makedirs(f"../11cell",exist_ok=True)
						os.chdir(f"../11cell")
						shutil.copy(f"../../{yfl}",yfl)
						shutil.copy(f"../../{sfl}",sfl)
						shutil.copy(f"../../{pfl}",pfl)
						file=open(yfl,"r")
						l=file.readlines()
						file.close()
						for k in range(len(l)):
							if "element:" in l[k]:
								l[k]=f"{l[k].rstrip()} {elen}\n"
							elif "lattice:" in l[k]:
								l[k]=f"{l[k].rstrip()} {struct}\n"
							elif "lattice_constant:" in l[k]:
								l[k]=f"{l[k].rstrip()} {const}\n"
							elif "mass:" in l[k]:
								l[k]=f"{l[k].rstrip()} {mass}\n"
							elif "pair_coeff:" in l[k]:
								l[k]=f"{l[k].rstrip()} {elen}'\n"
							elif "repeat:" in l[k]:
								l[k+1]=f"{l[k+1].rstrip()} {a}\n"
								l[k+2]=f"{l[k+2].rstrip()} {a}\n"
								l[k+3]=f"{l[k+3].rstrip()} {a}\n"
							elif "temperature:" in l[k]:
								l[k+1]=f"{l[k+1].rstrip()} {Temp1}\n"
								l[k+2]=f"{l[k+2].rstrip()} {Temp2}\n"
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
						while True:
							if id[0] not in subprocess.run(["squeue", "-j", id[0]], capture_output=True, text=True).stdout:
								if id[1] not in subprocess.run(["squeue", "-j", id[1]], capture_output=True, text=True).stdout:
									break
						time.sleep(60)
						subprocess.run(["python",pfl])
						sys.exit("Job done")
						break
					
			if s2:
				if ("System melted" in subprocess.run(["cat", "time.out"], capture_output=True, text=True).stdout) or (os.path.exists("temp.dat") and int(float(open("temp.dat","r").readlines()[0].strip())) == int(Temp)):
					for j in range(len(id)):
						subprocess.run(["scancel", id[j]])
					Temp=int(Temp-tpm/2)
					if Temp <= 0:
						sys.exit("T below 0")
					os.makedirs(f"../{Temp}",exist_ok=True)
					os.chdir(f"../{Temp}")
					shutil.copy(f"../../{yfl}",yfl)
					shutil.copy(f"../../{sfl}",sfl)
					shutil.copy(f"../../{pfl}",pfl)
					file=open(yfl,"r")
					l=file.readlines()
					file.close()
					for k in range(len(l)):
						if "element:" in l[k]:
							l[k]=f"{l[k].rstrip()} {elen}\n"
						elif "lattice:" in l[k]:
							l[k]=f"{l[k].rstrip()} {struct}\n" 
						elif "lattice_constant:" in l[k]:
							l[k]=f"{l[k].rstrip()} {const}\n"
						elif "mass:" in l[k]:
							l[k]=f"{l[k].rstrip()} {mass}\n"
						elif "pair_coeff:" in l[k]:
							l[k]=f"{l[k].rstrip()} {elen}'\n"
						elif "repeat:" in l[k]:
							l[k+1]=f"{l[k+1].rstrip()} {a}\n"
							l[k+2]=f"{l[k+2].rstrip()} {a}\n"
							l[k+3]=f"{l[k+3].rstrip()} {a}\n"
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
					Temp=int(Temp+int(tpm/3)+1)
					os.makedirs(f"../{Temp}",exist_ok=True)
					os.chdir(f"../{Temp}")
					shutil.copy(f"../../{yfl}",yfl)
					shutil.copy(f"../../{sfl}",sfl)
					shutil.copy(f"../../{pfl}",pfl)
					file=open(yfl,"r")
					l=file.readlines()
					file.close()
					for k in range(len(l)):
						if "element:" in l[k]:
							l[k]=f"{l[k].rstrip()} {elen}\n"
						elif "lattice:" in l[k]:
							l[k]=f"{l[k].rstrip()} {struct}\n"
						elif "lattice_constant:" in l[k]:
							l[k]=f"{l[k].rstrip()} {const}\n"
						elif "mass:" in l[k]:
							l[k]=f"{l[k].rstrip()} {mass}\n"
						elif "pair_coeff:" in l[k]:
							l[k]=f"{l[k].rstrip()} {elen}'\n"
						elif "repeat:" in l[k]:
							l[k+1]=f"{l[k+1].rstrip()} {a}\n"
							l[k+2]=f"{l[k+2].rstrip()} {a}\n"
							l[k+3]=f"{l[k+3].rstrip()} {a}\n"
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
