import numpy as np
import matplotlib.pyplot as plt
import glob

st, sfe, sferr = np.loadtxt(glob.glob("ts-*-solid-*-0/temperature_sweep.dat")[0], unpack=True)
lt, lfe, lferr = np.loadtxt(glob.glob("ts-*-liquid-*-0/temperature_sweep.dat")[0], unpack=True)
plt.plot(st, sfe, color="#E53935", label="solid")
plt.plot(lt, lfe, color="#0097A7", label="liquid")
plt.xlabel("Temperature (K)", fontsize=12)
plt.ylabel("F (ev/atom)", fontsize=12)
plt.legend()
plt.savefig("tm.png", dpi=300, bbox_inches="tight")
args = np.argsort(np.abs(sfe-lfe))
np.savetxt('temp.dat',[st[args[0]]], fmt='%.2f')
