# capstone-melting

The Python scripts in the RDF folder can be used to simulate the interaction of the atoms in an a x b x c supercell.

The mdrun.py simulates the trajectory of the atoms. While running the code, it will ask to enter the supercell size (a x b x c), the element and the temperatures. The temperature can be enter manually (m) or a range (r) and number of temperature be chosen. 

The mdplot.py transforms the trajectory to a radial distribution function plots it for every temperature seperatly. It also combines every plots and displace in 2D, 3D and contour.

In the folder melting_temperature there are folders for the elements and on for binary systems. 

For the elements, if the phase at room temperature is the same as the phase bordering the liquid phase, the files in no_lattice can be used. You must run the ele.sh and it will ask you which element you want to analyse. It will optimize the temperature range to search with a supercell 5x5x5 and then find the temperature with a supercell of 11x11x11. 

If the phase at room temperature is different than the phase bordering the liquid phase, then you can use the files in the file lattice. If you run the ele.sh, it will ask for the element, crystal structure and lattice constant a. It will run as described before.

To analyse binary systems, you must run the files in binary. If you run the ele.sh, it will ask for the first element, the second element, the melting temperature and the lattice constant a. If you want to analyse something else than cubic systems, you can add its lattice constant. To construct the supercell you must edit the line 36. There you must edit the Wyckoff parameters, spacegroup and if necessary the lattice parameter.