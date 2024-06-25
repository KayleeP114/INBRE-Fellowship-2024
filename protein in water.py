from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

#uploading pdb file
pdb = PDBFile('1AKI.pdb')

#specify the forcefield(reproduces molecular geometry and selected properties)
forcefield = ForceField('amber14-all.xml','amber14/tip3pfb.xml')

#clean up(removing crystal water structures and adding in H)
modeller = Modeller(pdb.topology, pdb.positions)
modeller.deletewater()
residues=modeller.addHydrogens(forcefield)

#adding solvate(simulates realistic environment)
modeller.addSolvent(forcefield, padding=1.0*nanometer)

#setting up the molecular system(this includes specifying the atoms, their positions, and their interactions with eachother)
#setting up integrator(determines how the system evolves throughout the simulation)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulaiton = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

#local energy minimization(helps us to find energetically favorable configurations)
print('Minimizing energy')
simulation.minimizeEnergy()

#setup reporting
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simultion.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(StateDataReport("md_log.txt", 100, step=True, potentialEnergy=True, temperature=True, volume=True))

#NVT equillibration(equilibrate temperature by running the sim a certain number of times)
print("Running NVT")
simulation.step(10000)

#NPT prodution MD(adding in a barostat to control the pressure in NPT ensemble)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
simulation.context.reinitialize(preserveState=True)

print("Running NPT")
simulation.step(10000)

#Analysis(making plots and setting axis)
import numpy as np
import matplotlib.pyplot as plt 
data = np.loadtxt("md_log.txt", delimiter=',')

step = data[:,0]
potential_energy = data[:,1]
temperature = data[:,2]
volume = data[:,3]

plt.plot(step, potential_energy)
plt.xlabel("Step")
plt.ylabel("Potential energy (kJ/mol)")
plt.show()
plt.plot(step, temperature)
plt.xlabel("Step")
plt.ylabel("Temperature (K)")
plt.show()
plt.plot(step, volume)
plt.xlabel("Step")
plt.ylabel("Volume (nm^3)")
plt.show()