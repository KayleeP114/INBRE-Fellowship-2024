from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import numpy as np
import matplotlib.pyplot as plt


# Load the PDB file
pdb = PDBFile('1AKI.pdb')


# Specify the forcefield
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')


# Clean up: remove water and add missing hydrogens
modeller = Modeller(pdb.topology, pdb.positions)
modeller.deleteWater()
modeller.addHydrogens(forcefield)


# Solvate the protein with water and ions
modeller.addSolvent(forcefield, padding=1.0*nanometer)


# Create the system
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)


# Add a barostat to control pressure
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))


# Set up the integrator
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)


# Create the simulation object
simulation = Simulation(modeller.topology, system, integrator)


# Set the initial positions
simulation.context.setPositions(modeller.positions)


# Minimize the energy
print("Minimizing energy")
simulation.minimizeEnergy()


# Set up reporters to output data
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
       potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(StateDataReporter("md_log.txt", 100, step=True,
       potentialEnergy=True, temperature=True, volume=True))


# NVT equilibration
print("Running NVT")
simulation.step(10000)


# Reinitialize with barostat for NPT production MD
simulation.context.reinitialize(preserveState=True)


print("Running NPT")
simulation.step(10000)


# Basic analysis
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