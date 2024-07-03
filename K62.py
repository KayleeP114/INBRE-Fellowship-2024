from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

# Load PDB file
pdb = PDBFile('K62.pdb')  # Ensure this file exists and has water molecules

# Create a force field object
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Define the periodic box dimensions (assuming a cubic box with a side length of 3 nm)
box_size = 3.0 * nanometers
pdb.topology.setPeriodicBoxVectors(Vec3(10, 0, 0), Vec3(0, 10, 0))

# Create a system object from the topology
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# Create an integrator object
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Create a simulation object
simulation = Simulation(pdb.topology, system, integrator)

# Set the initial positions of the atoms
simulation.context.setPositions(pdb.positions)

# Minimize the energy
print('Minimizing energy...')
simulation.minimizeEnergy()

# Set up reporters to record simulation data
simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(DCDReporter('trajectory.dcd', 1000))

# Run the simulation
print('Running simulation...')
simulation.step(10000)  # Run for 10,000 steps

print('Simulation complete')