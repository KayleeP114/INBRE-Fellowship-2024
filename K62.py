from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

# Create a system using the force field (use either amber14 or charmm36)
pdb = PDBFile('K62.pdb')  # Load a PDB file
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

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