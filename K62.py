import parmed as pmd
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
# Load the CIF file with parmed
structure = pmd.load_file('')

# Check if periodic box dimensions are defined, if not set default dimensions
if pdb.topology.getPeriodicBoxVectors() is None:
    print("Periodic box dimensions not found. Setting default dimensions.")
    box_size = 10.0 * unit.nanometers  #Box size, adjust as needed
    vectors = (
        [box_size, 0.0 * unit.nanometers, 0.0 * unit.nanometers],
        [0.0 * unit.nanometers, box_size, 0.0 * unit.nanometers],
        [0.0 * unit.nanometers, 0.0 * unit.nanometers, box_size]
    )
    pdb.topology.setPeriodicBoxVectors(vectors)
# Load a force field
forcefield = app.ForceField('amber14-all.xml', 'tip3p.xml')
# Create the system
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds
)
# Set up the integrator
integrator = openmm.LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picoseconds,
    2.0 * unit.femtoseconds
)
# Create a simulation object
simulation = app.Simulation(pdb.topology, system, integrator)
# Set the initial positions
simulation.context.setPositions(pdb.positions)
# Minimize the energy
print("Minimizing energy...")
simulation.minimizeEnergy()
# Set up reporters to record data during the simulation
simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(app.DCDReporter('propka_trajectory.dcd', 1000))
simulation.reporters.append(app.CheckpointReporter('propka_checkpoint.chk', 1000))
# Run the simulation
print("Running simulation...")
simulation.step(10000)
print("Simulation complete.")