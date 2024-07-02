import os
from simtk.openmm.app import PDBFile, Simulation
from simtk.openmm import LangevinIntegrator, Platform

# Define input and output file paths
input_pdb = "K62.pdb"
output_pdb = "k62_no_H.pdb"

# Load the PDB file
pdb = PDBFile(input_pdb)

# Create a system (you'll need to define your force field here)
# For example, using Amber99SBildn force field:
from openmm.app import Amber99SBildn
forcefield = Amber99SBildn()
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)

# Set up the simulation
integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
platform = Platform.getPlatformByName("CPU")  # Use "CUDA" for GPU acceleration
simulation = Simulation(pdb.topology, system, integrator, platform)

# Set initial positions
simulation.context.setPositions(pdb.positions)

# Minimize energy
simulation.minimizeEnergy()

# Run the simulation (you can adjust the number of steps)
simulation.step(1000)

# Save the output PDB
simulation.saveState(output_pdb)
print(f"Hydrogens removed. Output saved to {output_pdb}")
