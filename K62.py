import parmed as pmd
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

# Load the CIF file with parmed ()
structure = pmd.load_file('')

# Convert parmed to Openmm toopology and positions
topology = structure.topology
positions = structure.positions
    
# Load a force field and prepare system
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)

# Set up integrator (with reduced timestep)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
  
# Set platform (use CUDA) (single, double, or mixed for properties choices)
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}

# Define simulation
simulation = Simulation(topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)

# Minimize the energy
print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=10000)

# Output reporters
num_snapshots = 1000
interval = int((100 * 1000000) / (num_snapshots * 2))
#100ns/1000snaps=0.1ns (note: 1ns=500000steps)

# Define steps (for 100ns sim w/ 2fs timestep)
num_steps = int(100 * 1000000 / 2)
#100ns*1000000ns/microsecond/2fs

# DCD Reporter
simulation.reporters.append(DCDReporter('trajectory.dcd', interval))

# Data report statements for energy, temperature, progress
simulation.reporters.append(StateDataReporter(sys.stdout, interval, step=True, time=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=num_steps, separator='\t'))

# Checkpoint reporter
simulation.reporters.append(CheckpointReporter('checkpoint.chk', interval))

#Run sim (w/ try-except to catch any NaN errors)
try:
    simulation.step(num_steps)
except Exception as e:
    print(f"Simulation failed with error: {e}")
    print("Check the initial configuration, constraints, and force field parameters.")

#Write final positions in CIF file
PDBxFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), open('output.cif', 'w'))