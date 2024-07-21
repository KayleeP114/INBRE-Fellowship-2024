import MDAnalysis as mda
from MDAnalysis.analysis.align import *
from Bio.PDB.MMCIFParser import MMCIFParser

# Define paths to your trajectory and CIF file
trajectory_file = 'H++_trajectory.dcd'
cif_file = 'H++_output.cif'

# Initialize MDAnalysis Universe
u = mda.Universe(cif_file, trajectory_file)

# Select chain C (assuming chainid 0 corresponds to chain A, etc.)
protein_chain = u.select_atoms('chainid 0')

# Initialize a BioPython structure from the CIF file
parser = MMCIFParser()
structure = parser.get_structure('protein', cif_file)

# Extract coordinates from the CIF for chain C
ref_atoms = []
for chain in structure[0]:
    if chain.id == 'C':
        for residue in chain:
            for atom in residue:
                ref_atoms.append(atom.coord)

# Perform alignment and superposition
alignment = superpose.Alibe(u.atoms.positions, ref_atoms)
alignment.run()

# Output aligned coordinates as a PDB trajectory
with mda.Writer('H++_aligned_trajectory.pdb', protein_chain.n_atoms) as pdb_writer:
    for ts in u.trajectory:
        protein_chain.positions = alignment.transform(protein_chain.positions)
        protein_chain.write(pdb_writer)

print("Alignment and output complete.")
