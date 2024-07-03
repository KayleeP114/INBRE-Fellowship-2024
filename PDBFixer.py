from pdbfixer import PDBFixer

# Define paths to input and output files
input_pdb_file = 'K62_pH4_propka.pdb'
output_pdb_file = 'Fixed_K62_pH4_propka.pdb'

# Initialize PDBFixer instance
fixer = PDBFixer(filename=input_pdb_file)

# Example of actions you might perform
# For example, you could add missing hydrogens
fixer.addMissingHydrogens()

# Apply fixes
fixer.findMissingResidues()
fixer.removeHeterogens(False)

# Write out the fixed PDB file
with open(output_pdb_file, 'w') as f:
    PDBFixer.writeFile(fixer.topology, fixer.positions, f)

print(f"Fixed PDB written to {output_pdb_file}")
