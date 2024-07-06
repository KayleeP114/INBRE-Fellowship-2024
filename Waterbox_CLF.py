from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller, ForceField
from openmm import unit

# Step 1: Read in the PDB file
fixer = PDBFixer(filename='k1_H++.pdb')

# Step 2: Add missing atoms and residues (excluding hydrogens)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Step 3: Prepare the force field and modeller
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(fixer.topology, fixer.positions)

# Step 4: Add a dodecahedral water box with padding of 10-12 Å
padding = 1.2 * unit.nanometers  # 12 Å
modeller.addSolvent(forcefield, model='tip3p', padding=padding, boxShape='dodecahedron')

# Step 5: Add ions to achieve 0.15 M NaCl concentration
modeller.addSolvent(forcefield, model='tip3p', ionicStrength=0.15 * unit.molar)

# Save the system as a CIF file
with open('output.cif', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)