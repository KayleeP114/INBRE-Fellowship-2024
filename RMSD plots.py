import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy as np
import matplotlib.pyplot as plt

# Load the reference and target PDB files
ref = mda.Universe("K62.pdb")
target = mda.Universe("K62_pH4_H++.pdb")

# Select the atoms for RMSD calculation
ref_atoms = ref.select_atoms("protein and name CA")
target_atoms = target.select_atoms("protein and name CA")

# Ensure the atoms are aligned for RMSD calculation
R = rms.RMSD(target_atoms, ref_atoms, select="protein and name CA")
R.run()

# Extract RMSD values
rmsd = R.rmsd[0, 2]  # Since there is only one frame, get the RMSD value

# Plot RMSD (though this plot will just show a single point)
plt.figure(figsize=(6, 4))
plt.plot([0], [rmsd], 'o', label="CA RMSD (Å)")
plt.xlabel("Frame")
plt.ylabel("RMSD (Å)")
plt.title("RMSD between Reference and Target Structure")
plt.legend()
plt.show()

print(f"RMSD: {rmsd} Å")