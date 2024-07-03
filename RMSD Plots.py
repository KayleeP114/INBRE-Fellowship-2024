import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

# Load the trajectory and the reference PDB file
u = mda.Universe("your_structure.pdb", "your_trajectory.dcd")

# Select the atoms for RMSD calculation (e.g., CA atoms of protein)
ref = mda.Universe("your_structure.pdb")
atom_selection = u.select_atoms("protein and name CA")

# Calculate RMSD
R = rms.RMSD(atom_selection, ref, select="protein and name CA")
R.run()

# Plot RMSD
time = R.rmsd[:, 0]  # Extract time or frame numbers
rmsd_values = R.rmsd[:, 2]  # Extract RMSD values

plt.figure(figsize=(10, 6))
plt.plot(time, rmsd_values, label="CA RMSD (Å)")
plt.xlabel("Frame")
plt.ylabel("RMSD (Å)")
plt.title("RMSD over Time")
plt.legend()
plt.grid(True)
plt.show()
