import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

### Set to dark mode
plt.style.use('dark_background')

### Load the PDB/PQR and DCD files from H++
u = mda.Universe('K62_pH4_propka.pqr', 'propka_trajectory.dcd')

### Select all atoms
all_atoms = u.select_atoms('all')
protein = u.select_atoms('protein')
reference = protein

### Protonation state data from predictions
protonation_states = {
    ('N+', 'A', 32): 7.91,
    ('ASP', 'A', 33): 2.16,
    ('ASP', 'A', 41): 3.92,
    ('HIS', 'A', 45): 6.41,
    ('ASP', 'A', 53): 4.03,
    ('TYR', 'A', 56): 10.57,
    ('HIS', 'A', 57): 6.91,
    ('ASP', 'A', 58): 3.87,
    ('TYR', 'A', 61): 11.15,
    ('TYR', 'A', 62): 10.30,
    ('GLU', 'A', 63): 4.25,
    ('LYS', 'A', 69): 10.64,
    ('ASP', 'A', 72): 3.35,
    ('HIS', 'A', 75): 6.27,
    ('ASP', 'A', 86): 3.45,
    ('GLU', 'A', 88): 4.83,
    ('TYR', 'A', 93): 12.45,
    ('TYR', 'A', 94): 12.82,
    ('ASP', 'A', 106): 2.45,
    ('LYS', 'A', 111): 10.45,
    ('ARG', 'A', 112): 12.68,
    ('TYR', 'A', 113): 10.22,
    ('ASP', 'A', 117): 4.23,
    ('TYR', 'A', 118): 10.33,
    ('ASP', 'A', 119): 3.89,
    ('TYR', 'A', 122): 10.83,
    ('HIS', 'A', 123): 5.34,
    ('TYR', 'A', 135): 9.70,
    ('TYR', 'A', 141): 10.59,
    ('CYS', 'A', 146): 99.99,
    ('LYS', 'A', 148): 10.59,
    ('ASP', 'A', 153): 3.52,
    ('ASP', 'A', 157): 3.23,
    ('GLU', 'A', 159): 4.70,
    ('TYR', 'A', 164): 11.19,
    ('TYR', 'A', 166): 11.80,
    ('TYR', 'A', 170): 10.26,
    ('ASP', 'A', 171): 2.63,
    ('GLU', 'A', 179): 3.60,
    ('ASP', 'A', 186): 4.20,
    ('TYR', 'A', 187): 10.49,
    ('TYR', 'A', 195): 10.27,
    ('CYS', 'A', 201): 99.99,
    ('TYR', 'A', 213): 13.03,
    ('ARG', 'A', 216): 12.70,
    ('ASP', 'A', 221): 2.95,
    ('LYS', 'A', 224): 9.89,
    ('TYR', 'A', 226): 10.61,
    ('CYS', 'A', 227): 99.99,
    ('ARG', 'A', 228): 12.84,
    ('LYS', 'A', 234): 10.43,
    ('CYS', 'A', 237): 99.99,
    ('ASP', 'A', 238): 3.34,
    ('TYR', 'A', 243): 10.34,
    ('TYR', 'A', 244): 9.77,
    ('ARG', 'A', 245): 13.63,
    ('ASP', 'A', 257): 2.89,
    ('ARG', 'A', 258): 12.96,
    ('CYS', 'A', 262): 99.99,
    ('ASP', 'A', 266): 3.68,
    ('CYS', 'A', 269): 10.09,
    ('ARG', 'A', 270): 12.32,
    ('CYS', 'A', 271): 99.99,
    ('C-', 'A', 272): 3.32
}

### RMSD Analysis
def calculate_rmsd(interval_ns=1):
    ### Align the trajectory to the reference structure
    align.AlignTraj(u, reference, select='protein and name CA', in_memory=True).run()
    R = rms.RMSD(u, reference, select='protein and name CA')
    R.run()
    return R.results.rmsd[::interval_ns]  

### Radius of Gyration Analysis
def calculate_radius_of_gyration(interval_ns=1):
    Rg = []
    for i, ts in enumerate(u.trajectory):
        if i % interval_ns == 0:
            Rg.append(all_atoms.radius_of_gyration())
    return np.array(Rg)

### Hydrogen Bond Analysis
def calculate_hydrogen_bonds(interval_ns=1):
    donors_sel = 'name N HN'
    acceptors_sel = 'name O'
    hydrogens_sel = 'name H'
    h = HBA(universe=u, donors_sel=donors_sel, acceptors_sel=acceptors_sel, hydrogens_sel=hydrogens_sel)
    h.run()
    return h.count_by_time()[::interval_ns]

#### RMSF Analysis
def calculate_rmsf(all_atoms):
    backbone_atoms = ['N', 'CA', 'C', 'O']
    filtered_atoms = [atom for atom in all_atoms if atom.name in backbone_atoms]
    rmsf = RMSF(filtered_atoms).run()
    return rmsf.rmsf

#### Generate a timestamp
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

#### Plot RMSD
def plot_rmsd(rmsd, interval_ns):
    print("Plotting RMSD...")
    time = np.arange(0, len(rmsd) * interval_ns, interval_ns)  # Adjust the time axis based on interval
    plt.figure()
    plt.plot(time, rmsd[:, 2], color='cyan')
    avg_rmsd = np.mean(rmsd[:, 2])
    plt.axhline(y=avg_rmsd, color='yellow', linestyle='--', label=f'Avg RMSD: {avg_rmsd:.2f} Å')
    step = max(1, len(time) // 10)  
    for i in range(0, len(time), step):  
        plt.text(time[i], rmsd[i, 2], f'{rmsd[i, 2]:.2f}', fontsize=8, ha='center')
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.title('RMSD Over Time for propka')
    plt.legend()
    plt.savefig(f'propka_rmsd_dark_{interval_ns}ns_{timestamp}.png')
    plt.close()
    print("RMSD plot saved.")

#### Plot Radius of Gyration
def plot_radius_of_gyration(Rg, interval_ns):
    print("Plotting Radius of Gyration...")
    time = np.arange(0, len(Rg) * interval_ns, interval_ns)
    plt.figure()
    plt.plot(time, Rg, color='cyan')
    step = max(1, len(time) // 10)  
    for i in range(0, len(time), step):  
        plt.text(time[i], Rg[i], f'{Rg[i]:.2f}', fontsize=8, ha='center')
    plt.xlabel('Time (ns)')
    plt.ylabel('Radius of Gyration (Å)')
    plt.title('Radius of Gyration Over Time for propka')
    plt.savefig(f'propka_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png')
    plt.close()
    print("Radius of Gyration plot saved.")

### Plot Hydrogen Bonds
def plot_hydrogen_bonds(hbond_counts, interval_ns):
    print("Plotting Hydrogen Bonds...")
    time = np.arange(0, len(hbond_counts) * interval_ns, interval_ns)
    plt.figure()
    plt.plot(time, hbond_counts, color='cyan')
    step = max(1, len(time) // 10)  
    for i in range(0, len(time), step):  
        plt.text(time[i], hbond_counts[i], f'{hbond_counts[i]}', fontsize=8, ha='center')
    plt.xlabel('Time (ns)')
    plt.ylabel('Number of Hydrogen Bonds')
    plt.title('Hydrogen Bonds Over Time for propka')
    plt.savefig(f'propka_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png')
    plt.close()
    print("Hydrogen Bonds plot saved.")

# Plot RMSF
def plot_rmsf(rmsf):
    print("Plotting RMSF...")
    residues = np.arange(len(rmsf))
    plt.figure()
    plt.plot(residues, rmsf, color='cyan')
    avg_rmsf = np.mean(rmsf)
    plt.axhline(y=avg_rmsf, color='yellow', linestyle='--', label=f'Avg RMSF: {avg_rmsf:.2f} Å')
    step = max(1, len(residues) // 10)  
    for i in range(0, len(residues), step):
        plt.text(residues[i], rmsf[i], f'{rmsf[i]:.2f}', fontsize=8, ha='center')
    plt.xlabel('Residue')
    plt.ylabel('RMSF (Å)')
    plt.title('RMSF for #propka')
    plt.legend()
    plt.savefig(f'propka_rmsf_dark_{timestamp}.png')
    plt.close()
    print("RMSF plot saved.")

#### Highlight Key Residues Based on pKa
def highlight_key_residues(protonation_states, pH=4.0):
    key_residues = {k: v for k, v in protonation_states.items() if v and abs(v - pH) < 1.0}
    print(f"Key residues near pH {pH}:")
    for res, pka in key_residues.items():
        print(f"Residue: {res}, pKa: {pka}")
    return key_residues

### Main function
def main(interval_ns=1):
    #### Print the number of frames
    num_frames = u.trajectory.n_frames
    print(f"Number of frames in the simulation: {num_frames}")

    #### Highlight key residues
    key_residues = highlight_key_residues(protonation_states)
    
    #### Calculate and plot RMSD
    rmsd = calculate_rmsd(interval_ns)
    plot_rmsd(rmsd, interval_ns)
    
    ### Calculate and plot Radius of Gyration
    Rg = calculate_radius_of_gyration(interval_ns)
    plot_radius_of_gyration(Rg, interval_ns)
    
    #### Calculate and plot Hydrogen Bonds
    hbond_counts = calculate_hydrogen_bonds(interval_ns)
    plot_hydrogen_bonds(hbond_counts, interval_ns)
    
    #### Calculate and plot RMSF
    rmsf = calculate_rmsf(all_atoms)
    plot_rmsf(rmsf)
    
    print(f"Analysis complete. Plots saved as 'propka_rmsd_dark_{interval_ns}ns_{timestamp}.png', 'propka_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png', 'propka_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png', and 'propka_rmsf_dark_{timestamp}.png'.")

if __name__ == "__main__":
    interval_ns = 1  #### Set time interval here
    main(interval_ns)