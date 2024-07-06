import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

### Set to dark mode
plt.style.use('dark_background')

### Load the PDB/PQR and DCD files
u = mda.Universe('K62_pH4_H++.pdb', 'H++_trajectory.dcd')

### Select all atoms
all_atoms = u.select_atoms('all')
protein = u.select_atoms('protein')
reference = protein

### Protonation state data from predictions
protonation_states = {
    ('NT', 'A', 32): 7.7,
    ('ASP', 'A', 33): 2.4,
    ('ASP', 'A', 41): 3.9,
    ('HIS', 'A', 45): 6.7,
    ('ASP', 'A', 53): 3.9,
    ('TYR', 'A', 56): 11.8,
    ('HIS', 'A', 57): 8.2,
    ('ASP', 'A', 58): 3.5,
    ('TYR', 'A', 61): None,
    ('TYR', 'A', 62): 10.6,
    ('GLU', 'A', 63): 3.7,
    ('LYS', 'A', 69): 11.5,
    ('ASP', 'A', 72): 3.2,
    ('HIS', 'A', 75): 6.8,
    ('ASP', 'A', 86): 2.0,
    ('GLU', 'A', 88): 5.1,
    ('TYR', 'A', 93): None,
    ('TYR', 'A', 94): None,
    ('ASP', 'A', 106): 2.4,
    ('LYS', 'A', 111): 10.2,
    ('ARG', 'A', 112): None,
    ('TYR', 'A', 113): 10.5,
    ('ASP', 'A', 117): 4.5,
    ('TYR', 'A', 118): 10.7,
    ('ASP', 'A', 119): 3.3,
    ('TYR', 'A', 122): None,
    ('HIS', 'A', 123): 6.1,
    ('TYR', 'A', 135): 9.1,
    ('TYR', 'A', 141): 10.8,
    ('LYS', 'A', 148): 11.0,
    ('ASP', 'A', 153): 3.8,
    ('ASP', 'A', 157): 3.6,
    ('GLU', 'A', 159): 4.3,
    ('TYR', 'A', 164): None,
    ('TYR', 'A', 166): None,
    ('TYR', 'A', 170): 10.3,
    ('ASP', 'A', 171): 2.5,
    ('GLU', 'A', 179): 3.6,
    ('ASP', 'A', 186): 4.6,
    ('TYR', 'A', 187): 11.4,
    ('TYR', 'A', 195): 11.0,
    ('TYR', 'A', 213): None,
    ('ARG', 'A', 216): None,
    ('ASP', 'A', 221): 2.7,
    ('LYS', 'A', 224): None,
    ('TYR', 'A', 226): 11.5,
    ('ARG', 'A', 228): None,
    ('LYS', 'A', 234): 10.5,
    ('ASP', 'A', 238): 3.0,
    ('TYR', 'A', 243): 10.9,
    ('TYR', 'A', 244): 10.9,
    ('ARG', 'A', 245): None,
    ('ASP', 'A', 257): 2.6,
    ('ARG', 'A', 258): None,
    ('ASP', 'A', 266): 3.8,
    ('CYS', 'A', 269): None,
    ('ARG', 'A', 270): None,
    ('CT', 'A', 272): 3.9
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
def calculate_rmsf():
    rmsf = RMSF(protein.select_atoms('name CA')).run()
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
    plt.title('RMSD Over Time for H++')
    plt.legend()
    plt.savefig(f'H++_rmsd_dark_{interval_ns}ns_{timestamp}.png')
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
    plt.title('Radius of Gyration Over Time for H++')
    plt.savefig(f'H++_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png')
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
    plt.title('Hydrogen Bonds Over Time for H++')
    plt.savefig(f'H++_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png')
    plt.close()
    print("Hydrogen Bonds plot saved.")

# Plot RMSF
def plot_rmsf(rmsf):
    print("Plotting RMSF...")
    residues = np.arange(len(rmsf))
    plt.figure()
    plt.plot(residues, rmsf, color='cyan', label='RMSF')
    avg_rmsf = np.mean(rmsf)
    plt.axhline(y=avg_rmsf, color='yellow', linestyle='--', label=f'Avg RMSF: {avg_rmsf:.2f} Å')
    step = max(1, len(residues) // 10)  
    for i in range(0, len(residues), step):
        plt.text(residues[i], rmsf[i], f'{rmsf[i]:.2f}', fontsize=8, ha='center')
    plt.xlabel('Residue')
    plt.ylabel('RMSF (Å)')
    plt.title('RMSF for H++')
    plt.legend()
    plt.savefig(f'H++_rmsf_dark_{timestamp}.png')
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
    rmsf = calculate_rmsf()
    plot_rmsf(rmsf)
    
    print(f"Analysis complete. Plots saved as 'H++_rmsd_dark_{interval_ns}ns_{timestamp}.png', 'H++_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png', 'H++_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png', and 'H++_rmsf_dark_{timestamp}.png'.")

if __name__ == "__main__":
    interval_ns = 1  #### Set time interval here
    main(interval_ns)