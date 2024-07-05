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
u = mda.Universe('K62_pH4_PypKa.pdb', 'pypka_trajectory.dcd')

### Select all atoms
all_atoms = u.select_atoms('all')
protein = u.select_atoms('protein')
reference = protein

### Protonation state data from predictions
protonation_states = {
    ('NTR', 'A', 32): 8.27,
    ('ASP', 'A', 33): 2.72,
    ('ASP', 'A', 41): 3.78,
    ('HIS', 'A', 45): 6.73,
    ('ASP', 'A', 53): 4.29,
    ('TYR', 'A', 56): 10.91,
    ('HIS', 'A', 57): 7.27,
    ('ASP', 'A', 58): 3.56,
    ('TYR', 'A', 61): None,
    ('TYR', 'A', 62): 10.25,
    ('GLU', 'A', 63): 3.93,
    ('LYS', 'A', 69): 11.56,
    ('ASP', 'A', 72): 3.02,
    ('HIS', 'A', 75): 6.45,
    ('ASP', 'A', 86): 2.64,
    ('GLU', 'A', 88): 4.68,
    ('TYR', 'A', 93): None,
    ('TYR', 'A', 94): None,
    ('ASP', 'A', 106): 2.63,
    ('LYS', 'A', 111): 10.50,
    ('TYR', 'A', 113): 9.95,
    ('ASP', 'A', 117): 3.61,
    ('TYR', 'A', 118): 10.52,
    ('ASP', 'A', 119): 3.53,
    ('TYR', 'A', 122): 11.87,
    ('HIS', 'A', 123): 6.31,
    ('TYR', 'A', 135): 9.22,
    ('TYR', 'A', 141): 10.11,
    ('LYS', 'A', 148): 11.25,
    ('ASP', 'A', 153): 3.45,
    ('ASP', 'A', 157): 3.56,
    ('GLU', 'A', 159): 4.28,
    ('TYR', 'A', 164): None,
    ('TYR', 'A', 166): 11.48,
    ('TYR', 'A', 170): 9.96,
    ('ASP', 'A', 171): 2.88,
    ('GLU', 'A', 179): 3.64,
    ('ASP', 'A', 186): 4.11,
    ('TYR', 'A', 187): 10.71,
    ('TYR', 'A', 195): 10.55,
    ('TYR', 'A', 213): None,
    ('ASP', 'A', 221): 2.45,
    ('LYS', 'A', 224): None,
    ('TYR', 'A', 226): 11.02,
    ('LYS', 'A', 234): 10.52,
    ('ASP', 'A', 238): 2.73,
    ('TYR', 'A', 243): 10.35,
    ('TYR', 'A', 244): 10.40,
    ('ASP', 'A', 257): 2.98,
    ('ASP', 'A', 266): 3.11,
    ('CYS', 'A', 269): 10.84,
    ('CTR', 'A', 272): 2.94
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
    filtered_atoms = [atom for atom in all_atoms if atom.atom_name in backbone_atoms]
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
    plt.title('RMSD Over Time for ####')
    plt.legend()
    plt.savefig(f'####_rmsd_dark_{interval_ns}ns_{timestamp}.png')
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
    plt.title('Radius of Gyration Over Time for ####')
    plt.savefig(f'####_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png')
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
    plt.title('Hydrogen Bonds Over Time for ####')
    plt.savefig(f'####_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png')
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
    plt.title('RMSF for #####')
    plt.legend()
    plt.savefig(f'#####_rmsf_dark_{timestamp}.png')
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
    
    print(f"Analysis complete. Plots saved as 'PypKa_rmsd_dark_{interval_ns}ns_{timestamp}.png', 'PypKa_radius_of_gyration_dark_{interval_ns}ns_{timestamp}.png', 'PypKa_hydrogen_bonds_dark_{interval_ns}ns_{timestamp}.png', and 'PypKa_rmsf_dark_{timestamp}.png'.")

if __name__ == "__main__":
    interval_ns = 1  #### Set time interval here
    main(interval_ns)