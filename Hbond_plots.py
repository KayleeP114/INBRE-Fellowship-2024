import numpy as np
import matplotlib.pyplot as plt

#Set default font size
plt.rcParams.update({'font.size':20})

# Function to plot and save hbond data
def plot_hbonds(data_files, output_file):
    all_bonds = []
    all_fractions = []
    labels = ['H++', 'PropKa', 'PypKa']

    for i, data_file in enumerate(data_files):
        bonds = []
        fractions = []
        with open(data_file, 'r') as file:
            for line in file:
                if not line.startswith('#'):
                    parts = line.split()
                    bond = ''.join(parts[:-1])
                    fraction = float(parts[-1])
                    bonds.append(bond)
                    fractions.append(fraction)
        all_bonds.append(bonds)
        all_fractions.append(fractions)

#Convert lists to numpy arrays
max_bonds_length = max(len(bonds) for bonds in all_bonds)
bar_width = 0.2

# Plot the data
plt.figure(figsize=(12, 8))

for i, (bonds, fractions) in enumerate(zip(all_bonds, all_fractions)):
    y_positions = np.arange(len(bonds))
    plt.barh(y_positions + i * bar_width, fractions, height=bar_width, label=label[i], align='center')

#Customize Plot
plt.yticks(y_positions + bar_width, all_bonds[0])
plt.xlabel('Interactions')
plt.ylabel('Fraction of Frames')
    
# Show plot
plt.show()

#Plot hbonds
plot_hbonds(['H_hbond.dat', 'propka_hbond.dat', 'pypka_hbond.dat'])