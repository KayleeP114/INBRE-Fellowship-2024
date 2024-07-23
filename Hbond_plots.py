import numpy as np
import matplotlib.pyplot as plt

#Set default font size
plt.rcParams.update({'font.size':20})

# Function to read hydrogen bond data
def read_hbond_data(data_file):
    bonds = []
    fractions = []
    with open(data_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.split()
                bonds = ''.join(parts[:-1])
                fraction = float(parts[-1])
                bonds.append(bond)
                fractions.append(fraction)
    return bonds, fractions

#Read data from files
bonds1, fractions1 = read_hbond_data('H_hbond.dat')
bonds2, fractions2 = read_hbond_data('propka_hbond.dat')
bonds3, fractions3 = read_hbond_data('pypka_hbond.dat')

#Convert lists to numpy arrays
bonds = np.array(bonds)
fractions = np.array(fractions)

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(bonds1, fractions1, label='H++', color='blue')
plt.plot(bonds2, fractions2, label='PropKa', color='gray')
plt.plot(bonds3, fractions3, label='PypKa', color='red')
plt.xlabel('Interactions')
plt.ylabel('Fraction of Frames')
    
# Show plot
plt.show()