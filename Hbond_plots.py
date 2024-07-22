import numpy as np
import matplotlib.pyplot as plt

#Set default font size
plt.rcParams.update({'font.size':20})

# Files for data inputs
file1 = 'H_rmsf.dat'
file2 = 'propka_rmsf.dat'
file3 = 'pypka_rmsf.dat'

#Convert data from file
def read_data(file):
    with open(file, 'r') as f:
        data = np.loadtxt(f)
    return data

#Read data from files
rmsf1 = read_data(file1)
rmsf2 = read_data(file2)
rmsf3 = read_data(file3)

#Assigning residue and rmsf values
residue = rmsf1[:, 0]
rmsf1_values = rmsf1[:, 1]
rmsf2_values = rmsf2[:, 1]
rmsf3_values = rmsf3[:, 1]

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(residue, rmsf1_values, label='H++', color='blue')
plt.plot(residue, rmsf2_values, label='PropKa', color='gray')
plt.plot(residue, rmsf3_values, label='PypKa', color='red')
plt.xlabel('Residue Number')
plt.ylabel('RMSF (nm)')
    
# Show plot
plt.show()