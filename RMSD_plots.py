import numpy as np
import matplotlib.pyplot as plt

#Set default font size
plt.rcParams.update({'font.size':20})

# Files for data inputs
file1 = 'H_rmsd.dat'
file2 = 'propka_rmsd.dat'
file3 = 'pypka_rmsd.dat'

#Convert data from file
def read_data(file):
    with open(file, 'r') as f:
        data = np.loadtxt(f)
    return data

#Read data from files
rmsd1 = read_data(file1)
rmsd2 = read_data(file2)
rmsd3 = read_data(file3)

#Assigning time and rmsd values
time = (0, 100)
rmsd1_values = rmsd1[:, 1]
rmsd2_values = rmsd2[:, 1]
rmsd3_values = rmsd3[:, 1]

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, rmsd1_values, label='H++', color='blue')
plt.plot(time, rmsd2_values, label='PropKa', color='gray')
plt.plot(time, rmsd3_values, label='PypKa', color='red')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
    
# Show plot
plt.show()