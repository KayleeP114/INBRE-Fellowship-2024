import numpy as np
import matplotlib.pyplot as plt

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
time = rmsd1[:, 0]
rmsd1_values = rmsd1[:, 1]
rmsd2_values = rmsd2[:, 1]
rmsd3_values = rmsd3[:, 1]

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, rmsd1_values, label='H++ RMSD', color='blue')
plt.plot(time, rmsd2_values, label='PropKa RMSD', color='gray')
plt.plot(time, rmsd3_values, label='PypKa RMSD', color='red')
plt.title('RMSD', fontsize=40, fontname='Times New Roman')
plt.xlabel('Time (ns)', fontsize=40, fontname='Times New Roman')
plt.ylabel('RMSD (nm)', fontsize=40, fontname='Times New Roman')
    
# Show plot
plt.show()