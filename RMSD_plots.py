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

#Convert frames to time in ns
frames_to_time = 0.1
time = rmsd1[:, 0] * frames_to_time

#Assigning rmsd values
rmsd1_values = rmsd1[:, 1]
rmsd2_values = rmsd2[:, 1]
rmsd3_values = rmsd3[:, 1]

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, rmsd1_values, label='H++', color='purple')
plt.plot(time, rmsd2_values, label='PropKa', color='gray')
plt.plot(time, rmsd3_values, label='PypKa', color='green')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
    
# Show plot
plt.show()

#Calculate Average 
average_rmsd1 = np.mean(rmsd1_values)
average_rmsd2 = np.mean(rmsd2_values)
average_rmsd3 = np.mean(rmsd3_values)

print(f'Average RMSD (1): {average_rmsd1:.3f} nm')
print(f'Average RMSD (2): {average_rmsd2:.3f} nm')
print(f'Average RMSD (3): {average_rmsd3:.3f} nm')