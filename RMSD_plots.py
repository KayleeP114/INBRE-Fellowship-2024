import pandas as pd
import matplotlib.pyplot as plt

# Files for data inputs
file1 = 'H_rmsd.dat'
file2 = 'propka_rmsd.dat'
file3 = 'pypka_rmsd.dat'

#Reading data from files
with open(file1, 'r') as f1:
    rmsd1 = f1.readlines()

with open(file2, 'r') as f2:
    rmsd2 = f2.readlines()

with open(file3, 'r') as f3:
    rmsd3 = f3.readlines()

#Setting time
time = rmsd1.iloc[:, 0]

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, rmsd1, label='RMSD 1', color='blue')
plt.plot(time, rmsd2, label='RMSD 2', color='gray')
plt.plot(time, rmsd3, label='RMSD 3', color='red')
plt.title('RMSD', fontsize=25, fontname='Times New Roman')
plt.xlabel('Time (ns)', fontsize=20, fontname='Times New Roman')
plt.ylabel('RMSD (nm)', fontsize=20, fontname='Times New Roman')
    
# Show plot
plt.show()