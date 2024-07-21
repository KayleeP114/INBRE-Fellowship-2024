import pandas as pd
import matplotlib.pyplot as plt

# Files for data inputs
rmsd1 = 'H_rmsd.dat'
rmsd2 = 'propka_rmsd.dat'
rmsd3 = 'pypka_rmsd.dat'

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