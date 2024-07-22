import numpy as np
import matplotlib.pyplot as plt

#Set default font size
plt.rcParams.update({'font.size':20})

# Files for data inputs
file1 = 'H_hbond.dat'
file2 = 'propka_hbond.dat'
file3 = 'pypka_hbond.dat'

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(bonds, fractions, label='H++', color='blue')
plt.plot(bonds, fractions, label='PropKa', color='gray')
plt.plot(bonds, fractions, label='PypKa', color='red')
plt.xlabel('Interactions')
plt.ylabel('Fraction of Frames')
    
# Show plot
plt.show()