import pandas as pd
import matplotlib.pyplot as plt

# Function to read data from a file
def read_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Parse the data, skipping the header
    data = []
    for line in lines[1:]:  # Skip the header line
        parts = line.split()
        # Select necessary columns: donor, acceptor, and fraction
        donor = parts[0]
        acceptor = parts[2]
        fraction = float(parts[3])
        data.append([f"{donor} --> {acceptor}", fraction])

    return pd.DataFrame(data, columns=['Interaction', 'Fraction'])

# Read data from each file
data1 = read_data('H_hbond.dat')
data2 = read_data('propka_hbond.dat')
data3 = read_data('pypka_hbond.dat')

# Plotting the data
plt.figure(figsize=(12, 8))

# Plot each dataset
plt.plot(data1['Interaction'], data1['Fraction'], marker='o', label='H++')
plt.plot(data2['Interaction'], data2['Fraction'], marker='x', label='PropKa')
plt.plot(data3['Interaction'], data3['Fraction'], marker='^', label='PypKa')

# Customize the plot
plt.xlabel('Interaction')
plt.ylabel('Fraction')

# Display the plot
plt.show()
