import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#set default font size
plt.rcParams.update({'font.size':20})

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

# Merge data on the Interaction column
combined_data = pd.merge(data1, data2, on='Interaction', how='outer', suffixes=('_1', '_2'))
combined_data = pd.merge(combined_data, data3, on='Interaction', how='outer')
combined_data.columns = ['Interaction', 'Fraction_1', 'Fraction_2', 'Fraction_3']

# Replace NaN values with 0 for plotting
combined_data.fillna(0, inplace=True)

# Plotting the data
plt.figure(figsize=(15, 10))

# Unique interactions
interactions = combined_data['Interaction']

# Create a bar plot
bar_width = 0.5
indices = np.arange(len(interactions))

# Plot bars for each file
plt.bar(indices, combined_data['Fraction_1'], bar_width, label='H++', color='blue')
plt.bar(indices + bar_width, combined_data['Fraction_2'], bar_width, label='PropKa', color='gray')
plt.bar(indices + 2 * bar_width, combined_data['Fraction_3'], bar_width, label='PypKa', color='red')

# Customize the plot
plt.xlabel('Interaction')
plt.ylabel('Fraction')

# Display the plot
plt.show()

