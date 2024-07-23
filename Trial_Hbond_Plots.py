import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

# Combine data into a single DataFrame for easier plotting
combined_data = pd.concat([
    data1.assign(Source='File 1'),
    data2.assign(Source='File 2'),
    data3.assign(Source='File 3')
])

# Plotting the data
plt.figure(figsize=(15, 10))

# Unique interactions
interactions = combined_data['Interaction'].unique()

# Create a bar plot
bar_width = 0.1
indices = np.arange(len(interactions))

# Plot bars for each file
plt.bar(indices, data1['Fraction'], bar_width, label='File 1')
plt.bar(indices + bar_width, data2['Fraction'], bar_width, label='File 2')
plt.bar(indices + 2 * bar_width, data3['Fraction'], bar_width, label='File 3')

# Customize the plot
plt.xticks(indices + bar_width, interactions, rotation=90)
plt.xlabel('Interaction')
plt.ylabel('Fraction of Frames')


# Display the plot
plt.show()
