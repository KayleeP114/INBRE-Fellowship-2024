import pandas as pd
import matplotlib.pyplot as plt

# Function to plot RMSD data
def plot_rmsd(H_rmsd.dat, title="RMSD", x_label="Time (ns)", y_label="RMSD (nm)", 
              line_color='blue', font_name='Times New Roman', font_size=20):
    # Read data from file
    data = pd.read_csv(file_path, delim_whitespace=True)
    
    # Assuming the file has two columns: time and rmsd
    time = data.iloc[:, 0]
    rmsd = data.iloc[:, 1]
    
    # Set the font properties
    plt.rc('font', family=font_name, size=font_size)
    
    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.plot(time, rmsd, color=line_color)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    # Show plot
    plt.show()