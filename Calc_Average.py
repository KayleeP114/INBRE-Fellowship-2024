def calculate_average_number_of_hydrogen_bonds(file_path):
    # Initialize variables
    total_bonds = 0
    total_fractions = 0

    # Open and read the file
    with open(file_path, 'r') as file:
        for line in file:
            # Skip comment lines or empty lines
            if line.startswith('#') or not line.strip():
                continue
            
            # Split the line into parts
            parts = line.split()
            # The last part is the fraction value
            fraction = float(parts[-1])
            
            # Approximate number of hydrogen bonds for this entry
            # Assuming fraction represents the number of bonds observed
            total_bonds += fraction
            total_fractions += 1
    
    # Calculate average
    if total_fractions == 0:
        return 0
    average_number_of_bonds = total_bonds / total_fractions
    return average_number_of_bonds

# Example usage
file_path = 'H_hbond.dat'  # Replace with your file path
average_bonds = calculate_average_number_of_hydrogen_bonds(file_path)
print(f"The average number of hydrogen bonds is: {average_bonds:.2f}")
