def calculate_average_hydrogen_bonds(file_path):
    # Initialize variables
    total_fraction = 0
    count = 0

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
            
            # Accumulate the total fraction and count
            total_fraction += fraction
            count += 1
    
    # Calculate average
    if count == 0:
        return 0
    average_fraction = total_fraction / count
    return average_fraction


# Example usage
file_path = 'propka_saltbr.dat'  
average = calculate_average_hydrogen_bonds(file_path)
print(f"The average number is: {average:.2f}")