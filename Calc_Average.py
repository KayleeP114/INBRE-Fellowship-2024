# average_number_of_hydrogen_bonds.py

def calculate_average_number_of_hydrogen_bonds(file_path, total_snapshots):
    # Initialize variables
    total_bonds = 0
    total_count = 0

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
            # Each fraction represents the bond count for one snapshot
            bonds_per_snapshot = fraction
            total_bonds += bonds_per_snapshot
            total_count += 1
    
    # Calculate average number of hydrogen bonds per snapshot
    if total_count == 0:
        return 0
    average_bonds_per_snapshot = total_bonds / total_count
    
    # Calculate average number of hydrogen bonds over all snapshots
    average_bonds_total = average_bonds_per_snapshot * total_snapshots
    return average_bonds_total

# Example usage
file_path = 'propka_hbond.dat'  
total_snapshots = 1000  
average_bonds = calculate_average_number_of_hydrogen_bonds(file_path, total_snapshots)
print(f"The average number is: {average_bonds:.2f}")
