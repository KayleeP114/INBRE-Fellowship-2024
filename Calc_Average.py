def calculate_average_number_of_hydrogen_bonds(file_path, total_snapshots):
    total_bonds = 0
    total_count = 0

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.split()
            fraction = float(parts[-1])
            
            bonds_per_snapshot = fraction
            total_bonds += bonds_per_snapshot
            total_count += 1
    
    if total_count == 0:
        return 0
    average_bonds_per_snapshot = total_bonds / total_count
    average_bonds_total = average_bonds_per_snapshot * total_snapshots
    return average_bonds_total

def filter_interactions(file_path, threshold_fraction=0.50):
    interactions = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.split()
            fraction = float(parts[-1])
            
            if fraction >= threshold_fraction:
                interactions.append(line.strip())
    
    return interactions

def main(file_path, total_snapshots):
    # Calculate the average number of hydrogen bonds
    average_bonds = calculate_average_number_of_hydrogen_bonds(file_path, total_snapshots)
    print(f"The average number of hydrogen bonds is: {average_bonds:.2f}")
    
    # Filter interactions present in 50% or more of the frames
    filtered_interactions = filter_interactions(file_path, threshold_fraction=0.50)
    print("\nInteractions present in 50% or more of the frames:")
    for interaction in filtered_interactions:
        print(interaction)

# Example usage
file_path = 'H_hbond.dat'  
total_snapshots = 1000  

if __name__ == "__main__":
    main(file_path, total_snapshots)

