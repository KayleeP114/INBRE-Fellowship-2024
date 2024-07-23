def filter_interactions(file_path, threshold_fraction=0.50):
    interactions = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.split()
            fraction = float(parts[-1])
            
            if fraction >= threshold_fraction:
                interactions.append(fraction)
    
    return interactions

def calculate_average_bonds_from_filtered(interactions, total_snapshots):
    if not interactions:
        return 0
    
    # Calculate the average fraction from the filtered interactions
    average_fraction = sum(interactions) / len(interactions)
    # Calculate the average number of bonds across all snapshots
    average_bonds_total = average_fraction * total_snapshots
    return average_bonds_total

def main(file_path, total_snapshots):
    # Filter interactions present in 50% or more of the frames
    filtered_fractions = filter_interactions(file_path, threshold_fraction=0.50)
    
    # Calculate the average number of hydrogen bonds based on the filtered interactions
    average_bonds = calculate_average_bonds_from_filtered(filtered_fractions, total_snapshots)
    
    print(f"The average number of hydrogen bonds based on interactions present in 50% or more of the frames is: {average_bonds:.2f}")

# Example usage
file_path = 'propka_hbond.dat'  
total_snapshots = 1000  

if __name__ == "__main__":
    main(file_path, total_snapshots)
