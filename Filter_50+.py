def filter_interactions(file_path, threshold_fraction=0.50):
    interactions = []

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
            
            # Check if the fraction meets the threshold
            if fraction >= threshold_fraction:
                interactions.append(line.strip())
    
    return interactions

# Example usage
file_path = 'pypka_saltbr.dat'  
filtered_interactions = filter_interactions(file_path)

print("Interactions present in 50% or more of the frames:")
for interaction in filtered_interactions:
    print(interaction)