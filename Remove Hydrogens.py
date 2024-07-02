import os

input_file = "K62.pdb"
output_file = "K62 without H.pdb"

def remove_hydrogens(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                if not atom_name.startswith("H"):
                    outfile.write(line)

def parse_pdb(file_path):
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append(line)
    return atoms

def verify_removal(input_atoms, output_atoms):
    removed_atoms = [line for line in input_atoms if line not in output_atoms]
    incorrectly_removed = [line for line in removed_atoms if not line[12:16].strip().startswith('H')]
    return len(input_atoms), len(output_atoms), len(removed_atoms), len(incorrectly_removed)

if __name__ == "__main__":
    input_path = os.path.join(os.getenv("HOME"), "INBRE-Fellowship-2024",  input_file)
    output_path = os.path.join(os.getenv("HOME"), "INBRE-Fellowship-2024",  output_file)
    
    # Remove hydrogens
    remove_hydrogens(input_path, output_path)
    
    # Parse both input and output PDB files
    input_atoms = parse_pdb(input_path)
    output_atoms = parse_pdb(output_path)
    
    # Verify removal
    input_count, output_count, removed_count, incorrectly_removed_count = verify_removal(input_atoms, output_atoms)
    
    print(f"Hydrogens removed. Output saved to {output_path}")
    print(f"Total atoms in original file: {input_count}")
    print(f"Total atoms in output file: {output_count}")
    print(f"Number of atoms removed: {removed_count}")
    print(f"Number of non-hydrogen atoms incorrectly removed: {incorrectly_removed_count}")
    
    if incorrectly_removed_count == 0:
        print("Verification successful: No non-hydrogen atoms were incorrectly removed.")
    else:
        print("Verification failed: Some non-hydrogen atoms were incorrectly removed.")