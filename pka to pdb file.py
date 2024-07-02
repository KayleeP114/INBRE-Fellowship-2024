import re

def parse_propka_output(propka_output_file):
    pka_values = {}
    with open(propka_output_file, 'r') as file:
        for line in file:
            match = re.match(r'(\d+)\s+([A-Z]{3})\s+(\w)\s+([0-9.]+)', line)
            if match:
                residue_number = int(match.group(1))
                residue_name = match.group(2)
                chain_id = match.group(3)
                pka_value = float(match.group(4))
                pka_values[(chain_id, residue_number, residue_name)] = pka_value
    return pka_values

def update_pdb_protonation_states(original_pdb_file, pka_values, target_pH=7.0):
    modified_pdb_lines = []
    with open(original_pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                residue_number = int(line[22:26].strip())
                residue_name = line[17:20].strip()
                key = (chain_id, residue_number, residue_name)
                if key in pka_values:
                    pka_value = pka_values[key]
                    if residue_name in ['ASP', 'GLU']:
                        if pka_value > target_pH:
                            line = line[:17] + 'ASH' + line[20:]  # Protonated form of ASP/GLU
                        else:
                            line = line[:17] + 'ASP' + line[20:]  # Deprotonated form of ASP/GLU
                    elif residue_name in ['HIS']:
                        if pka_value > target_pH:
                            line = line[:17] + 'HIP' + line[20:]  # Protonated form of HIS
                        else:
                            line = line[:17] + 'HIS' + line[20:]  # Deprotonated form of HIS
                    elif residue_name in ['LYS']:
                        if pka_value > target_pH:
                            line = line[:17] + 'LYS' + line[20:]  # Protonated form of LYS
                        else:
                            line = line[:17] + 'LYN' + line[20:]  # Deprotonated form of LYS
                    elif residue_name in ['ARG']:
                        if pka_value > target_pH:
                            line = line[:17] + 'ARG' + line[20:]  # Protonated form of ARG
                        else:
                            line = line[:17] + 'ARN' + line[20:]  # Deprotonated form of ARG
            modified_pdb_lines.append(line)
    return modified_pdb_lines

def save_modified_pdb(modified_pdb_lines, output_pdb_file):
    with open(output_pdb_file, 'w') as file:
        for line in modified_pdb_lines:
            file.write(line)

# Paths to the input files
propka_output_file = 'K62_pH4_propka.pka'
original_pdb_file = 'Hydrogen_Free_K62.pdb'
output_pdb_file = 'K62_pH4_propka.pdb'

# Parse the PROPKA output file
pka_values = parse_propka_output(propka_output_file)

# Update the PDB file with the correct protonation states
modified_pdb_lines = update_pdb_protonation_states(original_pdb_file, pka_values, target_pH=7.0)

# Save the modified PDB file
save_modified_pdb(modified_pdb_lines, output_pdb_file)