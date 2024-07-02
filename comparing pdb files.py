from Bio import PDB

def extract_protonation_states(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    protonation_states = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = (chain.id, residue.id[1], residue.resname)
                if residue.resname in ['HIS', 'GLU', 'ASP', 'LYS', 'ARG']:
                    protonation_states[res_id] = residue.resname
    return protonation_states

def compare_protonation_states(pdb_files):
    protonation_states_list = [extract_protonation_states(pdb) for pdb in pdb_files]
    all_residues = set()
    for states in protonation_states_list:
        all_residues.update(states.keys())

    comparison = {res: [states.get(res, '-') for states in protonation_states_list] for res in all_residues}
    return comparison

pdb_files = ['K62_pH4_propka.pdb', 'K62_pH4_H++.pdb', 'K62_pH4_PypKa.pdb']
comparison = compare_protonation_states(pdb_files)

print(f"{'Residue':<20} {'File1':<10} {'File2':<10} {'File3':<10}")
for res, states in comparison.items():
    res_str = f"{res[0]} {res[1]} {res[2]}"
    print(f"{res_str:<20} {states[0]:<10} {states[1]:<10} {states[2]:<10}")