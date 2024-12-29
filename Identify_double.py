import sys

strFile = sys.argv[1]


def find_amide_bonds(file_path):
    # Initialize dictionaries and lists
    atom_info = {}  # Stores atom types by atom id
    co_bonds = []  # Stores pairs of bonded C and O atoms
    cn_bonds = []  # Stores pairs of bonded C and N atoms, where C is also bonded to O
    cc_bonds = []
    hn_bonds = []


    # Read file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Parse lines
    for line in lines:
        if line.startswith('ATOM'):
            parts = line.split()
            if len(parts) >= 4:
                atom_id, atom_type = parts[1], parts[2]
                atom_info[atom_id] = atom_type
        elif line.startswith('BOND'):
            parts = line.split()
            if len(parts) >= 3:
                atom1_id, atom2_id = parts[1], parts[2]
                if atom1_id in atom_info and atom2_id in atom_info:
                    if (atom_info[atom1_id].startswith('C') and atom_info[atom2_id].startswith('O')) or \
                       (atom_info[atom1_id].startswith('O') and atom_info[atom2_id].startswith('C')):
                        co_bonds.append((atom1_id, atom2_id))
                    if (atom_info[atom1_id].startswith('C') and atom_info[atom2_id].startswith('N')) or \
                       (atom_info[atom1_id].startswith('N') and atom_info[atom2_id].startswith('C')):
                        cn_bonds.append((atom1_id, atom2_id))
                    if (atom_info[atom1_id].startswith('C') and atom_info[atom2_id].startswith('C')):
                        cc_bonds.append((atom1_id, atom2_id))
                    if (atom_info[atom1_id].startswith('H') and atom_info[atom2_id].startswith('N')) or \
                       (atom_info[atom1_id].startswith('N') and atom_info[atom2_id].startswith('H')):
                        hn_bonds.append((atom1_id, atom2_id))

    # Filter CN bonds where C is also bonded to O
    c_with_o = set([pair[0] if atom_info[pair[0]].startswith('C') else pair[1] for pair in co_bonds])
    amide_bonds = [(pair[0], pair[1]) for pair in cn_bonds if pair[0] in c_with_o or pair[1] in c_with_o]

    if not amide_bonds:
        print("No amide bonds found.")
        return [], []

    # Extract atom types for C and N in amide bonds
    amide_atom_types = [(atom_info[pair[0]], atom_info[pair[1]]) for pair in amide_bonds]

    # Separate the types for C and N atoms
    NewType_C = [types[0] if types[0].startswith('C') else types[1] for types in amide_atom_types]
    NewType_N = [types[1] if types[1].startswith('N') else types[0] for types in amide_atom_types]

    # --- New Code Start ---
    # Extract the types of O atoms bonded to the amide C atoms
    # Extract the types of C atoms bonded to the amide C atoms
    O_types = []
    for bond in co_bonds:
        c_atom, o_atom = bond
        if c_atom.startswith('C') and atom_info[c_atom] in NewType_C:
            O_types.append(atom_info[o_atom])
        elif o_atom.startswith('C') and atom_info[o_atom] in NewType_C:
            O_types.append(atom_info[c_atom])

    C_to_NC_types = []
    for bond in cc_bonds:
        c2c_atom, nc_atom = bond
        if atom_info[c2c_atom] in NewType_C:
            C_to_NC_types.append(atom_info[nc_atom])
        elif atom_info[nc_atom] in NewType_C:
            C_to_NC_types.append(atom_info[c2c_atom])

    C2N_types = []
    for bond in cn_bonds:
        n_atom, c2n_atom = bond
        if n_atom.startswith('N') and atom_info[n_atom] in NewType_N and atom_info[c2n_atom] not in NewType_C:
            C2N_types.append(atom_info[c2n_atom])
        elif c2n_atom.startswith('N') and atom_info[c2n_atom] in NewType_N and atom_info[n_atom] not in NewType_C:
            C2N_types.append(atom_info[n_atom])

    H2N_types = []
    for bond in hn_bonds:
        n1_atom, h2n_atom = bond
        if n1_atom.startswith('N') and atom_info[n1_atom] in NewType_N:
            H2N_types.append(atom_info[h2n_atom])
        elif h2n_atom.startswith('N') and atom_info[h2n_atom] in NewType_N:
            H2N_types.append(atom_info[n1_atom])
    # --- New Code End ---

    return NewType_C, NewType_N, O_types, C_to_NC_types, C2N_types, H2N_types

# input .str file path
file_path = strFile
NewType_C, NewType_N, O_types, C_to_NC_types, C2N_types, H2N_types = find_amide_bonds(file_path)

print('Type_Carboxy_C_New: ' + (', '.join(NewType_C) if NewType_C else 'None'))
print('Type_Amine_N_New: ' + (', '.join(NewType_N) if NewType_N else 'None'))
print('Type_Oxygen_Bonded_to_C: ' + (', '.join(O_types) if O_types else 'None'))
print('Type_C_Bonded_to_NC: ' + (', '.join(C_to_NC_types) if C_to_NC_types else 'None'))
print('Type_C_Bonded_to_N: ' + (', '.join(C2N_types) if C2N_types else 'None'))
print('Type_H_Bonded_to_N: ' + (', '.join(H2N_types) if H2N_types else 'None'))

#print(f'Type_Carboxy_C_New:{NewType_C}')
#print(f'Type_Amine_N_New:{NewType_N}')



