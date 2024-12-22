import sys

strFile = sys.argv[1]


def find_amide_bonds(file_path):
    # Initialize dictionaries and lists
    atom_info = {}  # Stores atom types by atom id
    co_bonds = []  # Stores pairs of bonded C and O atoms
    cn_bonds = []  # Stores pairs of bonded C and N atoms, where C is also bonded to O

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

    return NewType_C, NewType_N

# input .str file path
file_path = strFile
NewType_C, NewType_N = find_amide_bonds(file_path)

print('Type_Carboxy_C_New: ' + (', '.join(NewType_C) if NewType_C else 'None'))
print('Type_Amine_N_New: ' + (', '.join(NewType_N) if NewType_N else 'None'))

#print(f'Type_Carboxy_C_New:{NewType_C}')
#print(f'Type_Amine_N_New:{NewType_N}')



