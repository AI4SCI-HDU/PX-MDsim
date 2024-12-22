from collections import defaultdict
import sys

strFile = sys.argv[1]

def parse_bond_lines(bond_lines):
    """Parse bond lines and create a connection map."""
    connection_map = defaultdict(set)
    for line in bond_lines:
        parts = line.split()
        if len(parts) >= 3:
            atom1, atom2 = parts[1], parts[2]
            connection_map[atom1].add(atom2)
            connection_map[atom2].add(atom1)
    return connection_map

def parse_atom_lines(atom_lines):
    """Parse atom lines and create a map of atom identifiers to their types."""
    atom_type_map = {}
    for line in atom_lines:
        parts = line.split()
        if len(parts) >= 3:
            atom_id, atom_type = parts[1], parts[2]
            atom_type_map[atom_id] = atom_type
    return atom_type_map

def identify_amine_groups(connection_map):
    """Identify amine N, amine H, and amine R."""
    amine_N = set()
    amine_H = set()
    amine_R = set()

    for atom, connections in connection_map.items():
        if atom.startswith('N') and any(conn.startswith('H') for conn in connections):
            amine_N.add(atom)

    for n_atom in amine_N:
        for connected_atom in connection_map[n_atom]:
            if connected_atom.startswith('H'):
                amine_H.add(connected_atom)
            else:
                amine_R.add(connected_atom)

    return amine_N, amine_H, amine_R

def check_and_assign_types(groups, atom_type_map):
    """Check types and assign to variables."""
    type_dict = {}
    error_occurred = False

    for group_name, atoms in groups.items():
        if not atoms:
            continue

        atom_types = [atom_type_map[atom] for atom in atoms if atom in atom_type_map]
        if len(set(atom_types)) == 1:
            type_dict[group_name] = atom_types[0]
        else:
            print(f"Error: Not all atoms in {group_name} have the same type.")
            error_occurred = True

    return type_dict, error_occurred

def main(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    bond_lines = [line for line in lines if line.startswith('BOND')]
    atom_lines = [line for line in lines if line.startswith('ATOM')]

    connection_map = parse_bond_lines(bond_lines)
    atom_type_map = parse_atom_lines(atom_lines)

    amine_N, amine_H, amine_R = identify_amine_groups(connection_map)
    identified_groups = {
        "Amine N": amine_N,
        "Amine H": amine_H,
        "Amine R": amine_R
    }

    types, error = check_and_assign_types(identified_groups, atom_type_map)
    if not error:
        Type_Amine_N = types.get("Amine N", None)
        Type_Amine_H = types.get("Amine H", None)
        Type_Amine_R = types.get("Amine R", None)

        print("Type_Amine_N:", Type_Amine_N)
        print("Type_Amine_H:", Type_Amine_H)
        print("Type_Amine_R:", Type_Amine_R)
    else:
        print("Types cannot be assigned due to an error.")

# Example usage
file_path = strFile
main(file_path)



