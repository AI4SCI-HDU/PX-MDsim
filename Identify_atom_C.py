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

def identify_carboxy_groups(connection_map):
    """Identify carboxy C, carboxy O, hydroxy O, hydroxy H, and carboxy R."""
    carboxy_C = set()
    carboxy_O = set()
    hydroxy_O = set()
    hydroxy_H = set()
    carboxy_R = set()

    # Identifying carboxy C
    for atom, connections in connection_map.items():
        if atom.startswith('C') and sum(1 for conn in connections if conn.startswith('O')) == 2:
            carboxy_C.add(atom)

    # Identifying carboxy O, hydroxy O, hydroxy H, and carboxy R
    for c_atom in carboxy_C:
        oxygens = [o for o in connection_map[c_atom] if o.startswith('O')]
        non_oxygens = [atom for atom in connection_map[c_atom] if not atom.startswith('O')]
        carboxy_R.update(non_oxygens)
        for o_atom in oxygens:
            other_connections = connection_map[o_atom] - {c_atom}
            if all(conn.startswith('C') for conn in other_connections):
                carboxy_O.add(o_atom)
            else:
                hydroxy_O.add(o_atom)
                for conn in other_connections:
                    if conn.startswith('H'):
                        hydroxy_H.add(conn)

    return carboxy_C, carboxy_O, hydroxy_O, hydroxy_H, carboxy_R

def check_and_assign_types(groups, atom_type_map):
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

    carboxy_C, carboxy_O, hydroxy_O, hydroxy_H, carboxy_R = identify_carboxy_groups(connection_map)
    identified_groups = {
        "Carboxy C": carboxy_C,
        "Carboxy O": carboxy_O,
        "Hydroxy O": hydroxy_O,
        "Hydroxy H": hydroxy_H,
        "Carboxy R": carboxy_R
    }

    types, error = check_and_assign_types(identified_groups, atom_type_map)
    if not error:
        Type_Carboxy_C = types.get("Carboxy C", None)
        Type_Carboxy_O = types.get("Carboxy O", None)
        Type_Hydroxy_O = types.get("Hydroxy O", None)
        Type_Hydroxy_H = types.get("Hydroxy H", None)
        Type_Carboxy_R = types.get("Carboxy R", None)

        print("Type_Carboxy_C:", Type_Carboxy_C)
        print("Type_Carboxy_O:", Type_Carboxy_O)
        print("Type_Hydroxy_O:", Type_Hydroxy_O)
        print("Type_Hydroxy_H:", Type_Hydroxy_H)
        print("Type_Carboxy_R:", Type_Carboxy_R)
    else:
        print("Types cannot be assigned due to an error.")

# Example usage
file_path = strFile
main(file_path)



