import argparse

def append_params_with_modifications(src_file, target_file):
    """
    Appends parameters from a CHARMM .str file to a GROMACS .itp file at corresponding sections,
    with each line indented by 5 spaces and replacing '!' with ';' in the copied content.
    """
    # Read the source file and extract parameters
    with open(src_file, 'r') as file:
        lines = file.readlines()

    # Flags for tracking current section in the .str file
    in_section = False
    # Dictionaries to store parameters by section
    params = {
        'BONDS': [],
        'ANGLES': [],
        'DIHEDRALS': [],
        'IMPROPERS': []
    }

    # Parsing the .str file
    for line in lines:
        if line.startswith('BONDS'):
            in_section = 'BONDS'
        elif line.startswith('ANGLES'):
            in_section = 'ANGLES'
        elif line.startswith('DIHEDRALS'):
            in_section = 'DIHEDRALS'
        elif line.startswith('IMPROPERS'):
            in_section = 'IMPROPERS'
        elif line.strip() in ['END', 'RETURN']:
            in_section = False  # End of section
        elif in_section and line.strip():
            # Modify and add line to the corresponding section
            modified_line = '     ' + line.replace('!', ';')
            params[in_section].append(modified_line)

    # Read the target file and identify sections to append parameters
    with open(target_file, 'r') as file:
        itp_lines = file.readlines()

    # Updated lines for the .itp file
    updated_itp_lines = []
    dihedral_section_count = 0  # Count dihedral sections

    # Parsing and updating the .itp file
    for line in itp_lines:
        updated_itp_lines.append(line)
        if '[ bondtypes ]' in line:
            updated_itp_lines.extend(params['BONDS'])
        elif '[ angletypes ]' in line:
            updated_itp_lines.extend(params['ANGLES'])
        elif '[ dihedraltypes ]' in line:
            dihedral_section_count += 1
            if dihedral_section_count == 1:
                updated_itp_lines.extend(params['DIHEDRALS'])
            elif dihedral_section_count == 2:
                updated_itp_lines.extend(params['IMPROPERS'])

    # Write the updated lines to the target file
    with open(target_file, 'w') as file:
        file.writelines(updated_itp_lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Append parameters from a CHARMM .str file to a GROMACS .itp file.")
    parser.add_argument("src_file", help="Path to the source .str file")
    parser.add_argument("target_file", help="Path to the target .itp file")

    args = parser.parse_args()

    append_params_with_modifications(args.src_file, args.target_file)