import re
import sys

strFile = sys.argv[1]
Res = sys.argv[2]

def read_input_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    return lines

def is_valid_line(line):
    return line.strip() and line.strip().startswith(("ATOM", "BOND", "IMPR"))

def parse_input(lines):
    atoms = []
    bonds = []
    impropers = []

    for line in lines:
        if is_valid_line(line):
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) >= 6:
                    _, atom, atom_type, charge, _, _ = parts
                    atoms.append([atom, atom_type, charge])
            elif line.startswith("BOND"):
                bond_match = re.match(r"BOND\s+(\S+)\s+(\S+)", line)
                if bond_match:
                    atom1, atom2 = bond_match.groups()
                    bonds.append(f"{atom1} {atom2}")
            elif line.startswith("IMPR"):
                impr_parts = line.split()[1:]  # 提取IMPR行的原子名部分
                num_atoms = len(impr_parts)
                if num_atoms == 4:
                    atom1, atom2, atom3, atom4 = impr_parts
                    impropers.append(f"{atom1} {atom2} {atom3} {atom4}")

    return atoms, bonds, impropers

def append_to_rtp_file(filename, text_to_append):
    with open(filename, 'a') as file:
        file.write("\n" + text_to_append)

def main(input_filename, output_filename):
    lines = read_input_file(input_filename)
    core_lines = [line for line in lines if is_valid_line(line)]
    atoms, bonds, impropers = parse_input(core_lines)

    title = Res
    output_text = f"[{title}]\n\n"
    output_text += "[atoms]\n"
    for i, atom_info in enumerate(atoms):
        atom, atom_type, charge = atom_info
        # 根据位置手动对齐
        output_text += f"    {atom:<5} {atom_type:<8} {' ' if float(charge) >= 0 else ''}{charge}  {i}\n"
    output_text += "[bonds]\n"
    output_text += "".join(f"    {bond}\n" for bond in bonds)
    output_text += "[impropers]\n"
    output_text += "".join(f"    {improper}\n" for improper in impropers)

    append_to_rtp_file("merged.rtp", output_text)

if __name__ == "__main__":
    input_filename = strFile
    output_filename = "merged.rtp"
    main(input_filename, output_filename)
