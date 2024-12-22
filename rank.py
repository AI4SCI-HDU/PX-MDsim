import sys

def process_gro(input_filepath):
    with open(input_filepath, 'r') as infile:
        lines = infile.readlines()

    header = lines[0:2]  # 前两行是头部
    atoms = lines[2:-1]  # 中间是原子信息
    footer = lines[-1]   # 最后一行是尾部

    output_lines = header
    current_molecule = None
    atom_index = 1

    for line in atoms:
        molecule_id = line[0:5].strip()  # 分子ID，例如 "1TMA"
        atom_name = line[10:15].strip()  # 原子名称，例如 "C"
        
        if molecule_id != current_molecule:
            current_molecule = molecule_id
            atom_index = 1
        
        # 格式化原子名称和编号
        formatted_atom_name = f"{atom_name}{atom_index}"
        atom_index += 1
        
        # 确保列对齐，不改变第三列序号
        new_line = f"{line[:10]}{formatted_atom_name:<5}{line[15:]}"
        output_lines.append(new_line)

    output_lines.append(footer)

    with open(input_filepath, 'w') as outfile:
        outfile.writelines(output_lines)
    print(f"Processed file saved to {input_filepath}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python process_gro.py <file_path>")
        sys.exit(1)
    
    input_filepath = sys.argv[1]
    process_gro(input_filepath)
