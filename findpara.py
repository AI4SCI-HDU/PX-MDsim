import re
import sys

def extract_error_lines_from_stderr(stderr):
    error_pattern = re.compile(r'ERROR \d+ \[file temp\.top, line (\d+)\]')
    error_line_numbers = []

    for line in stderr.splitlines():
        match = error_pattern.search(line)
        if match:
            error_line_numbers.append(int(match.group(1)))
    return error_line_numbers

def find_atoms_section(lines):
    start_idx = end_idx = None
    for i, line in enumerate(lines):
        if line.strip() == '[ atoms ]':
            start_idx = i + 1
        elif start_idx and line.strip() == '' and not end_idx:
            end_idx = i
            break
    return start_idx, end_idx

def create_atom_mapping(lines, atoms_start, atoms_end):
    atoms_section = lines[atoms_start:atoms_end]
    atom_mapping = {}
    for line in atoms_section:
        line = line.strip()
        if line and not line.startswith(';'):
            parts = line.split()
            if len(parts) > 2:
                serial_number = int(parts[0])
                atom_type = parts[1]
                atom_mapping[serial_number] = atom_type
    return atom_mapping

def replace_line_content(line, atom_mapping):
    parts = line.split()
    if parts[-1] == '1':
        parts = parts[:-1]  # Remove the last element if it's '1'
    replaced_content = ' '.join(atom_mapping.get(int(part), part) for part in parts)
    return replaced_content

if __name__ == "__main__":
    # 读取标准错误输入
    stderr_input = sys.stdin.read()

    # 提取错误行号
    error_line_numbers = extract_error_lines_from_stderr(stderr_input)

    # 读取 temp.top 文件
    with open('temp.top', 'r', encoding='utf-8') as file:
        lines = file.readlines()

    # 找到 [ atoms ] 部分
    atoms_start, atoms_end = find_atoms_section(lines)

    # 创建从序列号到原子类型的映射
    atom_mapping = create_atom_mapping(lines, atoms_start, atoms_end)

    # 替换错误行内容并去除重复行
    replaced_lines = {}
    seen_content = set()
    for line_num in error_line_numbers:
        original_line = lines[line_num - 1].strip()
        replaced_content = replace_line_content(original_line, atom_mapping)
        if replaced_content not in seen_content:
            replaced_lines[line_num] = replaced_content
            seen_content.add(replaced_content)

    # 输出替换后的行
    for line_num, replaced_content in replaced_lines.items():
        print(f"Line {line_num}: {replaced_content}")

