from collections import defaultdict

def parse_gro_with_fixed_positions_and_z_filter(gro_lines, z_lower_limit, z_upper_limit):
    system_atoms_fixed = []
    molecule_groups_fixed = defaultdict(list)
    frozen_atoms = []  # 存放Z轴在指定范围外的原子
    atom_index_start_pos = 15
    atom_index_end_pos = 20
    z_pos_start = 36  # 假设Z轴坐标从第36个字符开始
    z_pos_end = 44

    for line in gro_lines[2:-1]:  # 跳过头部和尾部的非原子行
        if len(line) > atom_index_end_pos:  # 确保行长度足够
            atom_index_str = line[atom_index_start_pos:atom_index_end_pos].strip()
            try:
                atom_index = int(atom_index_str)
                system_atoms_fixed.append(atom_index)
                molecule_name = ''.join(filter(str.isalpha, line[:10])).strip()
                molecule_groups_fixed[molecule_name].append(atom_index)
                
                # 提取Z轴坐标
                z_pos_str = line[z_pos_start:z_pos_end].strip()
                z_pos = float(z_pos_str)
                
                # 检查Z轴坐标是否在上下0.6nm之外，超出范围的原子加入 frozen 组
                if z_pos > z_upper_limit or z_pos < z_lower_limit:
                    frozen_atoms.append(atom_index)

            except ValueError:
                # 如果转换失败，可以在这里记录日志或者静默跳过
                pass
        else:
            # 处理行长度不足的情况，可以记录日志或跳过
            pass

    return system_atoms_fixed, molecule_groups_fixed, frozen_atoms

def generate_ndx_content_fixed(system_atoms, molecule_groups, frozen_atoms):
    def format_index_with_correct_spaces(indices):
        formatted_lines = []
        line = ""
        for i, index in enumerate(indices):
            if i % 15 == 0 and i != 0:
                formatted_lines.append(line.rstrip())
                line = ""
            space = " " * (5 - min(4, len(str(index))))
            line += f"{index}{space}"
        if line:
            formatted_lines.append(line.rstrip())
        return "\n".join(formatted_lines)

    ndx_content_final = ""
    if system_atoms:
        ndx_content_final += "[ System ]\n"
        ndx_content_final += format_index_with_correct_spaces(system_atoms) + "\n\n"
    if system_atoms:
        ndx_content_final += "[ Other ]\n"
        ndx_content_final += format_index_with_correct_spaces(system_atoms) + "\n"

    for molecule_name, atoms in molecule_groups.items():
        if atoms:
            ndx_content_final += f"\n[ {molecule_name} ]\n"
            ndx_content_final += format_index_with_correct_spaces(atoms) + "\n"

    if frozen_atoms:
        ndx_content_final += "\n[ frozen ]\n"
        ndx_content_final += format_index_with_correct_spaces(frozen_atoms) + "\n"

    return ndx_content_final.rstrip()

def remove_unwanted_zeros_and_empty_categories(ndx_content):
    lines = ndx_content.split('\n')
    cleaned_lines = [line for line in lines if line.strip() != '0']

    if cleaned_lines[-1].strip() == '0' and cleaned_lines[-2].strip() == '[  ]':
        cleaned_lines = cleaned_lines[:-2]
    elif cleaned_lines[-1].strip() == '0':
        cleaned_lines = cleaned_lines[:-1]

    return '\n'.join(cleaned_lines)

def remove_last_empty_category(ndx_content):
    lines = ndx_content.split('\n')
    if lines[-1].strip() == '[  ]':
        cleaned_content = '\n'.join(lines[:-1])
    else:
        cleaned_content = ndx_content
    return cleaned_content


# 读取.gro文件内容
gro_file_path = 'temp.gro'
with open(gro_file_path, 'r') as file:
    gro_lines = file.readlines()

# 获取盒子尺寸，通常是最后一行，解析Z轴大小
box_size_line = gro_lines[-1].strip().split()
z_box_size = float(box_size_line[2])  # 假设Z轴尺寸在第三个位置

# 定义Z轴的上下限（小于0.6nm和大于盒子高度减去0.6nm）
z_lower_limit = 1.0
z_upper_limit = z_box_size - 1.0

# 解析.gro文件，生成并格式化.ndx文件内容
system_atoms_fixed, molecule_groups_fixed, frozen_atoms = parse_gro_with_fixed_positions_and_z_filter(gro_lines, z_lower_limit, z_upper_limit)
ndx_content = generate_ndx_content_fixed(system_atoms_fixed, molecule_groups_fixed, frozen_atoms)
ndx_content = remove_unwanted_zeros_and_empty_categories(ndx_content)
ndx_content_final = remove_last_empty_category(ndx_content)

# 写入最终.ndx文件
output_ndx_file_path = 'init2.ndx'
with open(output_ndx_file_path, 'w') as file:
    file.write(ndx_content_final)
