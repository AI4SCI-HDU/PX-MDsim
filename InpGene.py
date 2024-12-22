import subprocess
import os
import sys

pdbFile1 = sys.argv[1]
DRUG1 = sys.argv[2]
pdbFile2 = sys.argv[3]
DRUG2 = sys.argv[4]

def generate_text(mole1_file, amount_number1, mole2_file, amount_number2, box_length, box_width, box_height):
    text = f"tolerance 2.0\n"
    text += f"filetype pdb\n"
    text += f"output temp.pdb\n\n"

    amount_float1 = float(amount_number1)
    amount_float2 = float(amount_number2)

    # 将 box_length, box_width, box_height 转换为浮点数
    box_length_float = float(box_length)
    box_width_float = float(box_width)
    box_height_float = float(box_height)

    text += f"structure {mole1_file}\n"
    text += f"  number {amount_float1:.0f}\n"
    text += f"  inside box 0. 0. 0. {box_length_float:.0f}. {box_width_float:.0f}. {box_height_float:.0f}.\n"
    text += "end structure\n\n"

    text += f"structure {mole2_file}\n"
    text += f"  number {amount_float2:.0f}\n"
    text += f"  inside box 0. 0. 0. {box_length_float:.0f}. {box_width_float:.0f}. {box_height_float:.0f}.\n"
    text += "end structure\n"

    return text

def save_to_file(text, filename):
    with open(filename, 'w') as file:
        file.write(text)

def valid_file(file_name):
    return file_name.endswith('.pdb')

def valid_number(input_str):
    try:
        float(input_str)
        return True
    except ValueError:
        return False

def modify_gro_box_size(gro_file, box_length, box_width, box_height):
    new_box_length = float(box_length) / 10
    new_box_width = float(box_width) / 10
    new_box_height = float(box_height) / 10

    # 读取 .gro 文件的所有行
    with open(gro_file, 'r') as file:
        lines = file.readlines()

    # 替换最后一行
    lines[-1] = f"   {new_box_length:.5f}   {new_box_width:.5f}   {new_box_height:.5f}\n"

    # 重新写入 .gro 文件
    with open(gro_file, 'w') as file:
        file.writelines(lines)

def add_residue_names(gro_file, residue_name1, residue_name2):
    with open(gro_file, 'r') as file:
        lines = file.readlines()

    total_atoms = int(lines[1].strip())  # 获取原子总数
    new_lines = [lines[0], lines[1]]  # 保留文件的前两行

    current_residue_name = residue_name1
    previous_residue_number = 0  # 上一个原子的分子序号

    for line in lines[2:2 + total_atoms]:
        residue_number = int(line[:5].strip())  # 当前原子的分子序号

        if residue_number < previous_residue_number:
            current_residue_name = residue_name2  # 切换残基名

        atom_name = line[10:15].strip()  # 获取原子名称
        rest_of_line = line[15:]  # 原子名称后的其余部分

        # 组装新的原子行，确保原子名称从第十五个字符位开始
        modified_line = f"{line[:5]}{current_residue_name}   {atom_name:>4}{rest_of_line}"
        new_lines.append(modified_line)
        previous_residue_number = residue_number

    box_size_line = lines[-1]
    new_lines.append(box_size_line)

    with open(gro_file, 'w') as file:
        file.writelines(new_lines)

if __name__ == "__main__":
    mole1_file = pdbFile1
    while not valid_file(mole1_file):
        print("文件名无效，请输入以 .pdb 结尾的文件名。")
        mole1_file = input("请输入交联单体1的pdb文件名：(xxx.pdb)")

    amount_number1 = input(f"请输入放入{DRUG1}的个数：")
    while not valid_number(amount_number1):
        print("无效的数字，请输入整数。")
        amount_number1 = input(f"请输入放入{DRUG1}的个数：")

    mole2_file = pdbFile2
    while not valid_file(mole2_file):
        print("文件名无效，请输入以 .pdb 结尾的文件名。")
        mole2_file = input("请输入交联单体2的pdb文件名：(xxx.pdb)")

    amount_number2 = input(f"请输入放入{DRUG2}的个数：")
    while not valid_number(amount_number2):
        print("无效的数字，请输入整数。")
        amount_number2 = input(f"请输入放入{DRUG2}的个数：")

    # 分别输入盒子的长、宽、高
    box_length = input("请输入盒子的长度（单位：埃）：")
    while not valid_number(box_length):
        print("无效的数字，请输入整数。")
        box_length = input("请输入盒子的长度（单位：埃）：")

    box_width = input("请输入盒子的宽度（单位：埃）：")
    while not valid_number(box_width):
        print("无效的数字，请输入整数。")
        box_width = input("请输入盒子的宽度（单位：埃）：")

    box_height = input("请输入盒子的高度（单位：埃）：")
    while not valid_number(box_height):
        print("无效的数字，请输入整数。")
        box_height = input("请输入盒子的高度（单位：埃）：")

    text = generate_text(mole1_file, amount_number1, mole2_file, amount_number2, box_length, box_width, box_height)
    save_to_file(text, "mixture.inp")
    print("成功生成mixture.inp")

    # 调用 packmol 运行生成的 inp 文件
    packmol_command = f"packmol < mixture.inp"
    result = subprocess.run(packmol_command, shell=True)
    if result.returncode != 0:
        print("packmol 命令执行失败。")
        sys.exit(1)
    else:
        print("成功运行 packmol")

    editconf_command = f"gmx editconf -f temp.pdb -o temp.gro"
    result = subprocess.run(editconf_command, shell=True)
    if result.returncode != 0:
        print("gmx editconf 命令执行失败。")
        sys.exit(1)
    else:
        print("成功运行 editconf")
        
    modify_gro_box_size("temp.gro", box_length, box_width, box_height)
    print("已更新 .gro 文件的盒子尺寸。")

    # 在gro文件中添加残基名
    residue_name1 = mole1_file.split('.')[0]
    residue_name2 = mole2_file.split('.')[0]
    add_residue_names("temp.gro", residue_name1, residue_name2)
