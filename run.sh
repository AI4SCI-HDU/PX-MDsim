#!/bin/bash

# 确保脚本在遇到任何错误时终止执行
set -e

# WELCOME
echo ""
echo "========================================"
echo "||                                    ||"
echo "||           Welcome to use           ||"
echo "||             PX-MDsim!              ||"
echo "||                                    ||"
echo "========================================"
echo ""

# 提示用户输入并读取输入值
#!/bin/bash

# 读取 Carboxyl monomer 的名称
read -p "Enter the Carboxyl monomer name（XXX）: " CarboxylName

# 为 Carboxyl monomer 自动填充后缀生成文件名
pdbFile1="${CarboxylName}.pdb"
mol2File1="${CarboxylName}.mol2"
strFile1="${CarboxylName}.str"

# 读取 Amino monomer 的名称
read -p "Enter the Amino monomer name（XXX）: " AminoName

# 为 Amino monomer 自动填充后缀生成文件名
pdbFile2="${AminoName}.pdb"
mol2File2="${AminoName}.mol2"
strFile2="${AminoName}.str"

# 生成 Double monomer 的 str 文件名
read -p "Enter the Double monomer str file name: " strFileDouble
strFileDouble="${strFileDouble}.str"

# 打印生成的文件名，以确认
echo "Please confirm the names of input files."
echo "Carboxyl monomer files: $pdbFile1, $mol2File1, $strFile1"
echo "Amino monomer files: $pdbFile2, $mol2File2, $strFile2"
echo "Double monomer str file: $strFileDouble"


# read -p "Enter the Carboxyl monomer pdb file name(XXX.pdb): " pdbFile1
# read -p "Enter the Carboxyl monomer mol2 file name(XXX.mol2): " mol2File1
# read -p "Enter the Carboxyl monomer str file name(XXX.str): " strFile1
# read -p "Enter the Amino monomer pdb file name(XXX.pdb): " pdbFile2
# read -p "Enter the Amino monomer mol2 file name(XXX.mol2): " mol2File2
# read -p "Enter the Amino monomer str file name(XXX.str): " strFile2
# read -p "Enter the Double monomer str file name: " strFileDouble
# CarboxylName=${pdbFile1:0:3}
# AminoName=${pdbFile2:0:3}

sleep 2
echo ""
echo "========================================"
echo "Residue name of Carboxyl set to: $CarboxylName"
echo "Residue name of AminoName set to: $AminoName"
echo "========================================"
echo ""
sleep 2

python3 process_str.py "$strFile1"
python3 process_str.py "$strFile2"
python3 process_str.py "$strFileDouble"


# 使用输入的值执行脚本
python3 InpGene.py "$pdbFile1" "$CarboxylName" "$pdbFile2" "$AminoName"
echo "****************************************"
echo "The initial system has been generated"
echo "****************************************"
echo ""
sleep 2

# 重新排序temp.gro （v1108）

python3 reordergro.py

###

python3 rank.py temp.gro
sleep 1

python3 appendrtp4.py "$strFile1" "$CarboxylName"
python3 appendrtp4.py "$strFile2" "$AminoName"
cp merged.rtp charmm36-mar2019.ff
echo ""
echo "****************************************"
echo "RTP has been added to the force field file"
echo "****************************************"
echo ""
sleep 2

python3 cgenff_charmm2gmx.py "$CarboxylName" "$mol2File1" "$strFile1" charmm36-mar2019.ff
python3 cgenff_charmm2gmx.py "$AminoName" "$mol2File2" "$strFile2" charmm36-mar2019.ff
echo ""
echo "****************************************"
echo "cgenff_charmm2gmx.py has been run"
echo "****************************************"
echo ""
sleep 2

python3 append_ffbonded.py "$CarboxylName"
python3 append_ffbonded.py "$AminoName"
cp ffbonded.itp charmm36-mar2019.ff/
echo ""
echo "****************************************"
echo "New parameters have been added to the force field"
echo "****************************************"
echo ""
sleep 2

gmx pdb2gmx -f temp.gro -o temp.gro -p temp.top -ff charmm36-mar2019 -water none
echo ""
echo "****************************************"
echo "top file has been added"
echo "****************************************"
echo ""
sleep 2
gmx grompp -f minim.mdp -c temp.gro -p temp.top -o em.tpr -maxwarn 99
gmx mdrun -v -deffnm em 
gmx grompp -f nvt.mdp -c em.gro -p temp.top -o temp.tpr -maxwarn 99
gmx mdrun  -v -deffnm temp  -pme gpu -ntmpi 1 -ntomp 8 -gpu_id 1


python3 ndx_generater.py
echo ""
echo "****************************************"
echo "ndx file has been generated"
echo "****************************************"
echo ""
sleep 2

##################################################################

echo "Start labeling functional groups..."
python3 Identify_atom_C.py "$strFile1"
output=$(python3 Identify_atom_C.py $strFile1)
# 变量的值被单引号或双引号包围
type_carboxy_c=$(echo "$output" | grep "Type_Carboxy_C:" | cut -d ' ' -f2)
type_carboxy_o=$(echo "$output" | grep "Type_Carboxy_O:" | cut -d ' ' -f2)
type_hydroxy_o=$(echo "$output" | grep "Type_Hydroxy_O:" | cut -d ' ' -f2)
type_hydroxy_h=$(echo "$output" | grep "Type_Hydroxy_H:" | cut -d ' ' -f2)
type_carboxy_r=$(echo "$output" | grep "Type_Carboxy_R:" | cut -d ' ' -f2)

# 检查是否成功提取所有值
if [ -z "$type_carboxy_c" ] || [ -z "$type_carboxy_o" ] || [ -z "$type_hydroxy_o" ] || [ -z "$type_hydroxy_h" ] || [ -z "$type_carboxy_r" ]; then
    echo "无法获取官能团的原子类型，请检查可识别的官能团是否存在。"
    exit 1
fi

# 使用sed替换PXLink.py中所有指定的字符串为提取到的相应值
sed -i "s/C_C/$type_carboxy_c/g" PXLink.py
sed -i "s/C_O/$type_carboxy_o/g" PXLink.py
sed -i "s/H_O/$type_hydroxy_o/g" PXLink.py
sed -i "s/H_H/$type_hydroxy_h/g" PXLink.py
sed -i "s/C_R/$type_carboxy_r/g" PXLink.py

sleep 2
python3 Identify_atom_N.py "$strFile2"

output=$(python3 Identify_atom_N.py $strFile2)

# 从输出中提取各种类型的值
type_amine_n=$(echo "$output" | grep "Type_Amine_N:" | cut -d ' ' -f2)
type_amine_h=$(echo "$output" | grep "Type_Amine_H:" | cut -d ' ' -f2)
type_amine_r=$(echo "$output" | grep "Type_Amine_R:" | cut -d ' ' -f2)

# 检查是否成功提取所有值
if [ -z "$type_amine_n" ] || [ -z "$type_amine_h" ] || [ -z "$type_amine_r" ]; then
    echo "无法获取官能团的原子类型，请检查可识别的官能团是否存在。"
    exit 1
fi

# 使用sed替换PXLink.py中所有指定的字符串为提取到的相应值
sed -i "s/A_N/$type_amine_n/g" PXLink.py
sed -i "s/A_H/$type_amine_h/g" PXLink.py
sed -i "s/A_R/$type_amine_r/g" PXLink.py


sleep 2
python3 Identify_double.py "$strFileDouble"

output=$(python3 Identify_double.py $strFileDouble)

# 从输出中提取NewType_C和NewType_N的值
new_type_c=$(echo "$output" | grep "Type_Carboxy_C_New:" | cut -d ' ' -f2)
new_type_n=$(echo "$output" | grep "Type_Amine_N_New:" | cut -d ' ' -f2)
new_type_o2c=$(echo "$output" | grep "Type_Oxygen_Bonded_to_C:" | cut -d ' ' -f2)
new_type_c2c=$(echo "$output" | grep "Type_C_Bonded_to_NC:" | cut -d ' ' -f2)
new_type_c2n=$(echo "$output" | grep "Type_C_Bonded_to_N:" | cut -d ' ' -f2)
new_type_h2n=$(echo "$output" | grep "Type_H_Bonded_to_N:" | cut -d ' ' -f2)
echo "NewType_C: $new_type_c"
echo "NewType_N: $new_type_n"
echo "NewType_o2c: $new_type_o2c"
echo "NewType_c2c: $new_type_c2c"
echo "NewType_c2n: $new_type_c2n"
echo "NewType_h2c: $new_type_h2n"

# 检查是否成功提取了需要的值
if [ -z "$new_type_c" ] || [ -z "$new_type_n" ] || [ -z "$new_type_o2c" ] || [ -z "$new_type_c2c" ] || [ -z "$new_type_c2n" ] || [ -z "$new_type_h2n" ]; then
    echo "无法获取官能团的原子类型，请检查可识别的官能团是否存在。"
    exit 1
fi

# 使用sed替换PXLink.py中所有指定的字符串为提取到的相应值
sed -i "s/N_C/$new_type_c/g" PXLink.py
sed -i "s/N_N/$new_type_n/g" PXLink.py

sed -i "s/MPD/$AminoName/g" PXLink.py
sed -i "s/TMA/$CarboxylName/g" PXLink.py
echo ""
echo "****************************************"
echo "Have updated names in PXLink"
echo "Prepared to crosslink..."
echo "****************************************"
echo ""
sleep 2

####################################################

output=$(python3 get_charge.py "$new_type_c" "$new_type_n" "$type_carboxy_c" "$type_amine_n" "$new_type_o2c" "$new_type_c2c" "$new_type_c2n" "$new_type_h2n" "$type_carboxy_o" "$type_carboxy_r" "$type_amine_r" "$type_amine_h")

# Extract the charge differences from the output
difference_c=$(echo $output | cut -d ' ' -f1)
difference_n=$(echo $output | cut -d ' ' -f2)
difference_o2c = $(echo $output | cut -d ' ' -f3)
difference_c2c = $(echo $output | cut -d ' ' -f4)
difference_c2n = $(echo $output | cut -d ' ' -f5)
difference_h2n = $(echo $output | cut -d ' ' -f6)

# Use sed to replace the first element of the self.delta_charge_C and self.delta_charge_N lists in PXLink.py
sed -i "s/self.delta_charge_C = \[0, 0, 0\]/self.delta_charge_C = \[$difference_c, $difference_o2c, $difference_c2c\]/g" PXLink.py
sed -i "s/self.delta_charge_N = \[0, 0, 0\]/self.delta_charge_N = \[$difference_n, $difference_h2n, $difference_c2n\]/g" PXLink.py

echo "Updated PXLink.py with new charge values."

####################################################

output=$(python3 count_atom.py "$strFile1")
value1=$((output - 2))
value2=$((output - 4))
value3=$((output - 6))

output2=$(python3 count_atom.py "$strFile2")
value4=$((output2 - 1))
value5=$((output2 - 2))

sed -i "s/\[16, 15, 14, 21, 19, 17, 15\]/\[$output2, $value4, $value5, $output, $value1, $value2, $value3\]/g" PXLink.py
sleep 1

python3 example_run_script.py


# python ndx_generater.py "$inputFile" "path/to/another/output.file"
# echo "ndx_generater.py completed."
# 
# ./cgenff_charmm2gmx.py "$DRUG" "$mol2File" "$strFile" charmm36.ff
# echo "cgenff_charmm2gmx.py completed."
# 
# echo "All scripts completed successfully."
# 