# count_atoms.py

import sys

def count_atoms(filename):
    try:
        atom_count = 0
        with open(filename, 'r') as file:
            for line in file:
                # 假设每个原子的数据行都以某个特定关键字开始，例如"ATOM"
                if line.strip().startswith("ATOM"):
                    atom_count += 1
        return atom_count
    except Exception as e:
        print(f"Error reading file: {e}")
        return -1  # 返回 -1 以指示错误

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_atoms.py <strname>")
        sys.exit(1)

    filename = sys.argv[1]
    number = count_atoms(filename)
    print(number)  # 输出原子数量，以便在 shell 中捕获
