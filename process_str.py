import sys

def process_str(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            if line.startswith('ATOM'):
                parts = line.split()[:6]
                file.write(' '.join(parts) + '\n')
            elif line.startswith('BOND'):
                parts = line.split()[:4]
                file.write(' '.join(parts) + '\n')
            else:
                file.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: failed to process str!")
        sys.exit(1)

    file_path = sys.argv[1]
    process_str(file_path)
    print(f"Processed str saved")
