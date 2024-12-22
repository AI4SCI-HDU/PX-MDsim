# Define the file path
file_path = "temp.gro"  # 替换为您的文件路径

# Initialize variables
a_line = None
a_index = None
m_value = None

try:
    # Step 1: Read the file into memory
    with open(file_path, 'r+', encoding='utf-8') as file:
        lines = file.readlines()  # Read all lines
        previous_letters = None  # To store the three letters from the previous line

        # Locate the first line meeting the condition (a) and extract m from (a - 1)
        for i in range(2, len(lines)):
            current_line = lines[i].strip()
            # Extract the first column
            parts = current_line.split()
            if len(parts) > 0:
                first_column = parts[0]  # Get the first column
                if len(first_column) > 3:
                    current_letters = first_column[-3:]  # Extract the last three letters
                    # Compare with the previous letters
                    if previous_letters is not None and current_letters != previous_letters:
                        a_line = current_line
                        a_index = i + 1  # 1-based index for the line
                        # Extract the number part of the previous line (a - 1)
                        previous_parts = lines[i - 1].strip().split()
                        if len(previous_parts) > 0:
                            m_value = int(''.join(filter(str.isdigit, previous_parts[0])))  # Extract numeric part
                        break  # Stop searching after the first match
                    # Update the previous letters
                    previous_letters = current_letters

        # Step 2: Modify lines from a_index to second last line
        if a_index is not None and m_value is not None:
            for i in range(a_index - 1, len(lines) - 1):  # From a_index to the second last line
                current_line = lines[i].strip()
                parts = current_line.split()
                if len(parts) > 0:
                    first_column = parts[0]  # Extract the first column
                    if first_column.isdigit() or any(c.isdigit() for c in first_column):
                        # Extract numeric part, add m_value, and reconstruct
                        numeric_part = int(''.join(filter(str.isdigit, first_column)))
                        updated_numeric_part = numeric_part + m_value
                        updated_first_column = str(updated_numeric_part) + first_column[-3:]
                        # Determine the number of spaces based on the length of updated number
                        if len(str(updated_numeric_part)) == 1:
                            leading_spaces = " " * 4
                        elif len(str(updated_numeric_part)) == 2:
                            leading_spaces = " " * 3
                        elif len(str(updated_numeric_part)) == 3:
                            leading_spaces = " " * 2
                        elif len(str(updated_numeric_part)) == 4:
                            leading_spaces = " " * 1
                        else:
                            leading_spaces = ""  # No spaces for five or more digits
                        lines[i] = leading_spaces + current_line.replace(first_column, updated_first_column, 1) + "\n"

            # Step 3: Overwrite the original file with modified lines
            file.seek(0)  # Move to the beginning of the file
            file.writelines(lines)
            file.truncate()  # Remove any remaining content in the original file

            # Output results
            print(f"First match found at line {a_index}: {a_line}")
            print(f"Number extracted from line {a_index - 1}: {m_value}")
            print("File has been updated successfully.")
        else:
            print("No matching line found or unable to extract m.")

except Exception as e:
    print(f"Error occurred: {e}")



