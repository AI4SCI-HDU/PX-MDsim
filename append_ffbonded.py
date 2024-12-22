import datetime
import sys

prm_path = f'{sys.argv[1]}.prm'


def copy_sections_with_timestamp_preserving_format(src_file, target_file, sections):
    """
    Copies sections from a source file to a target file, preserving the original formatting.
    Appends a timestamp and a custom note to each line. Ignores lines starting with ';'.
    """
    # Get current Beijing time
    beijing_time = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8)))
    timestamp = beijing_time.strftime("%Y-%m-%d %H:%M:%S")

    # Read the source file and extract sections
    with open(src_file, 'r') as file:
        src_lines = file.readlines()

    # Dictionary to store sections by name and occurrence
    src_sections = {section: [] for section in sections}
    current_section = None
    section_occurrences = {section: 0 for section in sections}

    for line in src_lines:
        # Check for section start and track occurrences
        if any(section in line for section in sections):
            current_section = [section for section in sections if section in line][0]
            section_occurrences[current_section] += 1
            src_sections[current_section].append((section_occurrences[current_section], []))
        elif line.startswith(';') or line.strip() == '':
            continue  # Ignore comments and empty lines
        elif current_section:
            modified_line = line.rstrip() + f" ;added by PX-MDsim {timestamp}\n"
            src_sections[current_section][-1][1].append(modified_line)

    # Read the target file and identify sections to append parameters
    with open(target_file, 'r') as file:
        target_lines = file.readlines()

    # Updated lines for the target file
    updated_target_lines = []
    section_counts = {section: 0 for section in sections}

    for line in target_lines:
        updated_target_lines.append(line)

        # Check if the line is a section header
        if any(section in line for section in sections):
            current_section = [section for section in sections if section in line][0]
            section_counts[current_section] += 1

            # Append the section from the source file if it's the correct occurrence
            for section_occurrence, content in src_sections[current_section]:
                if section_occurrence == section_counts[current_section]:
                    updated_target_lines.extend(content)

    # Write the updated lines to the target file
    with open(target_file, 'w') as file:
        file.writelines(updated_target_lines)


# File paths assuming the script and files are in the same directory
src_file_path = prm_path
target_file_path = 'ffbonded.itp'

# Define sections to be copied
sections_to_copy = ['[ bondtypes ]', '[ angletypes ]', '[ dihedraltypes ]']

# Execute the function
copy_sections_with_timestamp_preserving_format(src_file_path, target_file_path, sections_to_copy)




