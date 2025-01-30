import sys


def extract_charge(file_path, residue_type):
    """
    Extracts the charge of a given residue type from a .str file.

    Args:
    file_path (str): The path to the .str file.
    residue_type (str): The residue type whose charge value is to be extracted.

    Returns:
    float: The charge value of the residue, None if not found.
    """
    charge = None
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) > 3 and parts[2] == residue_type:
                charge = float(parts[3])
                break
    return charge


def main(residue_new_type_c, residue_new_type_n, residue_type_carboxy_c, residue_type_amine_n, residue_new_type_o2c, residue_new_type_c2c, residue_new_type_c2n, residue_new_type_h2n, residue_type_carboxy_o, residue_type_carboxy_r, residue_type_amine_r, residue_type_amine_h):
    # File paths
    path_two = "TMATAE.str"
    path_C = "TMA.str"
    path_N = "TAE.str"

    # Extract charges
    charge_new_type_c = extract_charge(path_two, residue_new_type_c)
    charge_new_type_n = extract_charge(path_two, residue_new_type_n)
    charge_type_carboxy_c = extract_charge(path_C, residue_type_carboxy_c)
    charge_type_amine_n = extract_charge(path_N, residue_type_amine_n)
    charge_new_type_o2c = extract_charge(path_two, residue_new_type_o2c)
    charge_new_type_c2c = extract_charge(path_two, residue_new_type_c2c)
    charge_new_type_c2n = extract_charge(path_two, residue_new_type_c2n)
    charge_new_type_h2n = extract_charge(path_two, residue_new_type_h2n)
    charge_type_carboxy_o = extract_charge(path_C, residue_type_carboxy_o)
    charge_type_carboxy_r = extract_charge(path_C, residue_type_carboxy_r)
    charge_type_amine_r = extract_charge(path_N, residue_type_amine_r)
    charge_type_amine_h = extract_charge(path_N, residue_type_amine_h)

    # Calculate differences
    difference_c = round(charge_new_type_c - charge_type_carboxy_c, 3)
    difference_n = round(charge_new_type_n - charge_type_amine_n, 3)
    difference_o2c = round(charge_new_type_o2c - charge_type_carboxy_o, 3)
    difference_c2c = round(charge_new_type_c2c - charge_type_carboxy_r, 3)
    difference_c2n = round(charge_new_type_c2n - charge_type_amine_r, 3)
    difference_h2n = round(charge_new_type_h2n - charge_type_amine_h, 3)



    return difference_c, difference_n, difference_o2c, difference_c2c, difference_c2n, difference_h2n


if __name__ == "__main__":
    # Expects four command-line arguments: the residue types
    if len(sys.argv) != 13:
        print("Usage: python script.py <new_type_c> <new_type_n>...")
        sys.exit(1)

    residue_new_type_c = sys.argv[1]
    residue_new_type_n = sys.argv[2]
    residue_type_carboxy_c = sys.argv[3]
    residue_type_amine_n = sys.argv[4]
    residue_new_type_o2c = sys.argv[5]
    residue_new_type_c2c = sys.argv[6]
    residue_new_type_c2n = sys.argv[7]
    residue_new_type_h2n = sys.argv[8]
    residue_type_carboxy_o = sys.argv[9]
    residue_type_carboxy_r = sys.argv[10]
    residue_type_amine_r = sys.argv[11]
    residue_type_amine_h = sys.argv[12]


    results = main(residue_new_type_c, residue_new_type_n, residue_type_carboxy_c, residue_type_amine_n, residue_new_type_o2c, residue_new_type_c2c, residue_new_type_c2n, residue_new_type_h2n, residue_type_carboxy_o, residue_type_carboxy_r, residue_type_amine_r, residue_type_amine_h)

    print(results[0], results[1], results[2], results[3], results[4], results[5])
# This script should be run with the proper command-line arguments to function correctly.
# It is designed to be reusable and adaptable to different sets of residue types as defined in any shell script.

