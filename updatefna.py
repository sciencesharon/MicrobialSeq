import os
import glob
import argparse

# Function to read strain info from a text file
def load_strain_info(file_path):
    strain_info = {}
    with open(file_path, "r") as f:
        for line in f:
            if line.strip():  # Skip empty lines
                parts = line.strip().split("\t")
                if len(parts) != 3:
                    print(f"Skipping malformed line: {line.strip()}")
                    continue
                refid, species, taxid = parts
                strain_info[refid] = {"species": species, "taxid": taxid}
    return strain_info

def update_fna_headers(strain_info_file, input_folder):
    # Load strain info
    strain_info = load_strain_info(strain_info_file)

    # Create the "cleaned" output directory inside the input folder
    output_folder = os.path.join(input_folder, "cleaned")
    os.makedirs(output_folder, exist_ok=True)

    # Process each .fna file
    for fna_file in glob.glob(os.path.join(input_folder, "*.fna")):
        filename = os.path.basename(fna_file)
        
        # Extract the RefID (up to the first period after the second underscore)
        refid = filename.split(".")[0].strip()

        if refid in strain_info:
            # Extract the strain name and taxid
            strain_name = strain_info[refid]["species"]
            taxid = strain_info[refid]["taxid"]
            
            # Open the input file and create an output file in the "cleaned" folder
            with open(fna_file, "r") as infile, open(os.path.join(output_folder, filename), "w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        # Replace the header line
                        new_header = f">{strain_name}|kraken:taxid|{taxid}\n"
                        outfile.write(new_header)
                    else:
                        # Write the sequence lines as is
                        outfile.write(line)
            print(f"Updated header for {filename}")
        else:
            print(f"Warning: RefID {refid} not found in strain info file.")

# Main function to parse arguments and run the script
def main():
    parser = argparse.ArgumentParser(description="Update headers in .fna files using strain info from a text file.")
    parser.add_argument("strain_info", help="Path to the strain info file (tab-separated).")
    parser.add_argument("input_folder", help="Path to the folder containing .fna files.")

    args = parser.parse_args()

    # Run the update function with the provided arguments
    update_fna_headers(args.strain_info, args.input_folder)

if __name__ == "__main__":
    main()
