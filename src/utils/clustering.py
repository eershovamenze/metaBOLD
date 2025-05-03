import csv
import os
import subprocess
import pandas as pd
from datetime import datetime
import glob

# -------------------------------
# Define project root and directories
# -------------------------------
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
raw_folder = os.path.join(project_root, "data", "processed")
processed_folder = os.path.join(project_root, "data", "processed")
reference_folder = os.path.join(project_root, "data", "reference")
tracking_file_path = os.path.join(raw_folder, "assigned_motu_files.txt")

# Ensure processed folder exists
os.makedirs(processed_folder, exist_ok=True)

# -------------------------------
# Define file paths for reference data
# -------------------------------
reference_fasta = os.path.join(reference_folder, "OTU_reference.fasta")
otu_csv_path = os.path.join(reference_folder, "OTU_database_updated.csv")

# -------------------------------
# Update system PATH for vsearch
# -------------------------------
os.environ["PATH"] += os.pathsep + "/usr/local/bin"

# -------------------------------
# Helper functions
# -------------------------------
def run_vsearch_on_sequence(seq_fasta, db_fasta, output_file, identity=0.97):
    command = [
        'vsearch',
        '--usearch_global', seq_fasta,
        '--db', db_fasta,
        '--id', str(identity),
        '--blast6out', output_file,
        '--maxaccepts', '1'
    ]
    subprocess.run(command, check=True)

def get_next_motu_id(otu_df):
    existing_ids = otu_df['MIN']
    numbers = []
    for motu in existing_ids:
        try:
            num = int(motu.split("_")[-1])
            numbers.append(num)
        except Exception:
            pass
    next_num = max(numbers) + 1 if numbers else 1
    return f"MZP_{next_num:09d}"

# -------------------------------
# Load OTU reference DataFrame
# -------------------------------
if os.path.exists(otu_csv_path):
    otu_df = pd.read_csv(otu_csv_path)
else:
    otu_df = pd.DataFrame(columns=["MIN", "sequence", "timestamp", "source"])

# -------------------------------
# Load or initialize tracking file (to record processed raw files)
# -------------------------------
if os.path.exists(tracking_file_path):
    with open(tracking_file_path, 'r') as f:
        processed_files = {line.strip() for line in f if line.strip()}
else:
    processed_files = set()




# -------------------------------
# Process each raw CSV file that hasn't been processed yet
# -------------------------------
raw_csv_files = glob.glob(os.path.join(raw_folder, "*motu_list.csv"))

for file_path in raw_csv_files:
    file_name = os.path.basename(file_path)
    if file_name in processed_files:
        print(f"Skipping already processed file: {file_name}")
        continue

    print(f"Processing file: {file_name}")
    source_file = os.path.splitext(file_name)[0]

    # Read the CSV file into a DataFrame
    new_seqs_df = pd.read_csv(file_path)
    new_seqs_df["MIN"] = ""  # New column for assigned MOTU IDs

    # Process each sequence (each row)
    for index, row in new_seqs_df.iterrows():
        # Create a unique temporary FASTA and output file for each sequence
        temp_fasta = os.path.join(project_root, f"temp_seq_{index}.fasta")
        output_file = os.path.join(project_root, f"match_output_{index}.txt")

        seq = row['sequence']
        seq_id = row['id']
        with open(temp_fasta, 'w') as f:
            f.write(f">{seq_id}\n{seq}\n")

        # Run vsearch using the reference FASTA file from the reference folder
        if os.path.exists(output_file):
            os.remove(output_file)
        run_vsearch_on_sequence(temp_fasta, reference_fasta, output_file, identity=0.97)

        # Parse the vsearch output for a match
        assigned_motu = None
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                lines = f.readlines()
                if lines:
                    parts = lines[0].strip().split("\t")
                    if len(parts) >= 2:
                        assigned_motu = parts[1]

        # If no match is found, assign a new MOTU ID and update the OTU reference
        if assigned_motu is None:
            assigned_motu = get_next_motu_id(otu_df)
            current_timestamp = datetime.now().isoformat()
            new_row = {
                'MIN': assigned_motu,
                'sequence': seq,
                'timestamp': current_timestamp,
                'source': source_file  # using base filename without extension
            }
            otu_df = pd.concat([otu_df, pd.DataFrame([new_row])], ignore_index=True)
            # Append the new OTU to the reference FASTA file
            with open(reference_fasta, 'a') as fasta_out:
                fasta_out.write(f">{assigned_motu}\n{seq}\n")

        # Update the new_seqs_df DataFrame with the assignment
        new_seqs_df.at[index, "MIN"] = assigned_motu
        print(f"Processed sequence {index} in {file_name}: Assigned to {assigned_motu}")

        # Remove temporary files
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
        if os.path.exists(output_file):
            os.remove(output_file)

    # Save updated OTU reference CSV
    otu_df.to_csv(otu_csv_path, index=False)

    # Save the processed CSV file in the processed folder
    processed_output_name = f"{source_file}_assigned.csv"
    processed_output_path = os.path.join(processed_folder, processed_output_name)
    new_seqs_df[['id', 'MIN']].to_csv(processed_output_path, index=False)

    # Append the processed file to the tracking file
    with open(tracking_file_path, 'a') as tf:
        tf.write(file_name + "\n")

    print(f"Finished processing {file_name} and saved to {processed_output_path}")

print("All unprocessed files have been processed.")
