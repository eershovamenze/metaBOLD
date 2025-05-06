#!/usr/bin/env python3
import sys, os, glob, subprocess, tempfile
from datetime import datetime
import pandas as pd

# ------------------------------------------------------------
# Define project root and directories
# ------------------------------------------------------------
project_root      = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
raw_folder        = os.path.join(project_root, "data", "processed")
processed_folder  = os.path.join(project_root, "data", "processed")
reference_folder  = os.path.join(project_root, "data", "reference")
tracking_file_path = os.path.join(raw_folder, "assigned_motu_files.txt")

# Ensure output/reference dirs exist
os.makedirs(processed_folder,  exist_ok=True)
os.makedirs(reference_folder,  exist_ok=True)

# ------------------------------------------------------------
# Define file paths for reference data
# ------------------------------------------------------------
reference_fasta = os.path.join(reference_folder, "OTU_reference.fasta")
otu_csv_path    = os.path.join(reference_folder, "OTU_database_updated.csv")

# ------------------------------------------------------------
# Update system PATH for vsearch
# ------------------------------------------------------------
os.environ["PATH"] += os.pathsep + "/usr/local/bin"

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def run_vsearch_on_sequence(seq_fasta, db_fasta, output_file, identity=0.97):
    cmd = [
        'vsearch', '--usearch_global', seq_fasta,
        '--db', db_fasta, '--id', str(identity),
        '--blast6out', output_file, '--maxaccepts', '1'
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def get_next_motu_id(otu_df):
    numbers = []
    for motu in otu_df['MIN']:
        try:
            numbers.append(int(motu.split("_")[-1]))
        except Exception:
            pass
    next_num = max(numbers) + 1 if numbers else 1
    return f"MZP_{next_num:09d}"

# ------------------------------------------------------------
def main():
    # Load OTU reference DataFrame
    otu_df = pd.read_csv(otu_csv_path) if os.path.exists(otu_csv_path) \
             else pd.DataFrame(columns=["MIN", "sequence", "timestamp", "source"])

    # Load or initialize tracking file
    if os.path.exists(tracking_file_path):
        with open(tracking_file_path) as f:
            processed_files = {ln.strip() for ln in f if ln.strip()}
    else:
        processed_files = set()

    # Process each raw CSV file that hasn't been processed yet
    raw_csv_files = glob.glob(os.path.join(raw_folder, "*motu_list.csv"))

    for file_path in raw_csv_files:
        file_name = os.path.basename(file_path)
        if file_name in processed_files:
            print(f"Skipping already processed file: {file_name}")
            continue

        print(f"Processing file: {file_name}")
        source_file = os.path.splitext(file_name)[0]

        new_seqs_df = pd.read_csv(file_path).assign(MIN="")

        for idx, row in new_seqs_df.iterrows():
            seq_id, seq = row['id'], row['sequence']

            # temp FASTA + output handled in /tmp
            with tempfile.NamedTemporaryFile('w+', suffix=".fasta", delete=False) as temp_fa, \
                 tempfile.NamedTemporaryFile('w+', delete=False)                as temp_out:

                temp_fa.write(f">{seq_id}\n{seq}\n")
                temp_fa.flush()

                run_vsearch_on_sequence(temp_fa.name, reference_fasta,
                                        temp_out.name, identity=0.97)

                # Parse vsearch output
                assigned_motu = None
                temp_out.seek(0)
                line = temp_out.readline().strip()
                if line:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        assigned_motu = parts[1]

            # If no match, mint new MIN
            if assigned_motu is None:
                assigned_motu = get_next_motu_id(otu_df)
                otu_df.loc[len(otu_df)] = {
                    'MIN': assigned_motu, 'sequence': seq,
                    'timestamp': datetime.now().isoformat(timespec="seconds"),
                    'source': source_file
                }
                with open(reference_fasta, 'a') as fa:
                    fa.write(f">{assigned_motu}\n{seq}\n")

            new_seqs_df.at[idx, "MIN"] = assigned_motu

            # clean up temp files
            os.unlink(temp_fa.name)
            os.unlink(temp_out.name)

            if idx % 500 == 0:
                print(f"  …{idx+1}/{len(new_seqs_df)} sequences done")

        # Save updated OTU DB
        otu_df.to_csv(otu_csv_path, index=False)

        # Save per-experiment assignment table
        out_csv = os.path.join(processed_folder, f"{source_file}_assigned.csv")
        new_seqs_df[['id', 'MIN']].to_csv(out_csv, index=False)
        print(f"  wrote {out_csv}")

        # record file as processed
        with open(tracking_file_path, 'a') as tf:
            tf.write(file_name + "\n")

    print("All unprocessed files have been processed.")

# ------------------------------------------------------------
if __name__ == "__main__":
    main()
