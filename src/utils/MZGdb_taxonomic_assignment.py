import os
import subprocess

# Determine project root (assumes this script is in src/utils/)
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(project_root)  # Optional: change working directory to the project root

# Construct absolute paths for each reference file
otu_to_id = os.path.join(project_root, "data", "reference", "OTU_reference.fasta")
reference_fasta = os.path.join(project_root, "external", "MZGdb", "MZGfasta-coi__MZGdbALL__o00__A.fasta")
taxonomy_file = os.path.join(project_root, "data", "reference", "MZGmothur-coi__MZGdbALL__o00__A.txt")

# Construct the MOTHUR command.
# Note: MOTHUR commands are typically given as a string starting with '#' when passed on the command-line.
mothur_command = (
    f">classify.seqs(fasta={otu_to_id}, "
    f"reference={reference_fasta}, "
    f"taxonomy={taxonomy_file})"
)
command = ["mothur", mothur_command]

print("Executing command:", " ".join(command))

# Run the command using subprocess
result = subprocess.run(command, capture_output=True, text=True)

# Check the result
if result.returncode == 0:
    print("Mothur command ran successfully!")
    print(result.stdout)
else:
    print("Error running Mothur command:")
    print(result.stderr)