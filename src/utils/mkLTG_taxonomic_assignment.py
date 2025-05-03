import os
import subprocess

os.environ["PATH"] = "/usr/local/ncbi/blast/bin" + os.pathsep + os.environ["PATH"]

# Determine the project root (assumes this script is in scripts/python/)
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(project_root)
print("Project root:", project_root)

# Build absolute paths for each file/directory required in the command.
mkLTG_script = os.path.join(project_root, "external/mkLTG/scripts/mkLTG.pl")
input_fasta = os.path.join(project_root, "data/reference/OTU_reference.fasta")
taxonomy_file = os.path.join(project_root, "external/COInr/COInr_for_vtam_2024_05_06_dbV5/COInr_for_vtam_taxonomy.tsv")
blast_db = os.path.join(project_root, "external/COInr/COInr_for_vtam_2024_05_06_dbV5/COInr_for_vtam")
outdir = os.path.join(project_root, "data/reference")
ltg_params = os.path.join(project_root, "external", "mkLTG", "params", "params_Xavier_Turon.tsv")
outname = "Taxonomy_MIN_global"

# Construct the command with the absolute paths.
command = [
    "perl",
    mkLTG_script,
    "-in", input_fasta,
    "-taxonomy", taxonomy_file,
    "-blast_db", blast_db,
    "-outdir", outdir,
    "-out_name", outname,
    "-ltg_params", ltg_params
]

# Run the command.
result = subprocess.run(command, capture_output=True, text=True)

# Check the result.
if result.returncode == 0:
    print("Perl script executed successfully!")
    print(result.stdout)
else:
    print("Error executing the Perl script:")
    print(result.stderr)