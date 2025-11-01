import subprocess
from pathlib import Path
import gzip
import shutil
import time
import csv


# Ensure that transcript IDs and peptide IDs match.
# If they are different, change the peptide FASTA headers so that the transcript ID is used as the identifier
# (You can also change transcriptome file so that peptide ID is used as the identifier).
# For example:
# zcat module2/Homo_sapiens.GRCh38.pep.all.fa.gz | \
# awk '/^>/{match($0, /transcript:([^ ]+)/, a); print ">"a[1]; next} {print}' | \
# gzip > module2/Homo_sapiens.GRCh38.pep.all.fa.gz.tmp && \
# mv module2/Homo_sapiens.GRCh38.pep.all.fa.gz.tmp module2/Homo_sapiens.GRCh38.pep.all.fa.gz


base_dir = Path.cwd()  # current working directory
csv_file = base_dir / "module1_ORF_stats.csv"

# Prepare CSV header
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Species", "Duration_sec", "Output_size_MB", "Num_ORFs"])

for species_folder in base_dir.iterdir():
    if not species_folder.is_dir():
        continue

    module1_folder = species_folder / "module1"
    transcript_files = list(module1_folder.glob("*.fa*"))
    if not transcript_files:
        print(f"No transcript file found in {module1_folder}")
        continue

    transcript_file = transcript_files[0]

    # Decompress gzipped transcript if needed
    if transcript_file.suffix == ".gz":
        transcript_fa = transcript_file.with_suffix("")  # remove .gz
        with gzip.open(transcript_file, "rt") as f_in, open(transcript_fa, "w") as f_out:
            f_out.write(f_in.read())
    else:
        transcript_fa = transcript_file

    print(f"Running MosaicProt detect_ORFs for {species_folder.name}...")
    start_time = time.time()

    cmd = [
        "mosaicprot",
        "detect_ORFs",
        "--transcriptome_file", str(transcript_fa),
        "--threshold", "30",
        "--output_file_type", "fasta"
    ]

    try:
        subprocess.run(cmd, cwd=module1_folder, check=True)
    except subprocess.CalledProcessError as e:
        print(f"MosaicProt failed for {species_folder.name}: {e}")
        continue

    duration = time.time() - start_time

    # Find the newly generated fasta file
    fasta_files = list(module1_folder.glob("*.fasta"))
    if not fasta_files:
        print(f"No ORF fasta output found for {species_folder.name}")
        continue

    default_output = fasta_files[0]  # pick the first fasta found
    desired_output = module1_folder / f"{species_folder.name}_ORFs.fasta"
    shutil.move(default_output, desired_output)

    # Calculate file size in MB
    output_size = desired_output.stat().st_size / (1024 * 1024)

    # Count number of ORFs
    def count_orfs(fasta_file):
        count = 0
        with open(fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
        return count

    num_orfs = count_orfs(desired_output)

    # Save stats to CSV
    with open(csv_file, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([species_folder.name, round(duration, 2),
                         round(output_size, 2), num_orfs])

    print(f"Finished {species_folder.name} in {round(duration,2)}s")
    print(f"Output: {num_orfs} ORFs, {round(output_size,2)} MB, saved to {desired_output}\n")
