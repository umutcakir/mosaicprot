import subprocess
from pathlib import Path
import gzip
import time
import csv

base_dir = Path.cwd()  # current working directory
csv_file = base_dir / "species_ORF_stats.csv"

# Prepare CSV header
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Species", "Duration_sec", "Refprot_size_MB", "Altprot_size_MB",
        "Num_ref_ORFs", "Num_alt_ORFs"
    ])

def count_orfs(fasta_file):
    count = 0
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

def clean_reference_fasta(ref_in, ref_out):

    with open(ref_in, "r") as fin, open(ref_out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if "transcript:" in line:
                    tid = line.split("transcript:")[1].split()[0]
                    fout.write(f">{tid}\n")
                else:
                    # fallback: use first word after '>'
                    fout.write(f">{line[1:].split()[0]}\n")
            else:
                fout.write(line)

for species_folder in base_dir.iterdir():
    if not species_folder.is_dir():
        continue

    module1_folder = species_folder / "module1"
    module2_folder = species_folder / "module2"
    module2_folder.mkdir(exist_ok=True)

    # Find ORFeome file in module1
    orfome_files = list(module1_folder.glob("*_ORFs.fasta"))
    if not orfome_files:
        print(f" No ORFeome file found for {species_folder.name}, skipping.")
        continue
    orfome_file = orfome_files[0]

    # Find reference proteome in module2
    ref_files = list(module2_folder.glob("*.pep*.fa*"))
    if not ref_files:
        print(f" No reference proteome found in {module2_folder}, skipping.")
        continue
    ref_file = ref_files[0]

    # Decompress reference proteome if gzipped
    if ref_file.suffix == ".gz":
        ref_fa = ref_file.with_suffix("")  # remove .gz
        print(f" Decompressing {ref_file.name}...")
        with gzip.open(ref_file, "rt") as f_in, open(ref_fa, "w") as f_out:
            f_out.write(f_in.read())
    else:
        ref_fa = ref_file

    # Clean headers
    cleaned_ref = module2_folder / "reference_cleaned.fasta"
    print(f"Cleaning reference FASTA headers for {species_folder.name}...")
    clean_reference_fasta(ref_fa, cleaned_ref)

    # Output files
    output_ref = module2_folder / "refprot.fasta"
    output_alt = module2_folder / "altprot.fasta"

    print(f" Running MosaicProt for {species_folder.name}...")
    start_time = time.time()

    cmd = [
        "mosaicprot",
        "separate_ORFs",
        "--ORFeome_file", str(orfome_file),
        "--reference_proteome_file", str(cleaned_ref),
        "--output_ref_file", str(output_ref),
        "--output_alt_file", str(output_alt)
    ]

    try:
        subprocess.run(cmd, cwd=module2_folder, check=True)
    except subprocess.CalledProcessError as e:
        print(f"MosaicProt failed for {species_folder.name}: {e}")
        continue

    duration = round(time.time() - start_time, 2)

    # Calculate sizes (MB)
    ref_size = round(output_ref.stat().st_size / (1024 * 1024), 2)
    alt_size = round(output_alt.stat().st_size / (1024 * 1024), 2)

    # Count ORFs
    num_ref_orfs = count_orfs(output_ref)
    num_alt_orfs = count_orfs(output_alt)

    # Write to CSV
    with open(csv_file, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            species_folder.name, duration,
            ref_size, alt_size,
            num_ref_orfs, num_alt_orfs
        ])

    print(f"Finished {species_folder.name} in {duration}s")
    print(f"     Refprot: {num_ref_orfs} ORFs, {ref_size} MB")
    print(f"     Altprot: {num_alt_orfs} ORFs, {alt_size} MB\n")
