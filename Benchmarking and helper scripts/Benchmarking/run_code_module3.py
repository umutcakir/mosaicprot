import subprocess
from pathlib import Path
import shutil
import time
import csv

base_dir = Path.cwd()  # current working directory
csv_file = base_dir / "module3_run_stats.csv"

# CSV header
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Species", "Processor", "Duration_sec", "Output_size_MB", "Num_transcripts"])

for species_folder in base_dir.iterdir():
    if not species_folder.is_dir():
        continue

    module1_folder = species_folder / "module1"
    module2_folder = species_folder / "module2"
    module3_folder = species_folder / "module3"
    module3_folder.mkdir(exist_ok=True)

    candidate_file = module3_folder / "candidate_list.txt"
    refprot_file = module2_folder / "refprot.fasta"
    altprot_file = module2_folder / "altprot.fasta"

    # Transcriptome file in module1
    transcript_files = list(module1_folder.glob("*.fa"))
    if not transcript_files:
        print(f"No transcriptome .fa file found for {species_folder.name}, skipping.")
        continue
    transcript_file = transcript_files[0]

    # Check required files
    if not candidate_file.exists() or not refprot_file.exists() or not altprot_file.exists():
        print(f"Missing required files for {species_folder.name}, skipping.")
        continue

    for proc in range(1, 129):  # processors 1 to 128
        proc_folder = module3_folder / str(proc)
        proc_folder.mkdir(exist_ok=True)

        print(f"Running Module3 for {species_folder.name}, processor {proc}...")

        cmd = [
            "mosaicprot",
            "simulate_chimeric_proteins",
            "--candidate_altProt_list", str(candidate_file),
            "--refProts_file", str(refprot_file),
            "--altProts_file", str(altprot_file),
            "--transcriptome_file", str(transcript_file),
            "--processor_num", str(proc)
        ]

        start_time = time.time()
        try:
            subprocess.run(cmd, cwd=module3_folder, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Module3 failed for {species_folder.name}, processor {proc}: {e}")
            continue

        duration = time.time() - start_time

        # Move and rename output
        output_file = module3_folder / "simulated_chimeric_proteins.fasta"
        if output_file.exists():
            dest_file = proc_folder / f"simulated_chimeric_proteins.fasta"
            shutil.move(str(output_file), dest_file)
        else:
            print(f"No output generated for {species_folder.name}, processor {proc}")
            continue

        # Calculate folder size in MB
        folder_size = sum(f.stat().st_size for f in proc_folder.glob("*")) / (1024 * 1024)

        # Count number of transcripts
        num_transcripts = 0
        with open(dest_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    num_transcripts += 1

        # Save stats to CSV
        with open(csv_file, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([species_folder.name, proc, round(duration,2), round(folder_size,2), num_transcripts])

        print(f"Finished {species_folder.name}, processor {proc}, duration: {round(duration,2)}s, size: {round(folder_size,2)}MB, transcripts: {num_transcripts}")