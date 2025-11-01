from pathlib import Path

base_dir = Path.cwd()  # current working directory

for species_folder in base_dir.iterdir():
    if not species_folder.is_dir():
        continue

    module2_folder = species_folder / "module2"
    module3_folder = species_folder / "module3"
    module3_folder.mkdir(exist_ok=True)

    altorf_file = module2_folder / "altprot.fasta"
    if not altorf_file.exists():
        print(f"No altprot.fasta found for {species_folder.name}, skipping.")
        continue

    # Read all ORF headers from altprot.fasta
    headers = []
    with open(altorf_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line[1:].strip())  # remove > and newline

    # Select every 1000th ORF (0-indexed)
    candidate_headers = headers[::1000]
    candidate_headers = candidate_headers[0:200]

    # Save to candidate_list.txt in module3
    candidate_file = module3_folder / "candidate_list.txt"
    with open(candidate_file, "w") as f:
        for header in candidate_headers:
            f.write(header + "\n")

    print(f"Candidate list created for {species_folder.name} with {len(candidate_headers)} ORFs.")