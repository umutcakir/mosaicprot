from Bio import SeqIO
from Bio.Align import PairwiseAligner
from multiprocessing import Pool, cpu_count
import os

# ===========================================
# USER PARAMETERS
# ===========================================
ALT_PROTS_FILE = "altProts.fasta"          # Input FASTA with alternative proteins
EXTERNAL_FILE = "OpenProt_db.fasta"        # External FASTA file to compare against
IDENTITY_THRESHOLD = 0.90                  # Minimum percent identity (0.90 = 90%)
THREADS = cpu_count() - 1                  # Number of CPU cores to use
# ===========================================


# Global variables for workers
aligner = None
external_prots = []
threshold = None


def init_worker(ext_prots, thr):
    """Initializer for worker processes (sets globals)."""
    global aligner, external_prots, threshold
    external_prots = ext_prots
    threshold = thr

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5


def compare_altprot(alt_record):
# "Compare one altProt against all sequences in external_prots."
    alt_id = alt_record.id
    alt_seq = str(alt_record.seq)
    for db_id, db_seq in external_prots:
        score = aligner.score(alt_seq, db_seq)
        max_possible = 2 * min(len(alt_seq), len(db_seq))
        identity = score / max_possible
        if identity >= threshold:
            return alt_id
    return None


def compare_fasta_files(alt_file, external_file, thr=0.9, threads=None):
# "Compare altProt sequences to external proteins and write matching IDs to a file."
    if threads is None:
        threads = cpu_count() - 1

    # OUTPUT FILE NAME
    alt_base = os.path.splitext(os.path.basename(alt_file))[0]
    ext_base = os.path.splitext(os.path.basename(external_file))[0]
    thr_label = str(int(thr * 100))
    output_file = f"{alt_base}_vs_{ext_base}_{thr_label}.txt"

    # LOAD SEQUENCES
    altprots = list(SeqIO.parse(alt_file, "fasta"))
    ext_prots = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(external_file, "fasta")]

    # RUN IN PARALLEL
    with Pool(threads, initializer=init_worker, initargs=(ext_prots, thr)) as pool:
        hits = list(filter(None, pool.map(compare_altprot, altprots)))

    # SAVE OUTPUT
    hits = sorted(set(hits))
    with open(output_file, "w") as out:
        for alt_id in hits:
            out.write(f"{alt_id}\n")

    print(f"{len(hits)} altProts with ≥{int(thr * 100)}% identity are detected.")
    print(f"List written to: {output_file}")

# MAIN EXECUTION
if __name__ == "__main__":
    compare_fasta_files(
        ALT_PROTS_FILE,
        EXTERNAL_FILE,
        thr=IDENTITY_THRESHOLD,
        threads=THREADS
    )

 

