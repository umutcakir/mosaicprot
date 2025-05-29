# What is MosaicProt?

The purpose of MosaicProt is to enable the de novo detection of chimeric peptide sequences produced by programmed ribosomal frameshifting (PRF). The tool can be used for the identification of individual PRF events and candidates for mosaic translation, i.e., sequences produced by multiple PRF events from a single transcript (hence the name MosaicProt). The tool can model potential chimeric sequences by analyzing the positional relationships between overlapping or closely spaced ORFs. It supports cases where altORFs:

- Are embedded within refORFs.
- Partially overlap with refORFs, spanning the border between a refORF and a UTR.
- Are located entirely in a UTR, with simulated “long” frameshift events bridging the non-overlapping regions up to 10 nucleotides apart.

The simulation engine systematically explores multiple reading frames and frameshift values (±1 and ±2 nucleotides) to generate all plausible chimeric peptide combinations between overlapping or adjacent ORFs. For each ORF pair, it performs up to 30 one-nucleotide stepwise iterations to model potential frameshift junctions, simulating chimeric sequences of 40 amino acids in length. Sequences containing premature stop codons are automatically excluded, ensuring that only continuous, MS-compatible peptide models are retained. MosaicProt is a command-line tool designed to:

- Detect chimeric peptides and proteins from all possible ORFs of a defined length range in a transcriptome. A user can provide either a spliced transcriptome or precursor sequences, depending on the scope of a study.
- Separate canonical (refProt) and non-canonical (altProt) sequences by comparing to a reference proteome. The reference proteome can be defined by a user or imported from a database.
- Simulate chimeric protein sequences by combining segments of refProts and altProts. It can also model chimeric sequences from overlapping altProts.

## Available Commands

The tool has three commands, which correspond to the modules described in Figure 1 of the published manuscript (the ORF detector module, the altProt/refProt separator module, and the chimeric modeler module):

### 1. detect_ORFs

This command detects all potential ORFs from a transcriptome FASTA file and translates them into corresponding amino acid sequences.

**Usage:**

> mosaicprot detect_ORFs --transcriptome_file <transcriptome.fasta> [--threshold 30] [--output_file_type fasta]


**Parameters:**

- `--transcriptome_file` (required): Input transcriptome FASTA file.
- `--threshold`: Minimum ORF length in amino acids (default: 30).
- `--output_file_type`: Output format (fasta or xml, default: fasta).

The input file `<transcriptome.fasta>` contains an annotated or a de novo sequenced transcriptome of any organism. The output file name is built by the addition of a prefix “_min_30aa_ORFs.fasta” to the input file name, if the threshold is 30. For example, input file name: “my_transcriptome.fasta”; output file name: “my_transcriptome_min_30aa_ORFs.fasta”. The prefix “min_30aa_ORFs” stands for regions (ORF) that are free from in-frame stop codons and can be in silico translated to altProts longer than 29 aa. If the output file format is chosen to be xml, it cannot be used with the next module, which accepts only FASTA files. This xml option was provided to enable export to other pipelines that require xml files as inputs.

### 2. separate_ORFs

This command separates ORF products into refProts and altProts based on a known or a user-defined reference proteome.

**Usage:**

> mosaicprot separate_ORFs --ORFeome_file <orfs.fasta> --reference_proteome_file <ref.fasta> [--output_alt_file altProts.fasta] [--output_ref_file refProts.fasta]

**Parameters:**

- `--ORFeome_file` (required): Input FASTA file with predicted ORF translations.
- `--reference_proteome_file` (required): Input FASTA file with reference proteins.
- `--output_alt_file`: Output file for alternative proteins (default: altProts.fasta).
- `--output_ref_file`: Output file for reference proteins (default: refProts.fasta).

The input file `<orfs.fasta>` is the output of module 1 (detect_ORFs). The input file `<ref.fasta>` contains a canonical proteome downloaded from a database or a user-defined canonical proteome. The default output file names are altProts.fasta and refProts.fasta.

### 3. simulate_chimeric_proteins

This command simulates chimeric proteins by combining segments of refProts and altProts.

**Usage:**

> mosaicprot simulate_chimeric_proteins \
--transcriptome_file <transcriptome.fasta> \
--refProts_file <refProts.fasta> \
--altProts_file <altProts.fasta> \
--candidate_altProt_list <altProt_candidates.txt> \
[--processor_num 1] \
[--repetition keep_first]

**Parameters:**

- `--transcriptome_file` (required): Transcriptome FASTA file.
- `--refProts_file` (required): FASTA file of reference proteins.
- `--altProts_file` (required): FASTA file of all alternative proteins.
- `--candidate_altProt_list` (required): List of identifiers of selected alternative proteins to include in simulation.
- `--processor_num`: Number of processors to use for parallel execution (default: 1).
- `--repetition`: Strategy to handle duplicates (default: keep_first).
  - `keep_first`: Keep only the first of duplicate chimeric proteins.
  - `keep_all`: Keep all duplicates.
  - `drop_all`: Remove all duplicates.

**Output:** A file named `simulated_chimeric_proteins.fasta` containing the simulated chimeric proteins.

The input file `<transcriptome.fasta>` is the same file that was used as input for the first module (detect_ORFs). The input files `<refProts.fasta>` and `<altProts.fasta>` are the output files of the second module (separate_ORFs). The input file `<altProt_candidates.txt>` contains a user-defined list of altProt identifiers separated by a new line. It may correspond to conserved altProts, MS-supported altProts, or both. Thus, it may contain a subset of altProt identifiers from the file `<altProts.fasta>` or identifiers of the entire set of altProts, depending on the scope of a study.

The number of processors to use for parallel execution (default: 1) can be user-defined based on the availability of processors in the system.

The flag [--repetition] helps deal with duplicated models. Such models arise from altORF-altORF pairs (as opposed to refORF-altORF pairs) where one of the ORFs is considered as a refORF and the other one as an altORF, and then the other way around. The “drop_all” option of handling such duplicated models can be useful when the goal is to exclude chimeric models generated from altORF-altORF pairs.

---

## MosaicProt Installation Guide

### Requirements

Before installing MosaicProt, ensure you have:

1. Python3.6 or higher  
> python --version

2. pip package manager  
> pip --version

3. System Requirements (recommended):  
- 4GB+ RAM for large datasets  
- Multi-core CPU for parallel processing  
- 500MB disk space

### Installation

1. Install from PyPI  
> pip install mosaicprot

2. For development/editable installation
> git clone https://github.com/aliyurtsevenn/mosaicprot.git

> cd mosaicprot

> pip install -e .


MosaicProt was developed to advance research on mosaic translation and programmed ribosomal frameshifting. It has enabled the discovery of chimeric proteins across various transcript types (mRNA, ncRNA, rRNA, tRNA) and is adaptable to any annotated or de novo sequenced transcriptome. For biological context and related studies, see our publication: Çakır et al. (2024, preprint).

---

## Citation

If you use MosaicProt in your research, please cite the following article:

> Umut Çakır, Noujoud Gabed, Ali Yurtseven, Igor Kryvoruchko (2025).  
> *A universal pipeline MosaicProt enables large-scale modeling and detection of chimeric protein sequences for studies on programmed ribosomal frameshifting.*  
> bioRxiv (Cold Spring Harbor Laboratory). https://doi.org/XXXXXXX

---

To report bugs, ask questions, or suggest features, feel free to open an issue on GitHub. Your feedback and citations help us improve and sustain this tool.
