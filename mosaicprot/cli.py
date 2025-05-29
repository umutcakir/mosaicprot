from Bio import SeqIO, BiopythonWarning
from Bio.Seq import Seq
import re
import pandas as pd
import math
import warnings
import os
import glob
from collections import Counter, defaultdict
from multiprocessing import Pool
import time
import argparse
import sys
import shutil
from functools import partial, reduce
from itertools import chain
from typing import Iterator


warnings.filterwarnings("ignore", category=SyntaxWarning)
warnings.filterwarnings("ignore", category=BiopythonWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)
warnings.filterwarnings("ignore", category=BiopythonWarning)


def frameshifter(sequence, frame):
    if frame > 0:
        return str(sequence)[frame - 1:]
    return str(sequence.reverse_complement())[abs(frame) - 1:]


def find_orfs(shifted_seq, min_aa_lim):
    translated_seq = Seq(shifted_seq).translate()
    orf_dict = {}
    stop_positions = [0] + [m.start() + 1 for m in re.finditer(r'\*', str(translated_seq))] + [len(translated_seq) + 1]
    translations = re.split(r'\*', str(translated_seq))
    for i, aa_seq in enumerate(translations):
        if len(aa_seq) >= min_aa_lim:
            orf_dict[aa_seq] = stop_positions[i:i + 2]
    return orf_dict


def update_positions(orf_dict, seq_len, frame):
    if frame > 0:
        for aa_seq, (start, end) in orf_dict.items():
            orf_dict[aa_seq] = [3 * start + frame, 3 * end - 4 + frame]
    else:
        for aa_seq, (start, end) in orf_dict.items():
            orf_dict[aa_seq] = [seq_len - 3 * end + 5 + frame, seq_len - 3 * start + 1 + frame]


def process_fasta_output(input_file, min_aa_lim, output_file, temp_file_name):
    print("Finding ORFs...")
    try:
        log_step(temp_file_name, f"Reading input file: {input_file}")
        with open(output_file, 'w') as file:
            frames = {1: '1F', 2: '2F', 3: '3F'}
            for record in SeqIO.parse(input_file, 'fasta'):
                log_step(temp_file_name, f"Processing ORFs for Record ID: {record.id}")
                for frame, frame_name in frames.items():
                    try:
                        orf_dict = find_orfs(frameshifter(record.seq, frame), min_aa_lim)
                        update_positions(orf_dict, len(record.seq), frame)
                        for aa_seq, (start, end) in orf_dict.items():
                            file.write(f'>{record.id}_{frame_name}_{start}-{end}_{abs(end - start + 1)}\n{aa_seq}\n')
                    except Exception as e:
                        log_step(temp_file_name, f"Error processing record {record.id}: {e}")
    except Exception as e:
        log_step(temp_file_name, f"Error in process_fasta_output: {e}")
        raise


def run_alt_orf(mrna_file, proteome_file, alt_output, ref_output, temp_file_name):
    print("Finding alternative ORFs...")
    SeqDict_mrna = SeqIO.to_dict(SeqIO.parse(open(mrna_file), "fasta"))
    SeqDict_prot = SeqIO.to_dict(SeqIO.parse(open(proteome_file), "fasta"))
    log_step(temp_file_name, f"mRNA and protein data has been read successfully...")

    with open(ref_output, "w") as ref_file, open(alt_output, "w") as alt_file:
        log_step(temp_file_name, f"Writing reference and alternative ORF files...")
        for prot_id, protein_seq in SeqDict_prot.items():
            alt_list = [k for k in SeqDict_mrna if k.startswith(prot_id)]
            for alt_id in alt_list:
                alt_seq = SeqDict_mrna[alt_id].seq.strip()
                if str(protein_seq.seq).strip() in str(alt_seq):
                    ref_file.write(f">{alt_id}\n{alt_seq}\n")
                else:
                    alt_file.write(f">{alt_id}\n{alt_seq}\n")


def generate_mosaic_proteins(input_file, refprots_file, altprots_file, transcript_file, temp_file_name, idx_num):

    files_to_remove = [
        f"clean_seq_{idx_num}.fasta",
        f"temp_file_for_seq_store_{idx_num}.txt",
        f"temp_file_for_seq_store_temp_{idx_num}.txt",
        f"your_final_result_{idx_num}.fasta",
        f"your_final_result_before_clean_{idx_num}.fasta",
        f"clean_final_result_{idx_num}.fasta",
        f"clean_upload_result_{idx_num}.fasta",
        os.path.join(".mosaic_prot_results", f"requested_file_{idx_num}.fasta"),
        os.path.join(".mosaic_prot_results", f"requested_file_filtered{idx_num}.fasta"),
        os.path.join(".mosaic_prot_results", f"temp_sequence_archieve_{idx_num}.fasta")
    ]

    if os.path.getsize(input_file) == 0:
        print(f"Skipping empty chunk {idx_num + 1}")
        log_step(temp_file_name, f"Skipping empty chunk {idx_num+ 1}")
        return

    if idx_num == 0:
        print("Started to find the chimeric proteins...")
    print(f"Waiting for processing the chunk number {idx_num + 1}...")
    results_directory = ".mosaic_prot_results"

    def getList(dict):
        return list(dict.keys())

    def ngram(seq: str, n: int) -> Iterator[str]:
        return (seq[i: i + n] for i in range(0, len(seq) - n + 1))

    def allngram(seq: str) -> set:
        lengths = range(len(seq))
        ngrams = map(partial(ngram, seq), lengths)
        return set(chain.from_iterable(ngrams))

    gencode = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}

    basepairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def translate(sequence):
        translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
        return translate

    def closestMultiple(n, x=3):
        if x > n:
            return x;
        z = (int)(x / 2);
        n = n + z;
        n = n - (n % x);
        return n;

    def join_string(list_string):
        string = '_'.join(list_string)

        return string

    def translate_all_frameshifted(sequence, skip_value=30, altprot_name="noname", refprot_name="noname",
                                   for_loop_number=90, length_second_part_raw=90, type_acc_to_position="not_specified"):
        for_loop_value = int(math.ceil((int(for_loop_number) / 3)))
        to_file_temp = ""
        for i in range(for_loop_value):
            length_second_part = int(length_second_part_raw) - (int(i) * 3)
            nt_seq_minus_1 = sequence[:(skip_value + i * 3)] + sequence[(skip_value + i * 3) - 1:((
                                                                                                          skip_value + i * 3) - 1 + length_second_part) + 3]
            to_file_temp = to_file_temp + str(
                ">" + altprot_name + "_" + refprot_name + "_" + "-1" + "_" + "iteration_" + str(i) + "\n" + translate(
                    nt_seq_minus_1) + "\n")
            nt_seq_minus_2 = sequence[:(skip_value + i * 3)] + sequence[(skip_value + i * 3) - 2:((
                                                                                                          skip_value + i * 3) - 2 + length_second_part) + 3]

            to_file_temp = to_file_temp + str(
                ">" + altprot_name + "_" + refprot_name + "_" + "-2" + "_" + "iteration_" + str(i) + "\n" + translate(
                    nt_seq_minus_2) + "\n")

            nt_seq_plus_1 = sequence[:(skip_value + i * 3)] + sequence[(skip_value + i * 3) + 1:((
                                                                                                         skip_value + i * 3) + 1) + length_second_part]
            to_file_temp = to_file_temp + str(
                ">" + altprot_name + "_" + refprot_name + "_" + "+1" + "_" + "iteration_" + str(i) + "\n" + translate(
                    nt_seq_plus_1) + "\n")

            nt_seq_plus_2 = sequence[:(skip_value + i * 3)] + sequence[(skip_value + i * 3) + 2:((
                                                                                                         skip_value + i * 3) + 2) + length_second_part]
            to_file_temp = to_file_temp + str(
                ">" + altprot_name + "_" + refprot_name + "_" + "+2" + "_" + "iteration_" + str(i) + "\n" + translate(
                    nt_seq_plus_2) + "\n")
            i = i + 1

        temp_file_for_seq_store_temp = open(
            os.path.join(results_directory, f'temp_file_for_seq_store_temp_{idx_num}.txt'), 'w+')
        temp_file_for_seq_store_temp.write(to_file_temp)
        temp_file_for_seq_store_temp.close()

        temp_file_for_seq_store = open(os.path.join(results_directory, f'temp_file_for_seq_store_{idx_num}.txt'), 'w+')
        to_file = ""
        for seq_record in SeqIO.parse(os.path.join(results_directory, f'temp_file_for_seq_store_temp_{idx_num}.txt'),
                                      format="fasta"):
            if len(seq_record) < 42:
                to_file = to_file + str(
                    ">" + str(seq_record.id) + "_" + str(type_acc_to_position) + "\n" + str(seq_record.seq) + "\n")
        temp_file_for_seq_store.write(to_file)
        temp_file_for_seq_store.close()

        return True

    def translate_all_frameshifted_for_within_embedded(sequence, skip_value=30, altprot_name="noname",
                                                       refprot_name="noname", for_loop_number=90,
                                                       length_second_part_raw=90, type_acc_to_position="not_specified"):
        for_loop_value = int(int(for_loop_number) / 3)
        to_file_temp = ""
        for i in range(for_loop_value):
            length_second_part = int(length_second_part_raw) - (int(i) * 3)
            nt_seq_minus_1 = sequence[:(skip_value + i * 3) - 1] + sequence[(skip_value + i * 3):((
                                                                                                          skip_value + i * 3) + length_second_part) + 3]
            to_file_temp = to_file_temp + str(
                ">" + altprot_name + "_" + refprot_name + "_" + "-1" + "_" + "iteration_" + str(i) + "\n" + translate(
                    nt_seq_minus_1) + "\n")
            nt_seq_minus_2 = sequence[:(skip_value + i * 3) - 2] + sequence[(skip_value + i * 3):((
                                                                                                          skip_value + i * 3) + length_second_part) + 3]

            to_file_temp = to_file_temp + str(
                ">" + altprot_name + "_" + refprot_name + "_" + "-2" + "_" + "iteration_" + str(i) + "\n" + translate(
                    nt_seq_minus_2) + "\n")

            nt_seq_plus_1 = sequence[:(skip_value + i * 3)] + sequence[(skip_value + i * 3) + 1:((
                                                                                                         skip_value + i * 3) + 1) + length_second_part]
            to_file_temp = to_file_temp + str(
                ">" + altprot_name + "_" + refprot_name + "_" + "+1" + "_" + "iteration_" + str(i) + "\n" + translate(
                    nt_seq_plus_1) + "\n")

            nt_seq_plus_2 = sequence[:(skip_value + i * 3)] + sequence[(skip_value + i * 3) + 2:((
                                                                                                         skip_value + i * 3) + 2) + length_second_part]
            to_file_temp = to_file_temp + str(
                ">" + altprot_name + "_" + refprot_name + "_" + "+2" + "_" + "iteration_" + str(i) + "\n" + translate(
                    nt_seq_plus_2) + "\n")
            i = i + 1

        temp_file_for_seq_store_temp = open(
            os.path.join(results_directory, f'temp_file_for_seq_store_temp_{idx_num}.txt'), 'w+')
        temp_file_for_seq_store_temp.write(to_file_temp)
        temp_file_for_seq_store_temp.close()

        temp_file_for_seq_store = open(os.path.join(results_directory, f'temp_file_for_seq_store_{idx_num}.txt'), 'w+')
        to_file = ""
        for seq_record in SeqIO.parse(os.path.join(results_directory, f'temp_file_for_seq_store_temp_{idx_num}.txt'),
                                      format="fasta"):
            if len(seq_record) < 42:
                to_file = to_file + str(
                    ">" + str(seq_record.id) + "_" + str(type_acc_to_position) + "\n" + str(seq_record.seq) + "\n")
        temp_file_for_seq_store.write(to_file)
        temp_file_for_seq_store.close()

        return True

    def append_aa_seqs_for_no_overlapped(firstsequence, secondsequence, width,
                                         fiveorthreeprime="five-or-three-not-specified", altprotname="noname",
                                         refprotname="noname", filename=os.path.join(results_directory,
                                                                                     f"temp_file_for_seq_store_{idx_num}.txt")):
        temp_file_for_seq_store = open(filename, 'w+')
        to_file = ""
        to_file = to_file + str(
            ">" + str(altprotname) + "_" + str(refprotname) + "_" + str(fiveorthreeprime) + "width_" + str(
                width) + "\n" + str(str(firstsequence) + str(secondsequence)) + "\n")
        temp_file_for_seq_store.write(to_file)
        temp_file_for_seq_store.close()

        shutil.copy(filename, f"clean_seq_{idx_num}.fasta")

    def determine_frame_from_file(common_seq,
                                  file=os.path.join(results_directory, f"temp_file_for_seq_store_{idx_num}.txt")):
        temp_file_for_seq_store = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
        frameshift_common = []
        for key, value in temp_file_for_seq_store.items():
            if str(common_seq) in str(value.seq):
                frameshift_common.append(key.split("_")[10] + "_" + "iteration")
        if len(frameshift_common) > 0:
            counter_object = Counter(frameshift_common)
            counter_df = pd.DataFrame.from_dict(counter_object, orient='index').reset_index()
            counter_df = counter_df.rename(columns={'index': 'frameshift', 0: 'count'})
            sorted_counter_df = counter_df.sort_values("count", ascending=False)
            list_of_frameshifts = sorted_counter_df["frameshift"].to_list()
        else:
            list_of_frameshifts = []
        return list_of_frameshifts

    def get_only_frame(list_of_frameshifts,
                       file=os.path.join(results_directory, f"temp_file_for_seq_store_{idx_num}.txt")):

        for list_of_frameshift in list_of_frameshifts:
            with open(file) as input_handle, open(f"clean_seq_{idx_num}.fasta", "a") as output_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    if list_of_frameshift.lower() in record.id.lower():
                        SeqIO.write(record, output_handle, "fasta")
        return True

    try:
        if os.path.exists("your_final_result.fasta"):
            os.remove("your_final_result.fasta")
    except Exception as e:
        log_step(temp_file_name, f"Error while trying to remove your_final_result.fasta: {e}")

    seq_dict_refprots = SeqIO.to_dict(SeqIO.parse(refprots_file, "fasta"))
    seq_dict_refprots_to_list = getList(seq_dict_refprots)
    seq_dict_altprots = SeqIO.to_dict(SeqIO.parse(altprots_file, "fasta"))
    seq_dict_altprots_to_list = getList(seq_dict_refprots)
    seq_dict_mrna = SeqIO.to_dict(SeqIO.parse(transcript_file, "fasta"))
    seq_dict_mrna_to_list = getList(seq_dict_mrna)

    log_file = ""

    seq_dict_refprots = dict(chain.from_iterable(d.items() for d in (seq_dict_refprots, seq_dict_altprots)))
    with open(input_file) as f:
        list_of_altProts = [line.strip() for line in f]
    stepnumber = 0
    skip_number_nts = 30
    for altprot in list_of_altProts:
        try:
            isthereanyerror = 0
            category_type = ""
            altprots_transcript = "_".join(altprot.split("_")[:-3])
            altprot_length_in_nt = altprot.split("_")[-1]
            altprot_start = int(altprot.split("_")[-2].split("-")[0])
            altprot_end = int(altprot.split("_")[-2].split("-")[1])
            refprot_transcript_list = list(
                filter(lambda x: x.startswith(altprots_transcript), seq_dict_refprots_to_list))
            altprot_altprot_overlap_list = list(filter(lambda x: x.startswith(altprots_transcript), list_of_altProts))
            altprot_altprot_overlap_list = list(filter(lambda a: a != str(altprot), altprot_altprot_overlap_list))
            refprot_transcript_list = refprot_transcript_list + altprot_altprot_overlap_list
            for refprot_transcript in refprot_transcript_list:
                try:
                    refprot_transcript = str(refprot_transcript)
                    refprot_start = int(refprot_transcript.split("_")[-2].split("-")[0])
                    refprot_end = int(refprot_transcript.split("_")[-2].split("-")[1])
                    altprot_in_aa_code = str(seq_dict_altprots[altprot].seq)
                    mrna_transcript_in_rna_code = " " + str(seq_dict_mrna[
                                                                altprots_transcript].seq)  # empty space in the beginning starts from the 1st index
                    altprot_frame = altprot.split("_")[-3]
                    ref_frame = refprot_transcript.split("_")[-3]
                    AS1 = altprot_start  # altprot 5' end
                    AS2 = altprot_end  # altprot 3' end
                    RS1 = refprot_start  # ref 5' end
                    RS2 = refprot_end  # ref 3' end
                    altprot_range = range(altprot_start, altprot_end);
                    refprot_range = range(refprot_start, refprot_end)
                    if (int(AS1) >= int(RS1)) and (int(AS1) >= int(RS2)) and (int(AS2 >= RS1)) and (
                            int(AS2) >= int(RS2)):
                        isthereanyerror += 1
                        category_type = "3' UTR side (No overlapped)"
                        log_step(temp_file_name, f"Category type = 3' UTR side (No overlapped)")
                        if int(altprot_start) - int(refprot_end) <= 10:
                            altprot_sequence_first_20_aa = str(seq_dict_altprots[altprot].seq)[:20]
                            refprot_sequence_second_20_aa = str(seq_dict_refprots[refprot_transcript].seq)[-20:]
                            append_aa_seqs_for_no_overlapped(firstsequence=refprot_sequence_second_20_aa,
                                                             secondsequence=altprot_sequence_first_20_aa,
                                                             altprotname=altprot, refprotname=refprot_transcript,
                                                             width=int(int(altprot_start) - int(refprot_end)),
                                                             fiveorthreeprime="3'_UTR_")

                    if (int(AS1) >= int(RS1)) and (int(AS1) <= int(RS2)) and (int(AS2) >= int(RS1)) and (
                            int(AS2) >= int(RS2)):
                        isthereanyerror += 1
                        category_type = "3' UTR side overlap with Reference ORF"
                        log_step(temp_file_name, f"category_type = 3' UTR side (No overlapped)")
                        overlapped_region = sorted(list(set(altprot_range) & set(refprot_range)))
                        if len(overlapped_region) > 0:

                            coordinate_acc_to_altprot = int(
                                int(closestMultiple(overlapped_region[0])) + int(altprot_frame[0])) - 3
                            coordinate_acc_to_refprot = int(int(closestMultiple(overlapped_region[0])) + int(
                                ref_frame[0])) - 3  # You may add -3 for fix borders.
                            refprot_non_range = range(refprot_end, altprot_end)
                            overlapped_region_non_range = sorted(list(set(altprot_range) & set(refprot_non_range)))
                            if int(overlapped_region[0]) - int(refprot_start) >= 30:
                                go_back = 30
                            else:
                                go_back = int(closestMultiple(int(overlapped_region[0]) - int(refprot_start)))
                            if int(overlapped_region[-1]) - int(overlapped_region[0]) >= 90:
                                for_loop_number = 90
                            else:
                                for_loop_number = int(overlapped_region[-1]) - int(overlapped_region[0]) + 6
                            if int(altprot_end) - int(altprot_start) >= 90:
                                length_second_part_raw = 90
                            else:
                                length_second_part_raw = int(altprot_end) - int(altprot_start)
                            if int(coordinate_acc_to_refprot - go_back) < 0:
                                go_back = go_back - 3
                            sequence = mrna_transcript_in_rna_code[(coordinate_acc_to_refprot - go_back):]
                            translate_all_frameshifted(sequence=sequence, altprot_name=altprot,
                                                       refprot_name=refprot_transcript, skip_value=int(go_back),
                                                       for_loop_number=int(for_loop_number),
                                                       type_acc_to_position="3'UTR_overlapped_with_CDS")
                            aa_seq_should_be_included = str(
                                translate(mrna_transcript_in_rna_code[coordinate_acc_to_altprot + 6:]))[:10]
                            list_of_frameshifts = determine_frame_from_file(common_seq=aa_seq_should_be_included)
                            if len(list_of_frameshifts) != 2:
                                log_file = log_file + str("Correct frameshifts could not determined for ") + str(
                                    altprot) + ", " + str(refprot_transcript) + str(". Detected frameshifts: ") + str(
                                    list_of_frameshifts) + str("\n")
                            get_only_frame(list_of_frameshifts=list_of_frameshifts)

                    if (int(AS1) >= int(RS1)) and (int(AS1) <= int(RS2)) and (int(AS2) >= int(RS1)) and (
                            int(AS2) <= int(RS2)):
                        isthereanyerror += 0.5
                        category_type = "Within (Embedded) Reference ORF"
                        log_step(temp_file_name, f"Category type = 'Within (Embedded) Reference ORF'")
                        overlapped_region = sorted(list(set(altprot_range) & set(refprot_range)))
                        if len(overlapped_region) > 0:

                            coordinate_acc_to_altprot = int(
                                int(closestMultiple(overlapped_region[0])) + int(altprot_frame[0])) - 3
                            coordinate_acc_to_refprot = int(int(closestMultiple(overlapped_region[0])) + int(
                                ref_frame[0])) - 3  # You may add -3 for fix borders.
                            refprot_non_range = range(refprot_start, altprot_start)
                            overlapped_region_non_range = sorted(list(set(refprot_range) & set(refprot_non_range)))
                            if int(overlapped_region[0]) - int(refprot_start) >= 30:
                                go_back = 30
                            else:
                                go_back = int(closestMultiple(int(overlapped_region[0]) - int(refprot_start)))
                            if int(overlapped_region[-1]) - int(overlapped_region[0]) >= 90:
                                for_loop_number = 90 - 30 + 3
                            else:
                                for_loop_number = int(overlapped_region[-1]) - int(overlapped_region[0]) - 30 + 3
                            if int(altprot_end) - int(altprot_start) >= 90:
                                length_second_part_raw = 90
                            else:
                                length_second_part_raw = int(altprot_end) - int(altprot_start)
                            if int(coordinate_acc_to_refprot - go_back) < 0:
                                go_back = go_back - 3
                            sequence = mrna_transcript_in_rna_code[(coordinate_acc_to_refprot - go_back):]
                            translate_all_frameshifted(sequence=sequence, altprot_name=altprot,
                                                       refprot_name=refprot_transcript, skip_value=int(go_back),
                                                       for_loop_number=int(for_loop_number),
                                                       type_acc_to_position="Within_5'_of_altprot")
                            aa_seq_should_be_included = str(
                                translate(mrna_transcript_in_rna_code[coordinate_acc_to_altprot + 6:]))[:10]
                            list_of_frameshifts = determine_frame_from_file(common_seq=aa_seq_should_be_included)
                            if len(list_of_frameshifts) != 2:
                                log_file = log_file + str("Correct frameshifts could not determined for ") + str(
                                    altprot) + ", " + str(refprot_transcript) + str(". Detected frameshifts: ") + str(
                                    list_of_frameshifts) + str("\n")
                            get_only_frame(list_of_frameshifts=list_of_frameshifts)

                    if (int(AS1) >= int(RS1)) and (int(AS1) <= int(RS2)) and (int(AS2) >= int(RS1)) and (
                            int(AS2) <= int(RS2)):
                        isthereanyerror += 0.5
                        category_type = "Within (Embedded) Reference ORF"
                        log_step(temp_file_name, f"Category type = 'Within (Embedded) Reference ORF'")
                        overlapped_region = sorted(list(set(altprot_range) & set(refprot_range)))
                        if len(overlapped_region) > 0:

                            coordinate_acc_to_altprot = int(
                                int(closestMultiple(overlapped_region[-1])) + int(altprot_frame[0])) - 3
                            coordinate_acc_to_refprot = int(int(closestMultiple(overlapped_region[0])) + int(
                                ref_frame[0])) - 3  # You may add -3 for fix borders.
                            refprot_non_range = range(altprot_end, refprot_end)
                            overlapped_region_non_range = sorted(list(set(refprot_range) & set(refprot_non_range)))
                            if int(overlapped_region[-1]) - int(overlapped_region[0]) >= 90:
                                go_back = 90
                            else:
                                go_back = int(closestMultiple(int(overlapped_region[-1]) - int(overlapped_region[0])))
                            if int(overlapped_region[-1]) - int(overlapped_region[0]) >= 90:
                                for_loop_number = 90 - 30 + 3
                            else:
                                for_loop_number = int(overlapped_region[-1]) - int(overlapped_region[0]) - 30 + 3
                            if int(refprot_end) - int(altprot_end) >= 90:
                                length_second_part_raw = 90
                            else:
                                length_second_part_raw = int(refprot_end) - int(altprot_end)
                            if int(coordinate_acc_to_altprot - go_back) < 0:
                                go_back = go_back - 3
                            sequence = mrna_transcript_in_rna_code[(coordinate_acc_to_altprot - go_back):]
                            translate_all_frameshifted(sequence=sequence, altprot_name=altprot,
                                                       refprot_name=refprot_transcript,
                                                       skip_value=int(go_back) - 60,
                                                       for_loop_number=int(for_loop_number),
                                                       type_acc_to_position="Within_3'_of_altprot")
                            aa_seq_should_be_included = str(translate(mrna_transcript_in_rna_code[
                                                                      coordinate_acc_to_refprot + (
                                                                              altprot_end - altprot_start + 1) - 57:]))[
                                                        :10]
                            list_of_frameshifts = determine_frame_from_file(common_seq=aa_seq_should_be_included)
                            if len(list_of_frameshifts) != 2:
                                log_file = log_file + str("Correct frameshifts could not determined for ") + str(
                                    altprot) + ", " + str(refprot_transcript) + str(". Detected frameshifts: ") + str(
                                    list_of_frameshifts) + str("\n")
                            get_only_frame(list_of_frameshifts=list_of_frameshifts)

                    if (int(AS1) <= int(RS1)) and (int(AS1) <= int(RS2)) and (int(AS2) >= int(RS1)) and (
                            int(AS2) <= int(RS2)):
                        isthereanyerror += 1
                        category_type = "5' UTR side overlap with Reference ORF"
                        log_step(temp_file_name, f"Category type = 5' UTR side overlap with Reference ORF")
                        overlapped_region = sorted(
                            list(set(range(int(altprot_range[0]), int(altprot_range[-1] + 4))) & set(refprot_range)))
                        if len(overlapped_region) > 0:

                            coordinate_acc_to_altprot = int(
                                int(closestMultiple(overlapped_region[0])) + int(altprot_frame[0])) - 3
                            coordinate_acc_to_refprot = int(int(closestMultiple(overlapped_region[0])) + int(
                                ref_frame[0])) - 3  # You may add -3 for fix borders.
                            refprot_non_range = range(altprot_start, refprot_start)
                            overlapped_region_non_range = sorted(list(set(altprot_range) & set(refprot_non_range)))
                            if int(overlapped_region_non_range[-1]) - int(overlapped_region_non_range[0]) >= 30:
                                go_back = 30
                            else:
                                go_back = int(closestMultiple(
                                    int(overlapped_region_non_range[-1]) - int(overlapped_region_non_range[0])))
                            if int(overlapped_region[-1]) - int(overlapped_region[0]) >= 90:
                                for_loop_number = 90
                            else:
                                for_loop_number = int(overlapped_region[-1]) - int(overlapped_region[0])
                            if int(refprot_end) - int(refprot_start) >= 90:
                                length_second_part_raw = 90
                            else:
                                length_second_part_raw = int(refprot_end) - int(refprot_start)
                            if int(coordinate_acc_to_altprot - go_back) < 0:
                                go_back = go_back - 3
                            sequence = mrna_transcript_in_rna_code[(coordinate_acc_to_altprot - go_back):]
                            translate_all_frameshifted(sequence=sequence, altprot_name=altprot,
                                                       refprot_name=refprot_transcript, skip_value=int(go_back),
                                                       for_loop_number=int(for_loop_number),
                                                       type_acc_to_position="5'UTR_overlapped_with_CDS")
                            aa_seq_should_be_included = str(
                                translate(mrna_transcript_in_rna_code[coordinate_acc_to_refprot + 6:]))[:10]
                            list_of_frameshifts = determine_frame_from_file(common_seq=aa_seq_should_be_included)
                            if len(list_of_frameshifts) != 2:
                                log_file = log_file + str("Correct frameshifts could not determined for ") + str(
                                    altprot) + ", " + str(refprot_transcript) + str(". Detected frameshifts: ") + str(
                                    list_of_frameshifts) + str("\n")
                            get_only_frame(list_of_frameshifts=list_of_frameshifts)

                    if (int(AS1) <= int(RS1)) and (int(AS1) <= int(RS2)) and (int(AS2) <= int(RS1)) and (
                            int(AS2) <= int(RS2)):
                        isthereanyerror += 1
                        category_type = "5' UTR side (No overlapped)"
                        log_step(temp_file_name, f"Category type = 5' UTR side (No overlapped)")
                        if int(refprot_start) - int(altprot_end) <= 10:
                            altprot_sequence_first_20_aa = str(seq_dict_altprots[altprot].seq)[-20:]
                            refprot_sequence_second_20_aa = str(seq_dict_refprots[refprot_transcript].seq)[:20]
                            append_aa_seqs_for_no_overlapped(firstsequence=altprot_sequence_first_20_aa,
                                                             secondsequence=refprot_sequence_second_20_aa,
                                                             altprotname=altprot, refprotname=refprot_transcript,
                                                             width=int(int(refprot_start) - int(altprot_end)),
                                                             fiveorthreeprime="5'_UTR_")

                    if (int(AS1) <= int(RS1)) and (int(AS1) <= int(RS2)) and (int(AS2) >= int(RS1)) and (
                            int(AS2) >= int(RS2)):
                        isthereanyerror += 0.5
                        category_type = "Spanned"
                        log_step(temp_file_name, f"Category type = Spanned - 5' UTR side of altProt")
                        overlapped_region = sorted(list(set(altprot_range) & set(refprot_range)))
                        if len(overlapped_region) > 0:

                            coordinate_acc_to_altprot = int(
                                int(closestMultiple(overlapped_region[0])) + int(altprot_frame[0])) - 3
                            coordinate_acc_to_refprot = int(int(closestMultiple(overlapped_region[0])) + int(
                                ref_frame[0])) - 3  # You may add -3 for fix borders.
                            refprot_non_range = range(altprot_start, refprot_start)
                            overlapped_region_non_range = sorted(list(set(altprot_range) & set(refprot_non_range)))
                            if int(overlapped_region_non_range[-1]) - int(overlapped_region_non_range[0]) >= 30:
                                go_back = 30
                            else:
                                go_back = int(closestMultiple(
                                    int(overlapped_region_non_range[-1]) - int(overlapped_region_non_range[0])))
                            if int(overlapped_region[-1]) - int(overlapped_region[0]) >= 90:
                                for_loop_number = 90
                            else:
                                for_loop_number = int(overlapped_region[-1]) - int(overlapped_region[0])
                            if int(refprot_end) - int(refprot_start) >= 90:
                                length_second_part_raw = 90
                            else:
                                length_second_part_raw = int(refprot_end) - int(refprot_start)
                            if int(coordinate_acc_to_altprot - go_back) < 0:
                                go_back = go_back - 3
                            sequence = mrna_transcript_in_rna_code[(coordinate_acc_to_altprot - go_back):]
                            translate_all_frameshifted(sequence=sequence, altprot_name=altprot,
                                                       refprot_name=refprot_transcript, skip_value=int(go_back),
                                                       for_loop_number=int(for_loop_number),
                                                       type_acc_to_position="Spanned_5'_of_altprot")
                            aa_seq_should_be_included = str(
                                translate(mrna_transcript_in_rna_code[coordinate_acc_to_refprot + 6:]))[:10]
                            list_of_frameshifts = determine_frame_from_file(common_seq=aa_seq_should_be_included)
                            if len(list_of_frameshifts) != 2:
                                log_file = log_file + str("Correct frameshifts could not determined for ") + str(
                                    altprot) + ", " + str(refprot_transcript) + str(". Detected frameshifts: ") + str(
                                    list_of_frameshifts) + str("\n")
                            get_only_frame(list_of_frameshifts=list_of_frameshifts)

                    if (int(AS1) <= int(RS1)) and (int(AS1) <= int(RS2)) and (int(AS2) >= int(RS1)) and (
                            int(AS2) >= int(RS2)):
                        isthereanyerror += 0.5
                        category_type = "Spanned"
                        log_step(temp_file_name, f"Category type = Spanned - 3' UTR side of altProt")
                        overlapped_region = sorted(list(set(altprot_range) & set(refprot_range)))
                        if len(overlapped_region) > 0:

                            coordinate_acc_to_altprot = int(
                                int(closestMultiple(overlapped_region[0])) + int(altprot_frame[0])) - 3
                            coordinate_acc_to_refprot = int(int(closestMultiple(overlapped_region[-1])) + int(
                                ref_frame[0])) - 3  # You may add -3 for fix borders.
                            refprot_non_range = range(refprot_end, altprot_end)
                            overlapped_region_non_range = sorted(list(set(altprot_range) & set(refprot_non_range)))
                            if int(overlapped_region[-1]) - int(refprot_start) >= 90:
                                go_back = 90
                            else:
                                go_back = int(closestMultiple(int(overlapped_region[-1]) - int(refprot_start)))
                            if int(overlapped_region[-1]) - int(overlapped_region[0]) >= 90:
                                for_loop_number = 90
                            else:
                                for_loop_number = int(overlapped_region[-1]) - int(overlapped_region[0])
                            if int(altprot_end) - int(refprot_end) >= 30:
                                length_second_part_raw = 90
                            else:
                                length_second_part_raw = int(altprot_end) - int(refprot_end)
                            if int(coordinate_acc_to_refprot - go_back) < 0:
                                go_back = go_back - 3
                            sequence = mrna_transcript_in_rna_code[(coordinate_acc_to_refprot - go_back):]
                            translate_all_frameshifted(sequence=sequence, altprot_name=altprot,
                                                       refprot_name=refprot_transcript,
                                                       skip_value=int(go_back) - 60,
                                                       for_loop_number=int(for_loop_number),
                                                       type_acc_to_position="Spanned_3'_of_altprot")
                            aa_seq_should_be_included = str(translate(mrna_transcript_in_rna_code[
                                                                      coordinate_acc_to_altprot + (
                                                                              refprot_end - refprot_start + 1) - 27:]))[
                                                        :10]
                            list_of_frameshifts = determine_frame_from_file(common_seq=aa_seq_should_be_included)
                            if len(list_of_frameshifts) != 2:
                                log_file = log_file + str("Correct frameshifts could not determined for ") + str(
                                    altprot) + ", " + str(refprot_transcript) + str(". Detected frameshifts: ") + str(
                                    list_of_frameshifts) + str("\n")
                            get_only_frame(list_of_frameshifts=list_of_frameshifts)

                    if isthereanyerror == 1:
                        continue
                    else:
                        log_step(temp_file_name,
                                 f"It may be an error or your altprot belongs to more than one category. Check these entries manually."
                                 f"altprot name is {altprot}")

                except Exception as e:
                    log_file = log_file + str(e) + "Error for altprot and refprot are: " + str(altprot) + ", " + str(
                        refprot_transcript)

        except Exception as e:
            refprot_transcript = None
            log_file = log_file + str(e) + "Error for altprot and refprot are: " + str(altprot) + ", " + str(
                refprot_transcript)

        try:
            if os.path.exists(f"clean_seq_{idx_num}.fasta"):
                with open(f"clean_seq_{idx_num}.fasta", "r") as src, open(f"your_final_result_{idx_num}.fasta",
                                                                          "a") as dst:
                    dst.write(src.read())

        except:
            log_step(temp_file_name, f"clean_seq_{idx_num}.fasta file is not found for {category_type}")

        try:
            if os.path.exists(f"clean_seq_{idx_num}.fasta"):
                os.remove(f"clean_seq_{idx_num}.fasta")

        except Exception as e:
            log_step(temp_file_name, f"Error while trying to remove clean_seq_{idx_num}.fasta: {e}")

        except Exception as e:
            log_step(temp_file_name, f"Error while trying to remove clean_seq.fasta: {e}")

        try:
            if os.path.exists(f"temp_file_for_seq_store_{idx_num}.txt"):
                os.remove(f"temp_file_for_seq_store_{idx_num}.txt")
        except Exception as e:
            log_step(temp_file_name, f"Error while trying to remove temp_file_for_seq_store_{idx_num}.txt: {e}")

        try:
            if os.path.exists(f"temp_file_for_seq_store_temp_{idx_num}.txt"):
                os.remove(f"temp_file_for_seq_store_temp_{idx_num}.txt")
        except Exception as e:
            log_step(temp_file_name, f"Error while trying to remove temp_file_for_seq_store_temp_{idx_num}.txt: {e}")

    input_file = f"your_final_result_{idx_num}.fasta"
    output_file = f"your_final_result_before_clean_{idx_num}.fasta"

    i = 1
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                line = line.strip() + f"_{i}\n"
                i += 1
            outfile.write(line)

    final_results = SeqIO.to_dict(SeqIO.parse(f"your_final_result_before_clean_{idx_num}.fasta", "fasta"))
    final_results_to_list = getList(final_results)
    to_file_for_final_result = ""
    for final_result in final_results_to_list:

        sequence_for_final_result = str(final_results[final_result].seq)
        if "*" not in sequence_for_final_result[10:]:
            if "*" in sequence_for_final_result[:10]:
                to_file_for_final_result = to_file_for_final_result + str(
                    ">" + str(final_results[final_result].id) + "\n" + str(
                        sequence_for_final_result.split("*")[1]) + "\n")
            else:
                to_file_for_final_result = to_file_for_final_result + str(
                    ">" + str(final_results[final_result].id) + "\n" + str(sequence_for_final_result) + "\n")
        else:
            to_file_for_final_result = to_file_for_final_result + str(
                ">" + str(final_results[final_result].id) + "\n" + str(sequence_for_final_result) + "\n")

    clean_final_result = open(f'clean_final_result_{idx_num}.fasta', 'w+')
    clean_final_result.write(to_file_for_final_result)
    clean_final_result.close()

    upload_results = SeqIO.to_dict(SeqIO.parse(f"clean_final_result_{idx_num}.fasta", "fasta"))
    upload_results_to_list = getList(upload_results)
    to_file_for_upload_result = ""

    for upload_result in upload_results_to_list:

        sequence_for_upload_result = str(upload_results[upload_result].seq)
        if "*" in sequence_for_upload_result:
            to_file_for_upload_result = to_file_for_upload_result + str(
                ">" + str(upload_results[upload_result].id) + "\n" + str(
                    sequence_for_upload_result.split("*")[0]) + "\n")
        else:
            to_file_for_upload_result = to_file_for_upload_result + str(
                ">" + str(upload_results[upload_result].id) + "\n" + str(sequence_for_upload_result) + "\n")

    clean_upload_result = open(f'clean_upload_result_{idx_num}.fasta', 'w+')
    clean_upload_result.write(to_file_for_upload_result)
    clean_upload_result.close()

    results_directory = ".mosaic_prot_results"

    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    input_file = f"clean_upload_result_{idx_num}.fasta"
    output_file = os.path.join(results_directory, f"requested_file_{idx_num}.fasta")
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if len(record.seq) >= 11:
                SeqIO.write(record, outfile, "fasta")


    try:
        if os.path.exists(f"clean_upload_result_{idx_num}.fasta"):
            os.remove(f"clean_upload_result_{idx_num}.fasta")
        else:
            log_step(temp_file_name, "clean_upload_result.fasta does not exist, so it cannot be removed.")
            print("clean_upload_result.fasta does not exist, so it cannot be removed.")
    except Exception as e:
        log_step(temp_file_name, f"Error while trying to remove clean_upload_result.fasta: {e}")
        print(f"Error while trying to remove clean_upload_result_{idx_num}.fasta: {e}")

    try:
        if os.path.exists(f"your_final_result_{idx_num}.fasta"):
            os.remove(f"your_final_result_{idx_num}.fasta")

    except Exception as e:
        log_step(temp_file_name, f"Error while trying to remove your_final_result_{idx_num}.fasta: {e}")
    try:
        if os.path.exists(f"your_final_result_before_clean_{idx_num}.fasta"):
            os.remove(f"your_final_result_before_clean_{idx_num}.fasta")
        else:
            log_step(temp_file_name,
                     f"your_final_result_before_clean_{idx_num}.fasta does not exist, so it cannot be removed.")
    except Exception as e:
        log_step(temp_file_name, f"Error while trying to remove your_final_result_before_clean_{idx_num}.fasta: {e}")
    try:
        if os.path.exists(f"clean_final_result_{idx_num}.fasta"):
            os.remove(f"clean_final_result_{idx_num}.fasta")
        else:
            log_step(temp_file_name, f"clean_final_result_{idx_num}.fasta does not exist, so it cannot be removed.")
    except Exception as e:
        log_step(temp_file_name, f"Error while trying to remove clean_final_result_{idx_num}.fasta: {e}")

    temp_path_ = os.path.join(".mosaic_prot_results", f"temp_sequence_archieve_{idx_num}.fasta")

    with open(temp_path_, "w") as outfile:
        for file_path in [altprots_file, refprots_file]:
            with open(file_path, "r") as infile:
                outfile.write(infile.read())

    results_directory = ".mosaic_prot_results"

    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    temp_path_ = os.path.join(results_directory, f"temp_sequence_archieve_{idx_num}.fasta")

    results_directory = ".mosaic_prot_results"

    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    temp_path_ = os.path.join(results_directory, f"temp_sequence_archieve_{idx_num}.fasta")
    backup_path = temp_path_ + ".withoutnumbers"  # to mimic Perl's backup

    shutil.copyfile(temp_path_, backup_path)

    updated_records = []
    counter = 1
    with open(temp_path_, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            record.id += f"|{counter:08d}"  # zero-padded to 8 digits
            record.description = record.id  # update full header line
            updated_records.append(record)
            counter += 1

    with open(temp_path_, "w") as outfile:
        SeqIO.write(updated_records, outfile, "fasta")

    sequence_archive = SeqIO.to_dict(SeqIO.parse(temp_path_, "fasta"))

    sequence_archive = SeqIO.to_dict(SeqIO.parse(temp_path_, "fasta"))
    results_directory = ".mosaic_prot_results"

    requested_file_path = os.path.join(results_directory, f"requested_file_{idx_num}.fasta")

    requested_file = SeqIO.to_dict(SeqIO.parse(requested_file_path, "fasta"))

    requested_file_to_list = getList(requested_file)

    to_file_for_requested_file_filtered = ""

    for requested_item in requested_file_to_list:
        first_item_id = join_string(str(requested_item).split("_")[0:5])
        second_item_id = join_string(str(requested_item).split("_")[5:10])
        first_item_seq = str(
            sequence_archive[str(list(filter(lambda x: x.startswith(first_item_id), sequence_archive))[0])].seq)
        second_item_seq = str(
            sequence_archive[str(list(filter(lambda x: x.startswith(second_item_id), sequence_archive))[0])].seq)
        requested_item_seq = str(requested_file[requested_item].seq)
        if requested_item_seq not in first_item_seq:
            if requested_item_seq not in second_item_seq:
                to_file_for_requested_file_filtered = to_file_for_requested_file_filtered + str(
                    ">" + requested_item + "\n" + str(requested_file[requested_item].seq) + "\n")

    results_directory = ".mosaic_prot_results"

    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    file_path = os.path.join(results_directory, f'requested_file_filtered{idx_num}.fasta')

    with open(file_path, 'w+') as file_for_requested_file_filtered:
        file_for_requested_file_filtered.write(to_file_for_requested_file_filtered)
    updated_records = []
    with open(file_path, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            record.id += "_CHIMERIC"
            record.description = record.id
            updated_records.append(record)

    with open(file_path, "w") as outfile:
        SeqIO.write(updated_records, outfile, "fasta")



def worker_process(args):
    """Worker function to process each chunk in parallel."""
    chunk, refprots_file, altprots_file, transcript_file, temp_file_name, idx = args
    generate_mosaic_proteins(
        input_file=chunk,
        refprots_file=refprots_file,
        altprots_file=altprots_file,
        transcript_file=transcript_file,
        temp_file_name=temp_file_name,
        idx_num=idx
    )


def split_file(input_file, output_dir, num_chunks):

    with open(input_file, 'r') as f:
        lines = f.readlines()

    lines_per_chunk = len(lines) // num_chunks
    remainder = len(lines) % num_chunks

    chunks = []
    start = 0
    for i in range(num_chunks):
        end = start + lines_per_chunk + (1 if i < remainder else 0)
        chunks.append(lines[start:end])
        start = end

    for i, chunk in enumerate(chunks):
        chunk_file = os.path.join(output_dir, f'chunk_{i}.txt')
        with open(chunk_file, 'w') as f:
            f.writelines(chunk)


def simulate_chimeric_proteins(input_file, refprots_file, altprots_file, transcript_file, processor_num, temp_file_name):
    print("Generating chimeric proteins...")

    with open(input_file, 'r') as f:
        lines = f.read().splitlines()

    transcript_ids = []
    for line in lines:
        parts = line.split("_")
        if len(parts) >= 2:
            transcript_ids.append("_".join(parts[:-3]))
        else:
            print(f"Skipping malformed line: {line}")

    results_directory = ".mosaic_prot_results"
    os.makedirs(results_directory, exist_ok=True)

    duplicate_ids = [item for item, count in Counter(transcript_ids).items() if count > 1]
    always_be_together = [line for line in lines if any(dup in line for dup in duplicate_ids)]
    should_be_splitted = set(lines) - set(always_be_together)

    list_of_altProts_path = os.path.join(results_directory, "list_of_altProts.txt")
    always_be_together_path = os.path.join(results_directory, "always_be_together.txt")
    with open(list_of_altProts_path, "w") as f:
        f.write("\n".join(should_be_splitted))
    with open(always_be_together_path, "w") as f:
        f.write("\n".join(always_be_together))

    list_of_chunks = []
    if processor_num == 1:
        combined_chunk = os.path.join(results_directory, "chunk_combined.txt")
        with open(combined_chunk, "w") as f:
            f.write("\n".join(always_be_together + list(should_be_splitted)))
        list_of_chunks = [combined_chunk]
    else:
        list_of_chunks.append(always_be_together_path)

        split_file(list_of_altProts_path, results_directory, processor_num - 1)

        list_of_chunks.extend(glob.glob(os.path.join(results_directory, 'chunk_*.txt')))

    with Pool(processor_num) as pool:
        args = [
            (chunk, refprots_file, altprots_file, transcript_file, temp_file_name, idx)
            for idx, chunk in enumerate(list_of_chunks)
        ]
        pool.starmap(generate_mosaic_proteins, args)
    print("Chimeric protein generation completed.")


def log_step(temp_file_name, message, overwrite=False):
    timestamp = time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime())
    mode = "w" if overwrite else "a"  # Overwrite or append
    with open(temp_file_name, mode) as temp_file:
        temp_file.write(f"{message}\n")


def format_time(seconds):
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    return f"{hours:02}:{minutes:02}:{secs:02}"


def clear_log(temp_file_name):
    with open(temp_file_name, "w") as temp_file:
        temp_file.write("")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        description="MosaicProt CLI: Detect ORFs, Find Alternative ORFs, Generate Mosaic Proteins",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_detect_orfs = subparsers.add_parser(
        'detect_ORFs',
        help="Detect ORFs from a transcriptome file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_detect_orfs.add_argument('--transcriptome_file', required=True, help="Input transcriptome file")
    parser_detect_orfs.add_argument('--threshold', type=int, default=30, help="ORF length threshold")
    parser_detect_orfs.add_argument('--output_file_type', choices=['fasta', 'xml'], default="fasta",
                                    help="Output file format")

    parser_find_alt_orfs = subparsers.add_parser(
        'separate_ORFs',
        help="Separate alternative and reference ORFs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_find_alt_orfs.add_argument('--ORFeome_file', required=True, help="Input ORFeome file")
    parser_find_alt_orfs.add_argument('--reference_proteome_file', required=True, help="Input reference proteome file")
    parser_find_alt_orfs.add_argument('--output_alt_file', type=str, default="altProts.fasta",
                                      help="Output file for alternative ORFs")
    parser_find_alt_orfs.add_argument('--output_ref_file', type=str, default="refProts.fasta",
                                      help="Output file for reference ORFs")

    parser_chimeric = subparsers.add_parser(
        'simulate_chimeric_proteins',
        help="Simulate chimeric proteins",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_chimeric.add_argument('--candidate_altProt_list', required=True,
                                 help="Candidate altProts for the simulation")
    parser_chimeric.add_argument('--refProts_file', default="refProts.fasta",
                                 help="Reference proteins file")
    parser_chimeric.add_argument('--altProts_file', default="altProts.fasta",
                                 help="Alternative proteins file")
    parser_chimeric.add_argument('--transcriptome_file', required=True, help="Transcriptome file")
    parser_chimeric.add_argument('--processor_num', type=int, default=1, help="Number of processors to use")
    parser_chimeric.add_argument('--repetition', choices=['keep_first', 'drop_all', 'keep_all'], default='keep_first',
                                 help="How to handle duplicates in the final simulated chimeric proteins FASTA file.")

    try:
        args = parser.parse_args()

        log_file_name = os.path.join(script_dir, f"logfile_{args.command}.temp")

        log_step(log_file_name, f"Starting MosaicProt CLI: {args.command}...", overwrite=True)

        start_time = time.time()

        if args.command == "detect_ORFs":
            log_step(log_file_name, "Processing ORFs...")
            transcriptome_basename = os.path.splitext(args.transcriptome_file)[0]
            output_file = f"{transcriptome_basename}_min_{args.threshold}aa_ORFs.{args.output_file_type}"
            process_fasta_output(args.transcriptome_file, args.threshold, output_file, log_file_name)
            log_step(log_file_name, f"ORFs processed and saved to {output_file}")

        elif args.command == "separate_ORFs":
            log_step(log_file_name, "Finding alternative ORFs...")
            run_alt_orf(args.ORFeome_file, args.reference_proteome_file, args.output_alt_file, args.output_ref_file,
                        log_file_name)
            log_step(log_file_name,
                     f"Alternative ORFs saved to {args.output_alt_file}, reference saved to {args.output_ref_file}")

        elif args.command == "simulate_chimeric_proteins":
            log_step(log_file_name, "Forming chimeric proteins...")
            simulate_chimeric_proteins(
                args.candidate_altProt_list,
                args.refProts_file,
                args.altProts_file,
                args.transcriptome_file,
                args.processor_num,
                log_file_name
            )

            results_directory = ".mosaic_prot_results"

            if not os.path.exists(results_directory):
                os.makedirs(results_directory)

            files_to_combine = glob.glob(os.path.join(results_directory, "requested_file_filtered*.fasta"))

            combined_file = os.path.join(results_directory, "simulated_chimeric_proteins.fasta")

            with open(combined_file, "w") as outfile:
                for file in files_to_combine:
                    with open(file, "r") as infile:
                        outfile.write(infile.read())

            current_directory = os.getcwd()
            final_combined_file = os.path.join(current_directory, "simulated_chimeric_proteins.fasta")
            shutil.move(combined_file, final_combined_file)

            def read_fasta(file_path):
                records = []
                with open(file_path, 'r') as f:
                    header = ''
                    sequence = []
                    for line in f:
                        line = line.strip()
                        if line.startswith('>'):
                            if header:
                                records.append((header, ''.join(sequence)))
                            header = line[1:]
                            sequence = []
                        else:
                            sequence.append(line)
                    if header:
                        records.append((header, ''.join(sequence)))
                return records

            def extract_tags_and_iteration(header):
                iteration_match = re.search(r'iteration_(\d+)', header)
                iteration = int(iteration_match.group(1)) if iteration_match else None

                header_part = header.split('_iteration_')[0]
                tag_pattern = re.compile(r'(Mtrun\w+(?:_\w+)+)')
                tags = tag_pattern.findall(header_part)

                return tags, iteration

            def find_duplicate_groups(records):
                grouped = defaultdict(list)
                for idx, (header, seq) in enumerate(records):
                    tags, iteration = extract_tags_and_iteration(header)
                    grouped_key = (seq, iteration)
                    grouped[grouped_key].append((header, tags, idx))  # Include index for ordering

                duplicate_groups = []
                for key in grouped:
                    headers_tags_indices = grouped[key]
                    adjacency = defaultdict(list)
                    for i in range(len(headers_tags_indices)):
                        h1, tags1, idx1 = headers_tags_indices[i]
                        for j in range(i + 1, len(headers_tags_indices)):
                            h2, tags2, idx2 = headers_tags_indices[j]
                            if set(tags1) & set(tags2):
                                adjacency[h1].append(h2)
                                adjacency[h2].append(h1)

                    visited = set()
                    for h, _, _ in headers_tags_indices:
                        if h not in visited:
                            stack = [h]
                            component = []
                            while stack:
                                node = stack.pop()
                                if node not in visited:
                                    visited.add(node)
                                    component.append(node)
                                    stack.extend(adjacency[node])
                            if len(component) > 1:
                                duplicate_groups.append(component)
                return duplicate_groups

            def handle_duplicates(records, duplicate_groups, option='keep_first'):
                all_duplicates = set()
                header_to_group = {}
                for group in duplicate_groups:
                    for header in group:
                        all_duplicates.add(header)
                        header_to_group[header] = group

                kept = []
                seen_groups = set()

                for header, seq in records:
                    if header in all_duplicates:
                        if option == 'keep_first':
                            group = header_to_group[header]
                            group_tuple = tuple(sorted(group))
                            if group_tuple not in seen_groups:
                                seen_groups.add(group_tuple)
                                kept.append((header, seq))
                        elif option == 'keep_all':
                            kept.append((header, seq))
                        elif option == 'drop_all':
                            continue
                    else:
                        kept.append((header, seq))
                return kept

            def write_fasta(records, output_path):
                with open(output_path, 'w') as f:
                    for header, seq in records:
                        f.write(f'>{header}\n{seq}\n')

            input_fasta = final_combined_file  # Replace with your input file path
            output_fasta = final_combined_file  # Replace with your desired output file path
            option = args.repetition  # Choose 'keep_first', 'keep_all', or 'drop_all'

            records = read_fasta(input_fasta)
            duplicate_groups = find_duplicate_groups(records)
            filtered_records = handle_duplicates(records, duplicate_groups, option)

            filtered_records = [
                ('_'.join(header.split('_')[:-2]), seq)
                for header, seq in filtered_records
            ]


            def count_categories(records):
                category_counts = defaultdict(int)

                for header, seq in records:
                    base_header = '_'.join(header.split('_')[:-2])

                    iteration_match = re.search(r"iteration_\d+_(.*)", base_header)
                    if iteration_match:
                        category = iteration_match.group(1)
                    else:
                        parts = base_header.split('_')
                        if "UTR" in parts:
                            idx = parts.index("UTR")
                            category = "_".join(parts[idx - 1:])  # e.g., "3'_UTR_width"
                        else:
                            category = parts[-1]

                    category_counts[category] += 1

                return category_counts

            category_counts = count_categories(filtered_records)
            print("Summary of Categories and Counts:")
            for category, count in category_counts.items():
                log_step(log_file_name, f" Number of chimeric proteins in {category} is found as: {count}")
                print(f"{category}: {count}")

            write_fasta(filtered_records, output_fasta)

            if os.path.exists(results_directory) and not os.listdir(results_directory):
                os.rmdir(results_directory)

            print(f"Combined file moved to {final_combined_file} and .mosaic_prot_results directory deleted.")

            for file in files_to_combine:
                os.remove(file)

            additional_files = [
                os.path.join(results_directory, "list_of_altProts.txt"),
                os.path.join(results_directory, "always_be_together_pre_sort.txt"),
                os.path.join(results_directory, "list_of_altProts_pre_sort.txt"),
            ]

            for file in additional_files:
                if os.path.exists(file):
                    os.remove(file)

            if os.path.exists("list_of_altProts.txt"):
                os.remove("list_of_altProts.txt")
            if os.path.exists("always_be_together_pre_sort.txt"):
                os.remove("always_be_together_pre_sort.txt")
            if os.path.exists("list_of_altProts_pre_sort.txt"):
                os.remove("list_of_altProts_pre_sort.txt")
            if os.path.exists(".mosaic_prot_results"):
                shutil.rmtree(".mosaic_prot_results")
            print(f"Combined {len(files_to_combine)} files into {combined_file} and deleted the original files.")

        elapsed_time = time.time() - start_time
        formatted_time = format_time(elapsed_time)
        log_step(log_file_name, f"Time taken for {args.command}: {formatted_time}")



    except Exception as e:
        log_step(log_file_name, f"Error occurred: {e}")
        print(f"Error: {e}")
        sys.exit(1)

    print(f"Process completed. Logs written to {log_file_name}")

def run():
    main()

if __name__ == "__main__":
    run()
