mosaicprot detect_ORFs --transcriptome_file transcriptome.fasta
mosaicprot separate_ORFs --ORFeome_file transcriptome_min_30aa_ORFs.fasta --reference_proteome_file refProt.fasta
mosaicprot simulate_chimeric_proteins --candidate_altProt_list candidate_list.txt --transcriptome_file transcriptome.fasta
