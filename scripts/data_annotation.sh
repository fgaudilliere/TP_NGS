data=/ifb/data/mydatalocal/data_tp_ngs



TransDecoder.LongOrfs -t $data/trinity_results/Trinity_RF.fasta


# Recover proteomic data from chosen reference species using Transdecoder

# Extracting the long open reading frames: 
# TransDecoder.LongOrfs -t target_transcripts.fasta
# Output file can be used to include homology searches as ORF retention criteria
# Output file: '${transcripts_file}.transdecoder_dir/longest_orfs.pep'
# Example of a BlastP search command:
# blastp -query transdecoder_dir/longest_orfs.pep  \
#    -db uniprot_sprot.fasta  -max_target_seqs 1 \
#    -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6
# The result can then be used with the following command:
# TransDecoder.Predict -t target_transcripts.fasta --retain_blastp_hits blastp.outfmt6
# Viewing ORFs on target transcripts: 
# java -jar $GENOMEVIEW/genomeview.jar transcripts.fasta transcripts.fasta.transdecoder.bed



# Use blast to detect the homologies between the transcripts and reference genes


# Write a script to create a tsv or csv annotation table