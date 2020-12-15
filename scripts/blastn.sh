data=/ifb/data/mydatalocal/data_tp_ngs


# Creating a blast reference database (in our case, human)
# dbtype: molecule type of target db, either nucl or prot
# in: input file name
# input_type: type of data in the input file, either asn1_bin', `asn1_txt', `blastdb', or `fasta'
# out: output name
# parse_seqids: Option to parse seqid for FASTA input if set, for all other input types seqids are parsed automatically
/softwares/ncbi-blast-2.10.1+/bin/makeblastdb -dbtype nucl -in $data/blastdb/human_db/Homo_sapiens.GRCh38.cds.fa -input_type fasta -out $data/blastdb/makeblastdb_output -parse_seqids

#blast our fasta sequences against the reference db
# db: blast database file name
# query: input file name
# evalue: expectation value E for saving hits
# outfmt: alignment view options. 6: Tabular
# out: output file name
# max_target_seqs: only one hit per contig (one annotation per gene)
/softwares/ncbi-blast-2.10.1+/bin/blastn -db $data/blastdb/makeblastdb_output -query $data/transdecoder_results/Trinity_RF.fasta.transdecoder.cds -evalue 1e-4 -outfmt 6 -out $data/blast_alignment/blast_alignment_output -max_target_seqs 1