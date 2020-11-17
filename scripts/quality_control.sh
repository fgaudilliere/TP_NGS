# Checking the quality of the data

# boucle for balayant chaque fichier et exécutant fastqc pour chaque fichier
# Créer une liste de fichiers avec ls ou utiliser *
# Donne overview de la qualité. Fichiers de sortie: html à download 
# Threads: nombre de coeurs à utiliser (4 coeurs sur la VM)

data=/ifb/data/mydatalocal/data_tp_ngs/
# data variable: path to file containing all data needed for the practical

data_fastq=$data/fastq_sequences/Lib*.fastq.gz
# data_fastq: list of the fastq files names


# Following lines: loop to assess the quality of the data using fastqc
for sequence_file in data_fastq
do
fastqc -t 4 $sequence_file
done


# Following loop: using Trimmomatic to clean the data
# Options used here: 
# threads: 8 cores used
# phreds33: one of two options, used because we were told to
# IlluminaClip: clipping sequences corresponding to Illumina adaptators (the reference sequences are given in the adapt.fasta file)
# Headcrop: cropping the beginning of reads by a given number of bases
# Minlen: discarding reads below a given size

# Structure of the Trimmomatic command: 
# Call program: java -jar <path to Trimmomatic> PE
# Give number of threads and phred option: -threads 8 -phred33
# Give input files: two files corresponding to paired reads (R1 and R2)
# Give desired output files names for: paired forward output, unpaired forward output (reads for which one of the paired reads was discarded), paired reverse output, unpaired reversed output
# Optional features: IlluminaClip with options and path to reference file for Illumina adaptator sequences, MINLEN with length below which to discard, HEADCROP with number of nt to crop at the beginning of reads
for number in 1 2 3 4 5 6
do

java -jar /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 -phred33 $data/fastq_sequences/Lib${number}_31_20_S${number}_R1_001.fastq.gz $data/fastq_sequences/Lib${number}_31_20_S${number}_R2_001.fastq.gz $data/trimmomatic_outputs/Lib${number}_output_forward_paired.fq.gz $data/trimmomatic_outputs/Lib${number}_output_forward_unpaired.fq.gz $data/trimmomatic_outputs/Lib${number}_output_reverse_paired.fq.gz $data/trimmomatic_outputs/Lib${number}_output_reverse_unpaired.fq.gz ILLUMINACLIP:$data/fastq_sequences/adapt.fasta:2:30:10 MINLEN:100 HEADCROP:9

done