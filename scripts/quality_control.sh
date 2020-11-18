data=/ifb/data/mydatalocal/data_tp_ngs/
# data variable: path to file containing all data needed for the practical

data_fastq=$data/fastq_sequences/Lib*.fastq.gz
# data_fastq: list of the fastq files names

mkdir $data/fastqc_results
# Creating a directory to put the fastqc results


# Following lines: loop to assess the quality of the data using fastqc
for sequence_file in data_fastq
do
fastqc -t 4 $sequence_file -o $data$fastqc_results
done


# Following loop: using Trimmomatic to clean the data
# Options used here: 
# threads: 8 cores used
# phreds33: one of two options, used because we were told to
# IlluminaClip: clipping sequences corresponding to Illumina adaptators 
# (the reference sequences are given in the adapt.fasta file)
# Headcrop: cropping the beginning of reads by a given number of bases
# Minlen: discarding reads below a given size


mkdir $data/trimmomatic_outputs
# Creating a directory to store the trimmomatic outputs

for number in 1 2 3 4 5 6
do
java -jar /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 \
-phred33 $data/fastq_sequences/Lib${number}_31_20_S${number}_R1_001.fastq.gz \
$data/fastq_sequences/Lib${number}_31_20_S${number}_R2_001.fastq.gz \
$data/trimmomatic_outputs/Lib${number}_output_forward_paired.fq.gz \
$data/trimmomatic_outputs/Lib${number}_output_forward_unpaired.fq.gz \
$data/trimmomatic_outputs/Lib${number}_output_reverse_paired.fq.gz \
$data/trimmomatic_outputs/Lib${number}_output_reverse_unpaired.fq.gz \
ILLUMINACLIP:$data/fastq_sequences/adapt.fasta:2:30:10 MINLEN:100 HEADCROP:9
done


# Now that we have trimmed the reads, we look at their quality again
data_trimmomatic=$data/trimmomatic_outputs/*.gz
mkdir $data/fastqc_trimmomatic
# Creating a directory to store the fastqc results

for sequence_file in $data_trimmomatic
do
fastqc -t 4 $sequence_file -o $data/fastqc_trimmomatic

done

