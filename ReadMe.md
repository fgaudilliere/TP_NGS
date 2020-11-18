---
output:
  html_document: default
  pdf_document: default
---
# NGS Practicals

This is the description of the code used in my NGS practicals in November 2020.

## Data download (download_data.sh)

To download the data used in these practicals, we use the following command:
```
wget -r --ftp-user=igfl-UE_NGS_2020 --ftp-password=UE_NGS_2020 ftp://sharegate-igfl.ens-lyon.fr/Projet_31_20_UE_NGS_2020/
```
The line specifies the adress of the directory containing what we want, as well as the username and password we need to get there.

The files are under the fastq.gz format, that is a compressed fastq format. fastq files contain read sequences as well as information on the quality of those reads.

After downloading it, I moved the data to a new directory: it is now in a specific directory called data_tp_ngs that will contain all the data for these practicals (one directory per type of files).

## Quality control (quality_control.sh)

### Running fastqc

To start, we run fastqc on our fastq.gz files: this program analyzes the contents of fastq files and returns a report on each file containing graphs summing up things like per base sequence quality, per tile sequence quality, per sequence quality scores, per base sequence content etc. 

#### Structure of the script

We build a for loop that goes through each fastq.gz file and runs fastqc on that file. 

#### fastqc parameters

The command to call fastqc is fairly simple and goes as follows:
```
fastqc -t 4 sequence_file_name.fastq.gz -o output_directory
```


The -t corresponds to threads, i.e. the number of cores to be dedicated to the task. The -o corresponds to the output directory to store the results in. 

#### Results

The output files are stored in a file called fastqc_results.

We can see in the outputs that we have per base sequence content anomalies at the beginning of the reads (on approximately 8-9 pairs) and that the end of reads often corresponds to the sequencing of Illumina adaptators. Therefore, the data is not clean, and we are going to need to trim it. 

### Running Trimmomatic

To trim the reads, we use Trimmomatic (NB: here we run it separately, but there is also an option in Trinity to run Trimmomatic).

#### Structure of the script

Like before, we build a for loop and for each iteration, we give Trimmomatic two paired files (forward and reverse reads corresponding to a given sample).

#### Trimmomatic parameters

- threads: number of cores to use.
- phreds33: one of two options, either phreds33 or phreds64.
- IlluminaClip: clips off sequences corresponding to Illumina adaptators (the reference sequences are given in a file called adapt.fasta).
- Headcrop: crops the beginning of reads by a given number of bases. Here we crop the 9 first bases (since they display per base sequence content anomalies).
- Minlen: discards reads below a given size. Here we discard reads below 100 bases (standard practice).

The command is structured as follows:
```
java -jar /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 -phred33 name_of_R1_sequence.fastq.gz name_of_R2.fastq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:adapt.fasta:2:30:10 MINLEN:100 HEADCROP:9
```
The unpaired outputs correspond to reads for which the corresponding paired reads were discarded in the process. 

#### Trimmomatic outputs

The outputs are stored in a specific directory called trimmomatic_ouputs, which can be found in the data directory. 


### Running fastqc on Trimmomatic results

Now that our reads have been trimmed, we run fastqc on Trimmomatic's output to make sure that everything went as planned and that our new data is cleaner. We repeat the operation described in the fastqc section.

## Data assembly (data_assembly.sh)

After running the quality control script, we have clean data files which we can use for data assembly. 

This is done using a program called Trinity.

### Structure of the script

We start by building two strings containing the names of 1) all the forward paired outputs from Trimmomatic, and 2) all the reverse paired outputs from Trimmomatic. 

We then feed these files to Trinity for data assembly. 

NB: this is what I would have done theoretically, but given that Trinity actually takes two days to run on the data used here, the teachers ran it ahead of the practicals and just sent us the results.

### A few tips on running demanding scripts

You can keep a log of everything a program prints while it runs with a command of the following type: 
```
./script.sh > log_file.txt
```

You can run a program in the background (so as to still be able to use the terminal while it's running) with the following command:
```
& ./script.sh &
```
You can see which processes are running in the background by typing:
```
ps
```
For more details (such as how many cores are busy), you can use: 
```
htop
```
If you want to kill a specific process using its PID number, do:
```
kill <PID number of the process>
```

The nohup command can be used to transfer messages relative to the execution of a program, which by default are displayed in the terminal, into a file created for the occasion. 
For example, typing: 
```
nohup./script.sh
```
will create a nohup.out file which is a log of the script's execution. Be careful though, if you run several scripts with nohup without specifying a specific file to write the log in, all logs will be written in the same nohup.out file one after the other. To specify a file for a given script, you can do: 
```
./script.sh > nohup.log_file.txt
```


The previous commands can be combined with:
```
./script.sh > & nohup.log_file.txt &
```

### Trinity parameters

- seqType: specifies the sequence format. Here we use fq (fastq) files. 
- left and right (for paired reads): enter reads R1 in left and reads R2 in right.
- SS_lib_type: orientation of the RNA-seq reads. Here it is RF. 
- CPU: number of CPU to use. Here it is 4.
- max_memory: maximum RAM to be dedicated to the process. Usually you should put in your maximum RAM minus 2 (to keep some calculation power in case you need to do something else simultaneously). 
- output: directory in which Trinity should store the outputs.

The command is structured as follows: 
```
Trinity --seqType fq --max_memory 14G --SS_lib_type RF --output $data/trinity_outputs/ --left list_of_paired_forward_files --right list_of_paired_reverse_files --CPU 4
```

### Output files

The output files are stored in a directory called trinity_outputs. We have two files: a .fasta and a .fasta.gene_trans_map. 

In the .fasta file, we have lines (starting with a '>') of the form TRINITY_DNX_cX_gX_iX, with TRINITY_DNX_cX identifying a cluster, TRINITY_DNX_cX_gX identifying a gene in the cluster and TRINITY_DNX_cX_gX_iX identifying an isoform of that gene. These lines are followed by another line with the sequence corresponding to the isoform labeled by the identifier line.  
The identifier lines can be used to count the number of genes, for instance with a line like: 
```
grep ">" Trinity_RF.fasta |cut -f1,2,3,4 -d "_" |sort |uniq |less
```
Basically, we select lines starting with '>' in the Trinity output file, cut them at each "_" symbol and keep the 1st, 2nd, 3rd and 4th parts, which we then sort; finally we only keep one occurrence of a given character string and we show the final result (using less). This operation allows us to keep the TRINITY_DNX_cX_gX part of the identifier lines (corresponding to a gene) and count how many different genes we have. 

The .fasta.gene_trans_map file contains a table with correspondences between genes and isoforms. 


## Reads alignment (salmon_alignment.sh)

To annotate the data, we use a program called salmon. Salmon is divided in several sub-programs, and we use two of them: salmon index to build an index on which to map the reads, and salmon quant to align and quantify the reads. 

### Salmon index

The command takes the following parameters:

- t: fasta file containing the assembled reads (Trinity output)
- i: directory in which to create the output
- p: number of cores to use in the calculation

The command line looks like this:
```
salmon index -t $data/trinity_results/Trinity_RF.fasta -i $data/salmon_index -p 4
```

The output is stored in a file called salmon_index. 

### Salmon alignment

We then use the index we just created to align our reads (the paired reads corresponding to the outputs of Trimmomatic) with salmon quant. 

#### Structure of the script

Just like in previous steps (fastqc or Trimmomatic), we need to process several files one after the other, so we create a for loop going through each pair of files (forward reads and reverse reads) and applying salmon quant to that pair of files.

#### salmon quant parameters

The parameters for salmon quant are:

-i: index to use (the salmon index we just built)
- l: library format. We put in A (for Automatic) so that salmon will look at the files we enter and determine it itself. 
- 1: forward reads.
- 2: reverse reads.
- validateMapping: option that is default in the newest versions of salmon, but here we need to indicate it. 
- o: directory in which to put output files.

The command line looks like this: 
```
salmon quant -i $data/salmon_index -l A -1  forward_paired_output.fq.gz -2 reverse_paired_output.fq.gz --validateMappings -o $data/salmon_alignment
```

#### salmon quant output

Normally, when aligning reads, a result is considered good when > 80% of reads are aligned. However, here, we only reach 40% of aligned reads. This is likely due to a problem during sequencing: inserts (RNA fragments) were too small, and therefore the forward and reverse reads overlap. 

#### Correcting the salmon quant output

To correct this, we can run salmon quant a second time, but this time we do as if the reads were not paired but single. The parameters are almost all the same except for -1 and -2 which are replaced by -r, and we only give salmon quant one file at a time instead of two. 

The command line looks like this:
```
salmon quant -i $data/salmon_index -l A -r forward_paired_output.fq.gz --validateMappings -o $data/salmon_alignment_single_end
```
