# Trinity parameters:
# Are sequences paired or unpaired reads? Paired reads. 
# Read1: left, Read2: right
# SS lib type: orientation des RNA-seq reads: RF ; fr first strand pour lexogen (seq de l'ARN correspond au brin reverse: R1 sur le brin reverse), correspond à rf dans Trinity
# seqType: fq
# CPU: 4

data=/ifb/data/mydatalocal/data_tp_ngs


# creating a list containing the names of the forward reads files and a list of the reverse reads files:
paired_forward_files=$(ls $data/trimmomatic_outputs/*forward_paired.fq.gz | paste -d "," -s)
paired_reverse_files=$(ls $data/trimmomatic_outputs/*reverse_paired.fq.gz | paste -d "," -s)

# running Trinity
Trinity --seqType fq --max_memory 14G --SS_lib_type RF --output $data/trinity_outputs/ --left $paired_forward_files --right $paired_reverse_files --CPU 4