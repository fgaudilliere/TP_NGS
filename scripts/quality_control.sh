# Checking the quality of the data

# boucle for balayant chaque fichier et exécutant fastqc pour chaque fichier
# Créer une liste de fichiers avec ls ou utiliser *
# Donne overview de la qualité. Fichiers de sortie: html à download 
# Threads: nombre de coeurs à utiliser (4 coeurs sur la VM)

data_fastq=/ifb/data/mydatalocal/fastq_sequences/Lib*.fastq.gz

for sequence_file in data_fastq
do

fastqc -t 4 $sequence_file

done