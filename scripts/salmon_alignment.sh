data=/ifb/data/mydatalocal/data_tp_ngs

# Creating a salmon index

#salmon index -t $data/trinity_results/Trinity_RF.fasta -i $data/salmon_index -p 4

# Alignment on created index

#for n in 1 2 3 4 5 6
#do
#salmon quant -i $data/salmon_index -l A -1  $data/trimmomatic_outputs/Lib${n}_output_forward_paired.fq.gz -2 $data/trimmomatic_outputs/Lib${n}_output_reverse_paired.fq.gz --validateMappings --gcBias -p 4 -o $data/salmon_alignment
#done

# Bons résultats : > 80% de reads alignés
# Ici: inserts trop petits, reads se chevauchent
# Réalignement avec salmon: ne fonctionne pas -> seulement 40% de reads alignés
# Correction : aligner avec des données single-end

# Alignement en single-end:

for n in 1 2 3 4 5 6
do
salmon quant -i $data/salmon_index -l A -r $data/trimmomatic_outputs/Lib${n}_output_forward_paired.fq.gz --validateMappings --gcBias -p 4 -o $data/salmon_alignment_single_end${n}_forward
salmon quant -i $data/salmon_index -l A -r $data/trimmomatic_outputs/Lib${n}_output_reverse_paired.fq.gz --validateMappings --gcBias -p 4 -o $data/salmon_alignment_single_end${n}_reverse
done
