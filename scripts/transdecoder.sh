data=/ifb/data/mydatalocal/data_tp_ngs
cd $data/transdecoder_results

# transdecoder: finding ORFs

TransDecoder.LongOrfs -t $data/trinity_results/Trinity_RF.fasta --gene_trans_map $data/trinity_results/Trinity_RF.fasta.gene_trans_map --output_dir $data/transdecoder_results/ -m 100 -S


# transdecoder: predicting coding regions

TransDecoder.Predict -t $data/trinity_results/Trinity_RF.fasta -O $data/transdecoder_results/ --single_best_only