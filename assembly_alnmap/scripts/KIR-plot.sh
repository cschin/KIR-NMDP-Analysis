#pgr-pbundle-bed2dist KIR_HPRC.bed KIR_HPRC.
#bash ../Immuannot/scripts/immuannot.sh -c assemble_results.fa -r Data-2023Oct27 -o immuannot_KIR-HPRC.out -t 20
# python get_gene_seq_from_immuannot.py assemble_results.fa immuannot_KIR-HPRC.out.gtf.gz > KIR-annotation.bed

pgr-pbundle-decomp assemble_results.fa -w 128 -r 12 --min-span 32 --bundle-length-cutoff 120 --min-branch-size 4 -i assemble_results.include assemble_results --min-cov 0 

pgr-pbundle-bed2dist assemble_results.bed assemble_results. -l

pgr-pbundle-bed2svg assemble_results.bed assemble_results. --html --ddg-file assemble_results.ddg --offsets assemble_results.offset --track-panel-width 800 --track-range 350000 --annotation-region-bedfile KIR-annotation.bed --annotation-region-stroke-width 5

#pgr-pbundle-bed2svg assemble_results_filtered.bed assemble_results_f. --html --ddg-file assemble_results.ddg --offsets assemble_results.offset --track-panel-width 800 --track-range 350000 --annotation-region-bedfile KIR-annotation.bed --annotation-region-stroke-width 5

#pgr-alnmap --preset overwrite -w 7 -k 13 -r 1 -m 1   assemble_results.fa kirnmdp.200.fasta prob_aln_test
