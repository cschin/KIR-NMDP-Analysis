#bash ../Immuannot/scripts/immuannot.sh -c assemble_results.fa -r Data-2023Oct27 -o immuannot_KIR-HPRC.out -t 20
# python get_gene_seq_from_immuannot.py assemble_results.fa immuannot_KIR-HPRC.out.gtf.gz > KIR-annotation.bed
pgr-pbundle-decomp HLA-ClassII_seq.fa.gz -w 80 -r 6 --min-span 8 --bundle-length-cutoff 250 --min-branch-size 8 HLA-ClassII_seq. --min-cov 0

pgr-pbundle-bed2dist HLA-ClassII_seq.bed HLA-ClassII_seq. -l

pgr-pbundle-bed2svg HLA-ClassII_seq.bed HLA-ClassII_seq. --html --ddg-file HLA-ClassII_seq.ddg --offsets HLA-ClassII_seq.offset --track-panel-width 800  --annotation-region-bedfile immuannot_MHC-c2_annotation_track.bed --annotation-region-stroke-width 5

