import os
import pgrtk

gene_names = "LILRA2 CDC42EP5 LAIR1 LAIR2 LENG8 LENG9 LILRB4 LILRA1 LILRA3 LILRA4 LILRA5 LILRA6 LILRB1 LILRB2 LILRB3 LILRB5 TTYH1".split()
for gene_name in gene_names:
    sdb = pgrtk.SeqIndexDB()
    sdb.load_from_fastx(f"results/genes/{gene_name}_DNA_rank1.fa")
    sinfo = sdb.seq_info.copy()
    lengths = []
    for sid in sinfo:
        ctg, _1, seq_len = sinfo[sid]
        lengths.append(seq_len)
    median_len = lengths[int(len(lengths)/2)] 
    target_len = int(median_len * 1.2);
    os.system(f"""
pgr-pbundle-decomp results/genes/{gene_name}_DNA_rank1.fa -w 80 -r 1 --min-span 12 --bundle-length-cutoff 10 --min-branch-size 0 --min-cov 0 results/genes/{gene_name}_DNA_rank1

pgr-pbundle-bed2dist results/genes/{gene_name}_DNA_rank1.bed results/genes/{gene_name}_DNA_rank1 -l

pgr-pbundle-bed2svg results/genes/{gene_name}_DNA_rank1.bed results/genes/{gene_name}_DNA_rank1. --html --ddg-file results/genes/{gene_name}_DNA_rank1.ddg --offsets results/genes/{gene_name}_DNA_rank1.offset --track-panel-width 800 --track-range {target_len} 
    """)

    sdb = pgrtk.SeqIndexDB()
    sdb.load_from_fastx(f"results/genes/{gene_name}_DNA.fa")
    sinfo = sdb.seq_info.copy()
    lengths = []
    for sid in sinfo:
        ctg, _1, seq_len = sinfo[sid]
        lengths.append(seq_len)
    median_len = lengths[int(len(lengths)/2)] 
    target_len = int(median_len * 1.2);
    os.system(f"""
pgr-pbundle-decomp results/genes/{gene_name}_DNA.fa -w 80 -r 1 --min-span 12 --bundle-length-cutoff 10 --min-branch-size 0 --min-cov 0 results/genes/{gene_name}_DNA

pgr-pbundle-bed2dist results/genes/{gene_name}_DNA.bed results/genes/{gene_name}_DNA -l

pgr-pbundle-bed2svg results/genes/{gene_name}_DNA.bed results/genes/{gene_name}_DNA. --html --ddg-file results/genes/{gene_name}_DNA.ddg --offsets results/genes/{gene_name}_DNA.offset --track-panel-width 800 --track-range {target_len} 
    """)
