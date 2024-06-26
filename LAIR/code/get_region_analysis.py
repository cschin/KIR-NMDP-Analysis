import os
import pgrtk

gene_names = "LILRA2 CDC42EP5 LAIR1 LAIR2 LENG8 LENG9 LILRB4 LILRA1 LILRA3 LILRA4 LILRA5 LILRA6 LILRB1 LILRB2 LILRB3 LILRB5 TTYH1".split()

annotation = []
with open("results/all_annotation.bed") as f:
    for r in f:
        r = r.strip().split("\t")
        annotation.append(r)

for gene_name in gene_names:
    
    f = open(f"results/annotation_{gene_name}.bed", "w")
    for r in annotation:
        if r[3] == gene_name:
            print("\t".join(r), file=f)
    f.close() 

    os.system(f"""
        pgr-pbundle-bed2svg LAIR_HPRC.bed results/genes/{gene_name}_region. --html --ddg-file LAIR_HPRC.ddg --offsets LAIR_HPRC.offset --track-panel-width 800 --track-range 650000 --annotation-region-bedfile results/annotation_{gene_name}.bed --annotation-region-stroke-width 8
    """)

