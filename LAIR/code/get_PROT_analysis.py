import os
import pgrtk

gene_names = "LILRA2 CDC42EP5 LAIR1 LAIR2 LENG8 LENG9 LILRB4 LILRA1 LILRA3 LILRA4 LILRA5 LILRA6 LILRB1 LILRB2 LILRB3 LILRB5 TTYH1".split()
for gene_name in gene_names:
    os.system(f"""
muscle -in results/genes/{gene_name}_PROT_rank1.fa -html -out results/genes/{gene_name}_PROT_rank1.html
    """)
    os.system(f"""
muscle -in results/genes/{gene_name}_PROT.fa -html -out results/genes/{gene_name}_PROT.html
    """)

