import os
import pgrtk

gene_names = ("LILRB3", "LILRA6", "LILRB5", 
              "LILRB2", "LILRA3", "LILRA5", 
              "LILRA4", "LAIR1", "TTYH1", 
              "LENG8", "LENG9", "CDC42EP5", 
              "LAIR2", "LILRA2", "LILRA1", 
              "LILRB1", "LILRB4")

f_count = open("unique_count.txt", "w")
for gene_name in gene_names:
    print(gene_name)

    sdb = pgrtk.SeqIndexDB()
    print(f"results/genes/{gene_name}_DNA.fa")
    sdb.load_from_fastx(f"results/genes/{gene_name}_DNA.fa")
    sinfo = sdb.seq_info.copy()
    seqs = []
    
    seqs = []    
    seq_set = set() 
    for sid in sinfo:
        ctg, _, slen = sinfo[sid]
        seq = sdb.get_seq_by_id(sid)
        if tuple(seq) not in seq_set and "_alt_" not in ctg:
            seq_set.add(tuple(seq))
            seqs.append((ctg, seq))

    f = open(f"results/genes/{gene_name}_DNA_unique.fa", "w")
    lengths = []
    for ctg, seq in seqs:
        print(f">{ctg}", file=f)
        print(pgrtk.u8_to_string(seq), file=f)
        lengths.append(len(seq))
    f.close()
    print(gene_name, "DNA", len(sinfo), len(seqs), file=f_count)
    median_len = lengths[int(len(lengths)/2)] 
    target_len = int(median_len * 1.2);
    w = 80
    m = 12
    os.system(f"""
pgr-pbundle-decomp results/genes/{gene_name}_DNA_unique.fa -w {w} -r 1 --min-span {m} --bundle-length-cutoff 10 --min-branch-size 0 --min-cov 0 results/genes/{gene_name}_DNA_unique

pgr-pbundle-bed2dist results/genes/{gene_name}_DNA_unique.bed results/genes/{gene_name}_DNA_unique -l

pgr-pbundle-bed2svg results/genes/{gene_name}_DNA_unique.bed results/genes/{gene_name}_DNA_unique. --html --ddg-file results/genes/{gene_name}_DNA_unique.ddg --offsets results/genes/{gene_name}_DNA_unique.offset --track-panel-width 800 --track-range {target_len} 
    """)



    sdb = pgrtk.SeqIndexDB()
    print(f"results/genes/{gene_name}_CDS.fa")
    sdb.load_from_fastx(f"results/genes/{gene_name}_CDS.fa")
    sinfo = sdb.seq_info.copy()
    seqs = []
    
    seqs = []    
    seq_set = set() 
    for sid in sinfo:
        ctg, _, slen = sinfo[sid]
        
        seq = sdb.get_seq_by_id(sid)
        if tuple(seq) not in seq_set and "_alt_" not in ctg:
            seq_set.add(tuple(seq))
            seqs.append((ctg, seq))

    f = open(f"results/genes/{gene_name}_CDS_unique.fa", "w")
    lengths = []
    for ctg, seq in seqs:
        print(f">{ctg}", file=f)
        print(pgrtk.u8_to_string(seq), file=f)
        lengths.append(len(seq))
    f.close()
    print(gene_name, "CDS", len(sinfo), len(seqs), file=f_count)
    median_len = lengths[int(len(lengths)/2)] 
    target_len = int(median_len * 1.2);
    m = 8
    w = 32
    os.system(f"""
pgr-pbundle-decomp results/genes/{gene_name}_CDS_unique.fa -w {w} -r 1 --min-span {m} --bundle-length-cutoff 0 --min-branch-size 0 --min-cov 0 results/genes/{gene_name}_CDS_unique

pgr-pbundle-bed2dist results/genes/{gene_name}_CDS_unique.bed results/genes/{gene_name}_CDS_unique -l

pgr-pbundle-bed2svg results/genes/{gene_name}_CDS_unique.bed results/genes/{gene_name}_CDS_unique. --html --ddg-file results/genes/{gene_name}_CDS_unique.ddg --offsets results/genes/{gene_name}_CDS_unique.offset --track-panel-width 800 --track-range {target_len} 
    """)




    sdb = pgrtk.SeqIndexDB()
    print(f"results/genes/{gene_name}_PROT.fa")
    sdb.load_from_fastx(f"results/genes/{gene_name}_PROT.fa")
    sinfo = sdb.seq_info.copy()
    seqs = []
    
    seqs = []    
    seq_set = set() 
    for sid in sinfo:
        ctg, _, slen = sinfo[sid]
        seq = sdb.get_seq_by_id(sid)
        if tuple(seq) not in seq_set and "_alt_" not in ctg:
            seq_set.add(tuple(seq))
            seqs.append((ctg, seq))

    f = open(f"results/genes/{gene_name}_PROT_unique.fa", "w")
    for ctg, seq in seqs:
        print(f">{ctg}", file=f)
        print(pgrtk.u8_to_string(seq), file=f)
    f.close()
    print(gene_name, "PROT", len(sinfo), len(seqs), file=f_count)
    os.system(f"""
muscle -in results/genes/{gene_name}_PROT_unique.fa -html -out results/genes/{gene_name}_PROT_unique.html
    """)

f_count.close()


