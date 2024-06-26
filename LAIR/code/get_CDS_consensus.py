import os
import pgrtk

gene_names = "LILRA2 CDC42EP5 LAIR1 LAIR2 LENG8 LENG9 LILRB4 LILRA1 LILRA3 LILRA4 LILRA5 LILRA6 LILRB1 LILRB2 LILRB3 LILRB5 TTYH1".split()


f = open("results/CDS_consensus.fa", "w")

for gene_name in gene_names:
    print(gene_name)

    sdb = pgrtk.SeqIndexDB()
    print(f"results/genes/{gene_name}_CDS.fa")
    sdb.load_from_fastx(f"results/genes/{gene_name}_CDS.fa")
    sinfo = sdb.seq_info.copy()
    seqs = []
    
    seq_lengths = []
    for sid in sinfo:
        _, _, slen = sinfo[sid]
        seq_lengths.append(slen)
        
    median = seq_lengths[int(len(seq_lengths)/2)]
    for sid in sinfo:
        seq = sdb.get_seq_by_id(sid)
        if abs(len(seq) - median) > median *0.05:
            continue
        seqs.append(seq)

    kmer_size = 23

    consensus=pgrtk.naive_dbg_consensus(seqs, kmer_size, 15)
    print(f">{gene_name}_CDS", file=f)
    print(pgrtk.u8_to_string(consensus), file=f)

    sdb = pgrtk.SeqIndexDB()
    sdb.load_from_fastx(f"results/genes/{gene_name}_CDS_rank1.fa")
    sinfo = sdb.seq_info.copy()
    seqs = []
    
    seq_lengths = []
    for sid in sinfo:
        _, _, slen = sinfo[sid]
        seq_lengths.append(slen)
        
    median = seq_lengths[int(len(seq_lengths)/2)]
    for sid in sinfo:
        seq = sdb.get_seq_by_id(sid)
        if abs(len(seq) - median) > median *0.05:
            continue
        seqs.append(seq)

    consensus=pgrtk.naive_dbg_consensus(seqs, kmer_size, 15)
    print(f">{gene_name}_CDS_rank1", file=f)
    print(pgrtk.u8_to_string(consensus), file=f)

f.close() 


