import pgrtk

sdb = pgrtk.SeqIndexDB()
sdb.load_from_fastx("./results/pg_seqs.000.fa", 80, 33, 1, 8)
sinfo = sdb.seq_info.copy()

for sid in sinfo:
    ctg, _1, seq_len = sinfo[sid]
    if seq_len < 20000:
        continue
    seq = sdb.get_seq_by_id(sid)
    print(f">{ctg}")
    print(pgrtk.u8_to_string(seq))
