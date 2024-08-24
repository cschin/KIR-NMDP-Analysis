import pgrtk
import os
import re

def partition_by_case(s):
    return re.findall(r'[A-Z]+|[a-z]+', s)

sdb = pgrtk.SeqIndexDB()
sdb.load_from_fastx("./LAIR_HPRC.fa")
sinfo = sdb.seq_info.copy()

out_bed_file = open("results/all_annotation.bed", "w")
out_bed_file_rank1 = open("results/all_annotation_rank1.bed", "w")
out_seq_files = {}
out_seq_files_rank1 = {}
out_cds_files = {}
out_cds_files_rank1 = {}
out_prot_files = {}
out_prot_files_rank1 = {}
out_gene_structure_files = {}
out_gene_structure_files_rank1 = {}

haplotype_to_seq = {}
for sid in sinfo:
    ctg, _1, seq_len = sinfo[sid]
    key = "_".join(ctg.split(".")[:2])
    seq = sdb.get_seq_by_id(sid)
    haplotype_to_seq.setdefault(key, [])
    haplotype_to_seq[key].append( (ctg, seq) )


for key in haplotype_to_seq:
    with open(f"results/target_seq_{key}.fa","w") as f:
        for ctg, seq in haplotype_to_seq[key]:
            print(f">{ctg}", file=f)
            print(pgrtk.u8_to_string(seq), file=f)
    os.system(f"""./miniprot results/target_seq_{key}.fa LAIR_protein_ref.fa -j 2 --trans --aln --max-intron-out 20000 -G 20000 --outs=0.975 --outc=0.8 --gff > results/LAIR_HPRC_annotation_{key}.gff""")
    with open(f"results/LAIR_HPRC_annotation_{key}.gff") as f:
        for r in f:
            r = r.strip().split("\t")
            if r[0] == "##ATN":
                ATN = r[1]
                seq = ATN
                seq = seq.replace("-", "")

                exon_introns = partition_by_case(seq)
                exons = enumerate([s for s in exon_introns if s[0].isupper()])
                introns = enumerate([s for s in exon_introns if s[0].islower()])

                seq = ATN.upper()
                continue

            if r[0] == "##PAF":
                ctg = r[6]
                continue

            if r[0] == "##STA":
                STA = r[1]
                continue

            if len(r) < 3:
                continue

            if r[2] == "mRNA":
                ATN = ATN.replace("a","").replace("t","").replace("g","").replace("c","").replace("-","")
                bgn, end = int(r[3]), int(r[4])
                strand = r[6]
                """ID=MP000007;Rank=2;Identity=0.6973;Positive=0.8016;Frameshift=2;StopCodon=1;Target=LILRA4"""
                annotation_dict = dict([_.split("=") for _ in r[8].split(" ")[0].split(";")])
                target = annotation_dict["Target"]
                rank = int(annotation_dict["Rank"])
                if rank == 1 and "Frameshift" not in annotation_dict:
                    print(r[0], r[3], r[4], target, "#000", r[8], file=out_bed_file_rank1, sep="\t") 
                    if target not in out_seq_files_rank1:
                        out_seq_files_rank1.setdefault(target, open(f"results/genes/{target}_DNA_rank1.fa", "w"))
                        out_cds_files_rank1.setdefault(target, open(f"results/genes/{target}_CDS_rank1.fa", "w"))
                        out_prot_files_rank1.setdefault(target, open(f"results/genes/{target}_PROT_rank1.fa", "w"))
                        out_gene_structure_files_rank1.setdefault(target, open(f"results/genes/{target}_Gene_Structure_rank1.txt", "w"))
                    out_seq_file_rank1 = out_seq_files_rank1[target]
                    DNA_seq = seq
                    print(">{}:{}-{} {}".format(ctg, r[3], r[4], r[8]), file=out_seq_file_rank1)
                    print(DNA_seq, file=out_seq_file_rank1)

                    out_cds_file_rank1 = out_cds_files_rank1[target]
                    print(">{}:{}-{} {}".format(ctg, r[3], r[4], r[8]), file=out_cds_file_rank1)
                    print(ATN, file=out_cds_file_rank1)

                    out_prot_file_rank1 = out_prot_files_rank1[target]
                    print(">{}:{}-{} {}".format(ctg, r[3], r[4], r[8]), file=out_prot_file_rank1)
                    print(STA, file=out_prot_file_rank1)

                    out_gene_structure_file_rank1 = out_gene_structure_files_rank1[target]
                    print("# {}:{}-{} {}".format(ctg, r[3], r[4], r[8]), file=out_gene_structure_file_rank1)
                    for rank, s in exons:
                        print("exon {} {}".format(rank, s), file=out_gene_structure_file_rank1)
                    for rank, s in introns:
                        print("introns {} {}".format(rank, s), file=out_gene_structure_file_rank1)



                print(r[0], r[3], r[4], target, "#000", r[8], file=out_bed_file, sep="\t") 
                if target not in out_seq_files:
                    out_seq_files.setdefault(target, open(f"results/genes/{target}_DNA.fa", "w"))
                    out_cds_files.setdefault(target, open(f"results/genes/{target}_CDS.fa", "w"))
                    out_prot_files.setdefault(target, open(f"results/genes/{target}_PROT.fa", "w"))
                    out_gene_structure_files.setdefault(target, open(f"results/genes/{target}_Gene_Structure.txt", "w"))
                out_seq_file = out_seq_files[target]
                DNA_seq = seq
                print(">{}:{}-{} {}".format(ctg, r[3], r[4], r[8]), file=out_seq_file)
                print(DNA_seq, file=out_seq_file)

                out_cds_file = out_cds_files[target]
                print(">{}:{}-{} {}".format(ctg, r[3], r[4], r[8]), file=out_cds_file)
                print(ATN, file=out_cds_file)

                out_prot_file = out_prot_files[target]
                print(">{}:{}-{} {}".format(ctg, r[3], r[4], r[8]), file=out_prot_file)
                print(STA, file=out_prot_file)

                out_gene_structure_file = out_gene_structure_files[target]
                print("# {}:{}-{} {}".format(ctg, r[3], r[4], r[8]), file=out_gene_structure_file)
                for rank, s in exons:
                    print("exon {} {}".format(rank, s), file=out_gene_structure_file)
                for rank, s in introns:
                    print("introns {} {}".format(rank, s), file=out_gene_structure_file)

out_bed_file.close()
out_bed_file_rank1.close()
[f.close() for f in out_seq_files.values()]
[f.close() for f in out_seq_files_rank1.values()]
[f.close() for f in out_cds_files.values()]
[f.close() for f in out_cds_files_rank1.values()]
