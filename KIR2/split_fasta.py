"""
Usage:
  split_fasta.py <sequences.fa.gz>

Options:
  -h --help     Show this screen.

"""
import pgrtk
import gzip
from docopt import docopt


"""
This simple script process the results generated by Immuanno, to get features to feed to GFE server.

Here is the script to generate the GTF fil woth immunannot

bash ../scripts/immuannot.sh -c HLA-ClassII_seq.fa.gz  -r refData-2023Jun05/ -o immuannot_MHC-c2.out -t 20
bash ../scripts/immuannot.sh -c KIR_HPRC.000.fa  -r refData-2023Jun05/ -o immuannot_KIR-HPRC.out -t 20
bash ../scripts/immuannot.sh -c assemble_results.fa.gz  -r refData-2023Jun05/ -o immuannot_KIR_asm.out -t 20

"""


if __name__ == "__main__":
    arguments = docopt(__doc__, version='get_feature_seq.py 0.0')

    sdb = pgrtk.SeqIndexDB()
    sdb.load_from_fastx(arguments['<sequences.fa.gz>'])
    sinfo = sdb.seq_info.copy()
    print(len(sinfo))
    for id_start in range(0, len(sinfo), 12):
        file = open(f"assemble_results-{id_start:03}.fa", "w")
        for id_ in range(id_start, id_start+12):
            if id_ in sinfo:
                print(sinfo[id_])
                seq = sdb.get_seq_by_id(id_)
                seq_name = sinfo[id_][0]
                seq = pgrtk.u8_to_string(seq)
                print(f">{seq_name}", file = file)
                print(f"{seq}", file = file)
        file.close()
