import os
import pgrtk
import requests
import json

gene_names = "LILRA2 CDC42EP5 LAIR1 LAIR2 LENG8 LENG9 LILRB4 LILRA1 LILRA3 LILRA4 LILRA5 LILRA6 LILRB1 LILRB2 LILRB3 LILRB5 TTYH1".split()
json_out = open("b12x_features.jsonl", "w") 
for gene_name in gene_names:
    f = open( f"results/genes/{gene_name}_Gene_Structure.txt" )
    for row in f:
        row = row.strip()
        if row[0] == "#":
            continue
        t, r, s = row.split()
        r = int(r)
        d = { "locus": gene_name, "term": t, "rank": r, "sequence": s.upper() }
        headers = {'Content-type': 'application/json'}
        r = requests.post('https://feature.b12x.org:443/features', json=d, headers=headers)
        out = r.text.replace("\n","")
        print(out, file=json_out)

json_out.close()
