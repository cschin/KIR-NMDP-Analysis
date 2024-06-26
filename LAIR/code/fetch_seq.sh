cat << EOF  | tr " " "\t" > regions_interest
LAIR hg38_tagged.fa chr19_hg38 54150000 54720000
LAIR_sub_1 hg38_tagged.fa chr19_hg38 54618706 54672913
LILRB1 hg38_tagged.fa chr19_hg38 54601705 54654041
EOF
mkdir -p results
\time -v pgr-fetch-seqs /wd/pgr-tk-demo-data/data/pgr-tk-HGRP-y1-evaluation-set-v0 \
    -r ./regions_interest > ./results/ROI_seq.fa
\time -v pgr-query /wd/pgr-tk-demo-data/data/pgr-tk-HGRP-y1-evaluation-set-v0 \
    ./results/ROI_seq.fa ./results/pg_seqs --merge-range-tol 100000
python filter.py > LAIR_HPRC.fa
