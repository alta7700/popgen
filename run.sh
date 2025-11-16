ALL_GENES_FILE=corona_genes_all.fasta
GENES_FILE_PREFIX=corona
GENE1=S
GENE2=N
GENE1_FN=${GENES_FILE_PREFIX}_${GENE1}
GENE2_FN=${GENES_FILE_PREFIX}_${GENE2}

GBLOCKS_PATH="$(dirname "$(which python)")/gblocks"

python fetch_genes.py --output $ALL_GENES_FILE --genes ${GENE1} ${GENE2}
python split_genes.py --input $ALL_GENES_FILE --genes ${GENE1} ${GENE2} --prefix ${GENES_FILE_PREFIX}
muscle5 -align ${GENE1_FN}.fasta -output ${GENE1_FN}_aligned.fasta
muscle5 -align ${GENE2_FN}.fasta -output ${GENE2_FN}_aligned.fasta
python gblocks_wrap.py ${GENE1_FN}_aligned.fasta --gblocks $GBLOCKS_PATH
python gblocks_wrap.py ${GENE2_FN}_aligned.fasta --gblocks $GBLOCKS_PATH
muscle5 -align ${GENE1_FN}_aligned.gblocks.fasta -output ${GENE1_FN}_aligned2.fasta
muscle5 -align ${GENE2_FN}_aligned.gblocks.fasta -output ${GENE2_FN}_aligned2.fasta

