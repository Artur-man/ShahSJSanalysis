# Processing RNA and ADT reads
cellranger count \
--id=bri793_cite \
--libraries=library.csv \
--transcriptome=/home/wk93w/yard/apps/cellranger/refdata-cellranger-GRCh38-3.0.0 \
--feature-ref=feature_ref_cite.csv

# Processing RNA and Hashtag reads
cellranger count \
--id=bri793 \
--libraries=library.csv \
--transcriptome=//home/wk93w/yard/apps/cellranger/refdata-cellranger-GRCh38-3.0.0 \
--feature-ref=feature_ref.csv

# Processing RNA and Hashtag reads
cellranger vdj \
--id=bri795 \
--fastqs=//home/wk93w/stevens_johnson/HKGL3DRXX/fastqs/BRI-795 \
--reference=/home/wk93w/yard/apps/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
--sample=BRI-795
