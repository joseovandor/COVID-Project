# Input data

Not tracked by git. The scripts download each file here when it is missing:

```
data/raw/GSE172274/GSE172274_raw_counts_GRCh38.p13_NCBI.tsv.gz
data/raw/GSE179277/GSE179277_all_adult_and_ped_samples_combined_counts_unfiltered.csv.gz
data/raw/GSE231409/GSE231409_BRAVE_RNASeq_counts.csv.gz
data/sc_rna_seq/Yoshida_airway.h5ad
```

Sample metadata are retrieved from GEO with GEOquery. The H5AD file is the
Yoshida et al. (2022) airway dataset distributed through CELLxGENE
(https://datasets.cellxgene.cziscience.com/f9efb73e-f116-46b5-a775-d4233e758024.h5ad).
If a download fails, place the file manually at the path above.
