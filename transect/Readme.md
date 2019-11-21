This folder includes data and scripts used for generating Figures 1-3 of the manuscript (and respective Supp. Figures).

There is one folder for each dataset (bacteria, fungi and oomycetes).

Those folders inlude the following files:

**bray_curtis_otu_table_norm.txt** : Table with sample to sample distances (Bray curtis). Note that tables have been normalized (CSS) before distance calculation.

**otu_table_filter.txt** :  Raw OTU-tables.

**rep_seqs.fasta** :  Representative sequences for all OTUs.

**taxonomy_otus.txt** :  Taxonomic assignment for all OTUs.

**shannon.txt** :  Alpha diversity (Shannon index) for all samples. Note, tables have been rarefied before calculation (1000 reads).

**observed_otus.txt** :  Alpha diversity (Observed OTUs) for all samples. Note, tables have been rarefied before calculation (1000 reads).

-----------------------------------
**"Description of sample names"**:

ALL sample names are build in the same way, e.g. spa3.RP.4S.V5.Y1

1. sample location, here **spa3** (see scripts/soil_props_apr10_2019.txt for details and alternative names)
2. compartment, here **RP** eq Rhizoplane, RS eq Rhizosphere, Soil & Root. 
3. sample type, here Single (**S**). Could also be pooled (P) or neighbouring plant (N)
4. primer, **V5** (bacteria), ITS1 (fungal) ITS1o oomycetal
5. Sampling year. **Y1** = 2015, Y2 = 2016, Y3 = 2017 
















