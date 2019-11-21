This folder includes data and scripts used for generating Figure 4 of the manuscript (and respective Supp. Figures).

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

ALL sample names are build in the same way, e.g. Italy.Root.It.It15.1.V5

1. sample location, here **Italy** otherwise Seden
2. compartment, here **Root** otherwise Soil 
3. soil origin, here Italy (**It**) otherwise sweden (Sw)
4. plant genotyoe, here Italy (**It**) otherwise sweden (Sw), plus sample identifier
5. primer **V5** = bacteria, ITS1 = fungi,  ITS1o = oomycetes 
















