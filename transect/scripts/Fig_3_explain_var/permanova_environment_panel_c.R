###
### nned to be called separately for each dataset (bacteria, fungi and oomycetes)
###

### need library("vegan")

library("vegan")


### Load data ###
mat.dist=read.table("../Fig1_dendogram_S2_barplots/bacteria/bac_avg_dist.txt", h=T, sep="\t", row.names = 1)
metadata=read.table("metadata_sites_no_years.txt", h=T)

# select samples
pat=".Y"

# find those samples in the distance matrix (tab_b)
aa_pat=grep(pat,rownames(mat.dist))

# find those samples in the propertie table (dump), by comparing rownames
aa_inter=intersect(rownames(mat.dist[aa_pat,]),rownames(metadata))

### loop into the fractions to test separately each fraction ###
fraction<- c("Soil","RS","RP", "Root")
for (f in 1:4){
  fac=fraction[f]
  result<-matrix(ncol=4,nrow=9)
  
  tab<-paste("t", fac, sep="")
  
  pat<-paste(fac, sep = ".")
  aa_pat=grep(pat,rownames(mat.dist))
  aa_inter=intersect(rownames(mat.dist[aa_pat,]),rownames(metadata))
  
  res=adonis(formula= mat.dist[aa_inter,aa_inter] ~ as.vector(metadata[aa_inter,"ph"])+ as.vector(metadata[aa_inter,"AvailableNO3"])+ as.vector(metadata[aa_inter,"Boron"])+ as.vector(metadata[aa_inter,"AvailableK"])+ as.vector(metadata[aa_inter,"Copper"])+ as.vector(metadata[aa_inter,"Manganesium"])+as.vector(metadata[aa_inter,"ReserveP"]) , data=as.data.frame(metadata[aa_inter,]))
  
  for (e in 1:9){
    
    result[e,1]<-res$aov.tab$Df[e]
    result[e,2]<-res$aov.tab$F.Model[e]
    result[e,3]<-res$aov.tab$`Pr(>F)`[e]
    result[e,4]<-res$aov.tab$R2[e]
  }
  assign(tab, result)
}

result_bac=cbind(tSoil,tRS,tRP,tRoot)

colnames(result_bac)=c("Soil_Df","Soil_F","Soil_Pr","Soil_R2","RS_Df", "RS_F","RS_Pr", "RS_R2", "RP_Df","RP_F", "RP_Pr","RP_R2","Root_Df","Root_F","Root_Pr","Root_R2")
row.names(result_bac)=c("pH", "AvailableNO3", "Boron", "AvailableK","Copper","Manganesium", "ReserveP","Residuals", "Total")
result_bac

### if you want to print results to a table
# write.table(result_bac,"results_env_bac_noyears.txt")



