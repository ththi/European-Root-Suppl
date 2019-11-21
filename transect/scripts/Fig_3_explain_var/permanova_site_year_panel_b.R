###
### nned to be called separately for each dataset (bacteria, fungi and oomycetes)
###

### need library("vegan")

library("vegan")

### Load data ###

mat.dist=read.table("../../oomycetal_data/bray_curtis_otu_table_norm.txt", h=T, sep="\t")
metadata=read.table("oomycete_table.txt", h=T)




### looping on the fractions to test the effect of site and year 
fraction<- c("Soil","RS","RP", "Root")
for (f in 1:4){
  fac=fraction[f]
  result<-matrix(ncol=4,nrow=5)
  
  tab<-paste("t", fac, sep="")
  
  pat<-paste(fac, "(1|2|3|4|)(P|S|)", sep = ".")
  aa_pat=grep(pat,rownames(mat.dist))
  aa_inter=intersect(rownames(mat.dist[aa_pat,]),rownames(metadata))
  
  res=adonis(formula= mat.dist[aa_inter,aa_inter] ~ as.vector(metadata[aa_inter,"Site"])*as.vector(metadata[aa_inter,"Experiment"]), data=as.data.frame(metadata[aa_inter,]), permutations = 9999)
  
  for (e in 1:5){
    
    result[e,1]<-res$aov.tab$Df[e]
    result[e,2]<-res$aov.tab$F.Model[e]
    result[e,3]<-res$aov.tab$`Pr(>F)`[e]
    result[e,4]<-res$aov.tab$R2[e]
  }
  assign(tab, result)
}

result_bac=cbind(tSoil,tRS,tRP,tRoot)
colnames(result_bac)=c("Soil_Df","Soil_F","Soil_Pr","Soil_R2","RS_Df", "RS_F","RS_Pr", "RS_R2", "RP_Df","RP_F", "RP_Pr","RP_R2","Root_Df","Root_F","Root_Pr","Root_R2")
row.names(result_bac)=c("Site", "Year", "Site:Year", "Residuals", "Total")
result_bac

###
###write.table(result_bac,"year_site_perFraction_oomyc_zotus.txt")

