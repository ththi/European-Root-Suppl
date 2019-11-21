###
### this script need to be called separately for bacteria, fungi and oomycetes
###

#tax_mat=read.table("fun_avg_tax_mat.txt",header=T)

tax_mat=read.table("bacteria/bac_avg_tax_mat.txt",header=T)

#tax_mat=read.table("oo_avg_tax_mat.txt",header=T)



#avg
col_soil=grep("(Soil)",colnames(tax_mat))
col_root=grep("(Root)",colnames(tax_mat))


#avg

use_tax=rev(order(rowMeans(tax_mat[,col_soil])))




sub_tax=tax_mat[use_tax,col_soil]

sub_tax2=tax_mat[use_tax,col_root]

col_swe=grep("swe",colnames(sub_tax))
col_no_swe=grep("swe",colnames(sub_tax),invert=T)

par(mar=c(12,4,4,4),mfrow=c(1,2))

boxplot(t(sub_tax[,col_swe]),las=2,col="blue",boxwex=0.3,at=seq(0.5,length(use_tax),1),ylim=c(-0.03,1),outline=F,notch=F,border="black")
boxplot(t(sub_tax[,col_no_swe]),las=2,col="grey",add=T,boxwex=0.3,at=seq(0.85,length(use_tax),1),outline=F,notch=F,border="black",xaxt="n")

new_mat_mean=matrix(NA, nrow(sub_tax),4)
new_mat_mean[,1]=rownames(sub_tax)
new_mat_mean[,2]=rowMeans(sub_tax[,col_swe])
new_mat_mean[,3]=rowMeans(sub_tax[,col_no_swe])
change_vals=abs(as.numeric(new_mat_mean[,2])-as.numeric(new_mat_mean[,3]))
new_mat_mean[,4]=change_vals/(as.numeric(new_mat_mean[,2])/100)


new_mat=matrix(NA, nrow(sub_tax),4)

new_mat[,1]=rownames(sub_tax)

### swe vs rest

for(i in 1:nrow(sub_tax)){

		
	w_res=wilcox.test(as.numeric(sub_tax[i,col_swe]),as.numeric(sub_tax[i,col_no_swe]),alternative="two.sided")
	new_mat[i,2]=w_res$p.value	
	
}

ad_li=p.adjust(new_mat[,2],method="fdr",n=nrow(new_mat))
ad_pos=which(ad_li<=0.05)
points(ad_pos-0.25,rep(-0.025,length(ad_pos)),pch=8,cex=0.5)

col_swe2=grep("swe",colnames(sub_tax2))
col_no_swe2=grep("swe",colnames(sub_tax2),invert=T)

boxplot(t(tax_mat[use_tax,col_soil]),las=2,col="brown",boxwex=0.3,at=seq(0.5,length(use_tax),1),ylim=c(-0.03,1),outline=F,notch=F,border="black")
# number of groups +1 for last one (when first se argument is above one)
boxplot(t(tax_mat[use_tax,col_root]),las=2,col="dark green",add=T,boxwex=0.3,at=seq(0.85,length(use_tax),1),outline=F,notch=F,border="black",xaxt="n")



### root vs soil

for(i in 1:nrow(sub_tax)){


	w_res=wilcox.test(as.numeric(sub_tax[i,c(col_swe,col_no_swe)]),as.numeric(sub_tax2[i,c(col_swe2,col_no_swe2)]),alternative="two.sided")	
	new_mat[i,4]=w_res$p.value	

}
ad_li=p.adjust(new_mat[,4],method="fdr",n=nrow(new_mat))
ad_pos=which(ad_li<=0.05)
points(ad_pos-0.25,rep(-0.025,length(ad_pos)),pch=8,cex=0.5)



