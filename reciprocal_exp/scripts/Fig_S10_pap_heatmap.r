### This script recapitulates Figures as seen in figure S10.
### It needs to be called separately for bacterial, fungal and oomycetal data.
### Therefore the name/path of the otu_table and taxonomy file need to be changed
### The heatmap has no description for its columns, except samplenames.
### Please see the publication for information of the colorcode.


###
### needs library("gplots")
###

library("gplots")

tab=read.table("../fungal_data/otu_tab_filter.txt",header=T)

colnames(tab)=gsub("Pul","pul",colnames(tab))

# !
colnames(tab)=gsub("DP.Soil","Soil",colnames(tab))

tax=read.table("../fungal_data/taxonomy_otus.txt",header=F,fill=T)

### really strange sample
tab=tab[,-38]

high_read=colSums(tab)>=1000
tab=tab[,high_read]


### for bacterial data
rem_1=grep("Chloroplast",tax[,4])
rem_2=grep("mitochondria",tax[,6])

rem=c(rem_1,rem_2)
rem_otu=tax[rem,1]
keep_otu=setdiff(rownames(tab),rem_otu)

tab=tab[keep_otu,]

tab_ra=sweep(tab,2,colSums(tab),"/")

col_1=grep("(Italy|Sweden).Root.(It|Sw)",colnames(tab_ra))


aa_thres=0.0001
xx=rowMeans(tab_ra[,col_1])>=aa_thres

tab_ra_red=tab_ra[xx,col_1]

tab_pap=tab_ra_red

tab_pap[tab_pap >=0.0001 ] <- 1


# do log transform

tab_l=log2(tab_ra_red)

# how to handle abscent samples ??? for later it is important to have really low RA and not zeros ! 
tab_l[is.infinite(as.matrix(tab_l))]<- log2(0.00001)


my_palette <- colorRampPalette(c("#033d9b","grey","yellow","orange","red3","magenta","magenta","red3","red3","red3","red3","red3"))(n = 200)	


col_it_it=grep("Italy.(Root|Soil).It.(It|Sw)",colnames(tab_pap))
col_it_sw=grep("Italy.(Root|Soil).Sw.(It|Sw)",colnames(tab_pap))

col_sw_sw=grep("Sweden.(Root|Soil).Sw.(It|Sw)",colnames(tab_pap))
col_sw_it=grep("Sweden.(Root|Soil).It.(It|Sw)",colnames(tab_pap))
	

it_soil=c(col_it_it,col_sw_it)
sw_soil=c(col_it_sw,col_sw_sw)

it_clim=c(col_it_it,col_it_sw)
sw_clim=c(col_sw_sw,col_sw_it)

new_res=cbind(rowMeans(tab_l[,col_it_it]),rowMeans(tab_l[,col_it_sw]) ,rowMeans(tab_l[,col_sw_it]) ,rowMeans(tab_l[,col_sw_sw]) )

test_val=rowMeans(new_res)
aa1=new_res[,1]>test_val & new_res[,3]>test_val & new_res[,2]<test_val & new_res[,4]<test_val
it_soil_enrich=rownames(new_res[aa1,])

aa2=new_res[,2]>test_val & new_res[,4]>test_val & new_res[,1]<test_val & new_res[,3]<test_val
sw_soil_enrich=rownames(new_res[aa2,])

aa3=new_res[,3]>test_val & new_res[,4]>test_val & new_res[,1]<test_val & new_res[,2]<test_val
sw_clim_enrich=rownames(new_res[aa3,])

aa4=new_res[,1]>test_val & new_res[,2]>test_val & new_res[,3]<test_val & new_res[,4]<test_val
it_clim_enrich=rownames(new_res[aa4,])

aa5=new_res[,1]>test_val & new_res[,2]<test_val & new_res[,3]<test_val & new_res[,4]<test_val
it_it_enrich=rownames(new_res[aa5,])

aa6=new_res[,4]>test_val & new_res[,1]<test_val & new_res[,2]<test_val & new_res[,3]<test_val
sw_sw_enrich=rownames(new_res[aa6,])

aa7=test_val> -9.96
univ_enrich=rownames(new_res[aa7,])


col_def=matrix("white",nrow=ncol(tab_pap),ncol=2,dimnames=list(colnames(tab_pap) ,c("def","dump") ))
xx=grep("(Sw|It).It",colnames(tab_pap))
col_def[xx,1]="red"
xx=grep("(Sw|It).Sw",colnames(tab_pap))
col_def[xx,1]="blue"


otu_def=matrix("white",nrow=nrow(tab_pap),ncol=3,dimnames=list(rownames(tab_pap) ,c("def","tax","wide") ))


otu_def[univ_enrich,1]="black"
otu_def[it_soil_enrich,1]="red"
otu_def[sw_soil_enrich,1]="blue"
otu_def[sw_clim_enrich,1]="light blue"
otu_def[it_clim_enrich,1]="pink"
otu_def[it_it_enrich,1]="red4"
otu_def[sw_sw_enrich,1]="dark blue"


tax_match=match(rownames(tab_pap),tax[,1])
otu_def[,2]=as.vector(tax[tax_match,4])


no_white=grep("white",otu_def[,1],invert=T)


dev.new(width=7,height=11.8)
heat_res=heatmap.2(as.matrix(tab_l[no_white,]),tracecol=F,mar=c(9,5),RowSideColors=otu_def[no_white,1],col=my_palette)



###

all_nam=names(table(otu_def[no_white,2]))
dev.new(width=10.41,height=11)
par(mfrow=c(4,5))

for(i in 1:12){

	aa_nam=names(rev(sort(table(otu_def[no_white,2])))[i])
	aa_match=grep(aa_nam,otu_def[no_white,2])
	aa_tab=table(otu_def[no_white[aa_match],1])
	leg_nam=paste(aa_nam,sum(aa_tab),sep=" ")

	pie(aa_tab,main=leg_nam,col=names(aa_tab))	
	
}

pie(table(otu_def[no_white,1]),main=paste("all:",sum(table(otu_def[no_white,1])),sep=" "),col= names(table(otu_def[no_white,1])))



# 

