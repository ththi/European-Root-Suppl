library("vegan")

tab_b=read.table("../../bacterial_data/bray_curtis_otu_table_norm.txt", sep="\t", header=T, check.names=F)


colnames(tab_b)=gsub("Pul","pul",colnames(tab_b))
rownames(tab_b)=gsub("Pul","pul",rownames(tab_b))

soilp=read.table("../soil_props_apr10_2019.txt",header=T)

dump=soilp
dump=cbind(dump,"year"=dump[,11])

####
#### this need to be changed according to input data, v5 for bacteria, ITS1 for fungi and ITS1o for oomycetes
####

rownames(dump)=gsub("XXX","V5",rownames(dump))
#rownames(dump)=gsub("XXX","ITS1",rownames(dump))
#rownames(dump)=gsub("XXX","ITS1o",rownames(dump))

dev.new(height=5.044,width=8.093)
par(mfrow=c(1,4))
plo_val=1

df_vals=c("x")
t_vals=c("x")


for(i in 1:8){

	if(i==1){pat="(Soil).(1|2|3|4)(P|S).*Y1"; plot_nam="Year 1"}
	if(i==2){pat="(Root).(1|2|3|4)(P|S).*Y1"}
	if(i==3){pat="(Soil).(1|2|3|4)(P|S).*Y2"; plot_nam="Year 2"}
	if(i==4){pat="(Root).(1|2|3|4)(P|S).*Y2"}
	if(i==5){pat="(Soil).(1|2|3|4)(P|S).*Y3"; plot_nam="Year 3"}
	if(i==6){pat="(Root).(1|2|3|4)(P|S).*Y3"}
	if(i==7){pat="(Soil).(1|2|3|4)(P|S).*Y"; plot_nam="All years"}
	if(i==8){pat="(Root).(1|2|3|4)(P|S).*Y"}

	aa_pat=grep(pat,rownames(tab_b))
	aa_inter=intersect(rownames(tab_b[aa_pat,]),rownames(dump))
	
	#site col 5 , fraction col 8

	res=betadisper(as.dist(tab_b[aa_inter,aa_inter]),group=dump[aa_inter,5])
	#plot(res)

	perm_res=permutest(res)	

	pval=perm_res$tab[1,6]

	res_dist=dist(res$centroids[,1:2])
	
	## permutest(res)

	nam=paste("res_dist_",i,sep="")
	assign(nam,res_dist)


	if(plo_val==1 && i <=6){
		boxplot(res_dist,xlim=c(0.5,2.5),ylim=c(0,1),main=plot_nam,col="brown",width=1.5)
		plo_val=2
		res_last=res_dist
	}
	else{
		if(i <= 6){
			boxplot(res_dist,at=2,add=T,col="dark green",width=1.5)
			plo_val=1
			t_res=t.test(res_dist,res_last)
			text(0.5,1,t_res$p.value,pos=4,cex=0.7)
		
			df_vals=c(df_vals,t_res$parameter)
			t_vals=c(t_vals,t_res$statistic)
		}
		
	}

}


### for figure done via t-test,

all_soil=c(as.vector(res_dist_1),as.vector(res_dist_3),as.vector(res_dist_5))
all_root=c(as.vector(res_dist_2),as.vector(res_dist_4),as.vector(res_dist_6))

boxplot(all_soil,xlim=c(0.5,2.5),ylim=c(0,1),main=plot_nam,col="brown",width=1.5)
boxplot(all_root,at=2,add=T,col="dark green",width=1.5)

t_res=t.test(all_soil,all_root)
text(0.5,1,t_res$p.value,pos=4,cex=0.7)



