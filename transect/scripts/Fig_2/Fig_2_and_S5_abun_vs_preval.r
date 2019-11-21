###
### this script needs to be called separately for bacteria_fungi and oomycetes
###

tab=read.table("../../bacterial_data/otu_tab_filter.txt",header=T)

colnames(tab)=gsub("Pul","pul",colnames(tab))

tax=read.table("../../bacterial_data/taxonomy_otus.txt",header=F,fill=T)

high_read=colSums(tab)>=1000
tab=tab[,high_read]

## only for bacterial data
##

rem_1=grep("Chloroplast",tax[,4])
rem_2=grep("mitochondria",tax[,6])
rem=c(rem_1,rem_2)
rem_otu=tax[rem,1]
keep_otu=setdiff(rownames(tab),rem_otu)

tab=tab[keep_otu,]

##
##
	
tab_ra=sweep(tab,2,colSums(tab),"/")
tab_c=tab_ra>0.001    #### important cutoff, when is something counted as present ?! orig 0.001


tab_l=log2(tab_ra)
tab_l[is.infinite(as.matrix(tab_l))]<-0


dev.new(width=3.68,height=11.8)
par(mfrow=c(4,1),mar=c(5,5,3,3))


otu_def=matrix(data=4,nrow=nrow(tab_ra),ncol=5,dimnames=list(rownames(tab_ra),c("soil","rs","rp","root","root_nei")))

otu_def_preval=matrix(data=0,nrow=nrow(tab_ra),ncol=4,dimnames=list(rownames(tab_ra),c("soil","rs","rp","root")))
otu_def_num_sites=matrix(data=0,nrow=nrow(tab_ra),ncol=6,dimnames=list(rownames(tab_ra),c("soil","rs","rp","root","root_n","root_y1")))
otu_def_avg_ra_log=matrix(data=0,nrow=nrow(tab_ra),ncol=4,dimnames=list(rownames(tab_ra),c("soil","rs","rp","root")))
otu_def_avg_ra=matrix(data=0,nrow=nrow(tab_ra),ncol=4,dimnames=list(rownames(tab_ra),c("soil","rs","rp","root")))
otu_def_year_pres=matrix(data=0,nrow=nrow(tab_ra),ncol=4,dimnames=list(rownames(tab_ra),c("soil","rs","rp","root")))

ra_dat=""

for(ii in 1:4){
	if(ii ==1){	
		col_1=grep("Soil.*(P|S).(ITS1|V5|ITS1o).Y1",colnames(tab))
		col_2=grep("Soil.*(P|S).(ITS1|V5|ITS1o).Y2",colnames(tab))
		col_3=grep("Soil.*(P|S).(ITS1|V5|ITS1o).Y3",colnames(tab))
		col_nei=grep("Soil.N.*.Y1",colnames(tab))
		plot_nam="Soil"
	}
	if(ii ==2){	
		col_1=grep("RS.*(P|S).(ITS1|V5|ITS1o).Y1",colnames(tab))
		col_2=grep("RS.*(P|S).(ITS1|V5|ITS1o).Y2",colnames(tab))
		col_3=grep("RS.*(P|S).(ITS1|V5|ITS1o).Y3",colnames(tab))
		col_nei=grep("RS.N.*.Y1",colnames(tab))
		plot_nam="RS"
	}
	if(ii ==3){	
		col_1=grep("RP.*(P|S).(ITS1|V5|ITS1o).Y1",colnames(tab))
		col_2=grep("RP.*(P|S).(ITS1|V5|ITS1o).Y2",colnames(tab))
		col_3=grep("RP.*(P|S).(ITS1|V5|ITS1o).Y3",colnames(tab))
		col_nei=grep("RP.N.*.Y1",colnames(tab))
		plot_nam="RP"
		
	}
	if(ii ==4){	
		col_1=grep("Root.*(P|S).(ITS1|V5|ITS1o).Y1",colnames(tab))
		col_2=grep("Root.*(P|S).(ITS1|V5|ITS1o).Y2",colnames(tab))
		col_3=grep("Root.*(P|S).(ITS1|V5|ITS1o).Y3",colnames(tab))
		col_nei=grep("Root.N.*.Y1",colnames(tab))
		plot_nam="Root"
		
	}
		


	plot(0,0.1,cex=0.1,xlim=c(0,1),ylim=c(-10,0),ylab="log2 (RA)",xlab="Prevalence")
	text(0.2,-0.5,plot_nam,pos=2)

	###
	aa_thres=mean(log2(0.001))
	

	aa_gen=0
	aa_spec=0
	aa_local=0
	aa_couy1=0
	aa_couy2=0
	aa_couy3=0
	list_gen=""
	list_spec=""
	list_local=""
	list_nei=""		
	
	### get number of max samples
			
	num_si_1=length(unique(sub("(Soil|RS|RP|Root).*","",colnames(tab_c[,col_1]))))
	num_si_2=length(unique(sub("(Soil|RS|RP|Root).*","",colnames(tab_c[,col_2]))))
	num_si_3=length(unique(sub("(Soil|RS|RP|Root).*","",colnames(tab_c[,col_3]))))

	num_si_nei=length(unique(sub("(Soil|RS|RP|Root).*","",colnames(tab_c[,col_nei]))))
	
	for(i in 1:nrow(tab_ra)){
	
		if(sum(tab_c[i,c(col_1,col_2,col_3)])>5){   #### min 5 samples present !!!!		
			
		
			aa_ra=rowSums(tab_ra[i,c(col_1,col_2,col_3)])/sum(tab_c[i,c(col_1,col_2,col_3)])

			aa_prev=sum(tab_c[i,c(col_1,col_2,col_3)])/ncol(tab_c[,c(col_1,col_2,col_3)])
			
			aa_ra_nei=rowSums(tab_ra[i,c(col_nei)]/sum(tab_c[i,c(col_nei)]))			

			aa_prev_nei=sum(tab_c[i,c(col_nei)])/ncol(tab_c[,c(col_nei)])			


			# prev for all samples	

			aa_prev1=sum(tab_c[i,col_1])/ncol(tab_c[,col_1])
			aa_prev2=sum(tab_c[i,col_2])/ncol(tab_c[,col_2])
			aa_prev3=sum(tab_c[i,col_3])/ncol(tab_c[,col_3])

			aa_size=sum(aa_prev1>0,aa_prev2>0,aa_prev3>0)
		
			mean_prev=sum(c(aa_prev1,aa_prev2,aa_prev3))/aa_size

			
			### get only thos samples where OTU is present using which

			pres_1=which(tab_c[i,col_1]==TRUE)
			pres_2=which(tab_c[i,col_2]==TRUE)
			pres_3=which(tab_c[i,col_3]==TRUE)

			pres_col_nei=which(tab_c[i,col_nei]==TRUE)
			
			
			### get substring, eg. remove(!) everything that is indicated by sub (keep only site name)

			num_pres_si_1=length(unique(sub("(Soil|RS|RP|Root).*","",colnames(tab_c[,col_1])[pres_1] )))
			num_pres_si_2=length(unique(sub("(Soil|RS|RP|Root).*","",colnames(tab_c[,col_2])[pres_2] )))
			num_pres_si_3=length(unique(sub("(Soil|RS|RP|Root).*","",colnames(tab_c[,col_3])[pres_3] )))



			mean_prev_si=((num_pres_si_1/num_si_1)+(num_pres_si_2/num_si_2)+(num_pres_si_3/num_si_3))/aa_size

			otu_def_preval[i,ii]=mean_prev
			otu_def_num_sites[i,ii]=mean_prev_si
			otu_def_avg_ra_log[i,ii]=log2(aa_ra)
			otu_def_avg_ra[i,ii]=aa_ra
			otu_def_year_pres[i,ii]=aa_size

			### for N OTUs
			
			num_pres_nei=length(unique(sub("(Soil|RS|RP|Root).*","",colnames(tab_c[,col_nei[pres_col_nei]]))))
			mean_prev_nei=(num_pres_nei/num_si_nei)	

			if(ii==4){otu_def_num_sites[i,5]=mean_prev_nei}
			if(ii==4){otu_def_num_sites[i,6]=num_pres_si_1/num_si_1}
		
			if(mean_prev_si >= 0.8 & log2(aa_ra)>= aa_thres){
				acolor="red"
				#aa_gen=aa_gen+1
				aa_gen=aa_gen+aa_ra
				list_gen=c(list_gen,rownames(tab_ra[i,]))
				otu_def[i,ii]=1
			}
			else if(mean_prev_si > 0.2 & mean_prev < 0.8 & log2(aa_ra)>= aa_thres){
				acolor="orange"
				#aa_spec=aa_spec+1
				aa_spec=aa_spec+aa_ra
				list_spec=c(list_spec,rownames(tab_ra[i,]))
				otu_def[i,ii]=2
			}
			else if(mean_prev_si <= 0.2  & log2(aa_ra)>= aa_thres){
				acolor="dark blue"
				#aa_local=aa_local+1
				aa_local=aa_local+aa_ra
				list_local=c(list_local,rownames(tab_ra[i,]))
				otu_def[i,ii]=3
			}		
			else{
				acolor="grey"
				otu_def[i,ii]=4
			
			}
			if(aa_size==1){
				points(mean_prev_si,log2(aa_ra),pch=3,col=acolor)
				aa_couy1=aa_couy1+1
					
			}
			if(aa_size==2){
				points(mean_prev_si,log2(aa_ra),pch=2,col=acolor)
				aa_couy2=aa_couy2+1
				
			}
			if(aa_size==3){
				points(mean_prev_si,log2(aa_ra),pch=1,col=acolor)
				aa_couy3=aa_couy3+1
				
				
			}
			
			# mark those who are bigger 0.8 preval in Neighb. samples (not used in publication)	
			if(mean_prev_nei >0.8){ 

				list_nei=c(list_nei,rownames(tab_ra[i,]))
				if(ii==4){otu_def[i,5]=1}
			}
			if(mean_prev_nei <0.8 & mean_prev_nei >0.2){ 
				
				if(ii==4){otu_def[i,5]=2}
				
			}
			if( mean_prev_nei <0.2){ 
				
				if(ii==4){otu_def[i,5]=3}
				
			}

			
		}
		

	}
	
	### correlation / lm

	not_zero=otu_def_num_sites[,ii]!=0
	#plot(otu_def_num_sites[not_zero,ii],otu_def_avg_ra_log[not_zero,ii])	
	abline(lm(otu_def_avg_ra_log[not_zero,ii]~otu_def_num_sites[not_zero,ii]))
	lm_val=lm(otu_def_avg_ra_log[not_zero,ii]~otu_def_num_sites[not_zero,ii])
	lm_pval=summary(lm_val)$coefficients[2,4]
	lm_adr=summary(lm_val)$adj.r.squared
	text(1,-1,paste("r2",lm_adr,sep="\t"),pos=2)
	text(1,-2,paste("pval",lm_pval,sep="\t"),pos=2)


}



dev.new(width=4,height=4.2)
aa_cou=0	
for(ii in 1:4){
		
		#if(i==1){aa_col="red"}
		#if(i==2){aa_col="orange"}
		#if(i==3){aa_col="dark blue"}


	for(i in 1:3){

		if(i==1){aa_col="dark blue";aa_group=3}
		if(i==2){aa_col="orange";aa_group=2}
		if(i==3){aa_col="red";aa_group=1}

		aa_cou=aa_cou+1
		
		aa=otu_def[,ii]==aa_group
		aa_logs=otu_def_avg_ra_log[aa,ii]
		if(i == 1 & ii == 1){	
			boxplot(aa_logs,ylim=c(-10,0),col=aa_col,medcol="white",at=aa_cou,xlim=c(0.5,12.5),boxwex=1.2,outwex=0.001,outline=F)
		}else{
			boxplot(aa_logs,col=aa_col,medcol="white",at=aa_cou,add=T,boxwex=1.2,outwex=0.001,outline=F)
		}		

	}


}
axis(1,at=c(2,5,8,11),labels=c("Soil","RS","RP","Root"),las=2)

### group stats

#aa_soil=otu_def[,1]==3
#aa_rs=otu_def[,2]==3
#aa_rp=otu_def[,3]==3
#aa_root=otu_def[,4]==3	
#aa_all=c(otu_def_avg_ra_log[aa_soil,1],otu_def_avg_ra_log[aa_rs,2],otu_def_avg_ra_log[aa_rp,3],otu_def_avg_ra_log[aa_root,4])

### library("dunn.test")

#kruskal.test(aa_all,c(rep(1,sum(aa_soil)),rep(2,sum(aa_rs)),rep(3,sum(aa_rp)),rep(4,sum(aa_root))))
#dunn.test(aa_all,c(rep(1,sum(aa_soil)),rep(2,sum(aa_rs)),rep(3,sum(aa_rp)),rep(4,sum(aa_root))),method="bh")

#or
#pairwise.wilcox.test(aa_all,c(rep(1,sum(aa_soil)),rep(2,sum(aa_rs)),rep(3,sum(aa_rp)),rep(4,sum(aa_root))),p.adjust.method="BH")


