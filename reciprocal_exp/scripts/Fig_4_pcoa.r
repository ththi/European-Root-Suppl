### This script recapitulates Figure 4, panel b. It is not the original script, therefore axis might be mirrored.
### The output includ one figure with data for bacteria, fungi and oomycetes
###

tab_b=read.table("../bacterial_data/bray_curtis_otu_table_norm.txt", sep="\t", header=T, check.names=F)
tab_f=read.table("../fungal_data/bray_curtis_otu_table_norm.txt", sep="\t", header=T, check.names=F)
tab_o=read.table("../oomycetal_data/bray_curtis_otu_table_norm.txt", sep="\t", header=T, check.names=F)


cou=0
dev.new(width=12.9,height=6.9)
par(mfrow=c(2,3),mar=c(8,8,2,8),xpd=T)



for(i in 1:2){


	aa_pat="(Italy|Sweden).(Root|Soil|DP.Soil).(It|Sw)"

	for(ii in 1:3){

		if(ii==1){
			bc_bac=tab_b

			plot_nam="Bacteria"	
		}
		if(ii==2){ 
			bc_bac=tab_f

			plot_nam="Fungi"
		}
		if(ii==3){ 
			bc_bac=tab_o

			plot_nam="Oomycetes"
		}

		aa=grep(aa_pat,colnames(bc_bac))
		
		if(ii==1){bac_num_sam=length(aa)}
		if(ii==2){fun_num_sam=length(aa)}
		if(ii==3){oo_num_sam=length(aa)}	

		bc_bac=bc_bac[aa,aa]

		pcoa_bac <- cmdscale(bc_bac, k=2, eig=T)
		bac_v1=format(100*pcoa_bac$eig[1]/sum(pcoa_bac$eig),digits=4)
		bac_v2=format(100*pcoa_bac$eig[2]/sum(pcoa_bac$eig),digits=4)

		plot(pcoa_bac$points[,1],pcoa_bac$points[,2],cex=0.5,col="white",xlab=bac_v1,ylab=bac_v2,main=plot_nam)

		

		if(i==1){

			aa_ro_it=grep("(Sweden|Italy).Root.(It|Sw).It",row.names(bc_bac),perl=T)
			aa_ro_sw=grep("(Sweden|Italy).Root.(It|Sw).Sw",row.names(bc_bac),perl=T)
			aa_so_it=grep("(Sweden|Italy).Soil.(It|Sw).It",row.names(bc_bac),perl=T)
			aa_so_sw=grep("(Sweden|Italy).Soil.(It|Sw).Sw",row.names(bc_bac),perl=T)
			


			points(pcoa_bac$points[aa_ro_it,1],pcoa_bac$points[aa_ro_it,2],col="red",pch=0,cex=1.5)
			points(pcoa_bac$points[aa_ro_sw,1],pcoa_bac$points[aa_ro_sw,2],col="blue",pch=0,cex=1.5)
			points(pcoa_bac$points[aa_so_it,1],pcoa_bac$points[aa_so_it,2],col="red",pch=15,cex=1.5)
			points(pcoa_bac$points[aa_so_sw,1],pcoa_bac$points[aa_so_sw,2],col="blue",pch=15,cex=1.5)
			

			if(ii==3){
				legend(0.5,0.2,c("It1 plants","Sw4 plants"),fill=c("red","blue"))
				legend(0.5,-0.1,c("Soil","Whole Root"),pch=c(15,0))
				
			}
		
			
		}

		if(i==2){
			
			aa_sw_it=grep("Sweden.(Soil|Root).It.",row.names(bc_bac),perl=T)
			aa_sw_sw=grep("Sweden.(Soil|Root).Sw.",row.names(bc_bac),perl=T)
			aa_it_it=grep("Italy.(Soil|Root).It.",row.names(bc_bac),perl=T)
			aa_it_sw=grep("Italy.(Soil|Root).Sw.",row.names(bc_bac),perl=T)
			


			points(pcoa_bac$points[aa_sw_it,1],pcoa_bac$points[aa_sw_it,2],col="black",pch=1,cex=1.5)
			points(pcoa_bac$points[aa_sw_sw,1],pcoa_bac$points[aa_sw_sw,2],col="black",pch=2,cex=1.5)
			points(pcoa_bac$points[aa_it_it,1],pcoa_bac$points[aa_it_it,2],col="black",pch=19,cex=1.5)
			points(pcoa_bac$points[aa_it_sw,1],pcoa_bac$points[aa_it_sw,2],col="black",pch=17,cex=1.5)

			if(ii==3){
				legend(0.5,0.2,c("It1 location","Sw4 location"),fill=c("black","white"))
				legend(0.5,-0.1,c("It1 soil","Sw4 soil"),pch=c(1,2))
			}
			
		}
		

		
	

	}# for ii 1:3 bac fun oo




}# i 1:4 soil, rs, rp ro



