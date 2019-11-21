

tab_b=read.table("../bacterial_data/bray_curtis_otu_table_norm.txt", sep="\t", header=T, check.names=F)
tab_f=read.table("../fungal_data/bray_curtis_otu_table_norm.txt", sep="\t", header=T, check.names=F)
tab_o=read.table("../oomycetal_data/bray_curtis_otu_table_norm.txt", sep="\t", header=T, check.names=F)


cou=0
dev.new(width=11,height=8)
par(mfrow=c(3,3),mar=c(3,3,2,10),xpd=T)



for(i in 1:3){


	aa_pat="(Soil|RS|RP|Root).(1|2|3|4)(P|S)"

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

			aa_spa1=grep("(spa1).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_spa2=grep("(spa2).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_spa3=grep("(spa3).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)

			aa_fra1=grep("(SD).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_fra2=grep("(tou1).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_fra3=grep("(tou2).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)

			aa_ger1=grep("(tub2).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_ger2=grep("(tub3).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_ger3=grep("(tub4).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_ger4=grep("(Pulh|pulh).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_ger5=grep("geyen.(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_ger6=grep("ham.(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)

			aa_swe1=grep("(swe1).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_swe2=grep("(swe2).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_swe3=grep("(swe3).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)
			aa_swe4=grep("(swe4).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)

			aa_ita=grep("(Ita).(Soil|RS|Root|RP).*Y",row.names(bc_bac),perl=T)

			point_val=1.5		

			points(pcoa_bac$points[aa_spa1,1],pcoa_bac$points[aa_spa1,2],col="tomato4",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_spa2,1],pcoa_bac$points[aa_spa2,2],col="tomato1",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_spa3,1],pcoa_bac$points[aa_spa3,2],col="firebrick3",pch=19,lwd=0.75,cex=point_val,bg="white")

			points(pcoa_bac$points[aa_fra1,1],pcoa_bac$points[aa_fra1,2],col="yellow",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_fra2,1],pcoa_bac$points[aa_fra2,2],col="orange",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_fra3,1],pcoa_bac$points[aa_fra3,2],col="chocolate",pch=19,lwd=0.75,cex=point_val,bg="white")


			points(pcoa_bac$points[aa_ger1,1],pcoa_bac$points[aa_ger1,2],col="green",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_ger2,1],pcoa_bac$points[aa_ger2,2],col="lightcyan",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_ger3,1],pcoa_bac$points[aa_ger3,2],col="palegreen2",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_ger4,1],pcoa_bac$points[aa_ger4,2],col="olivedrab1",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_ger5,1],pcoa_bac$points[aa_ger5,2],col="seagreen",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_ger6,1],pcoa_bac$points[aa_ger6,2],col="darkslategrey",pch=19,lwd=0.75,cex=point_val,bg="white")

			points(pcoa_bac$points[aa_swe1,1],pcoa_bac$points[aa_swe1,2],col="skyblue",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_swe2,1],pcoa_bac$points[aa_swe2,2],col="deepskyblue",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_swe3,1],pcoa_bac$points[aa_swe3,2],col="royalblue1",pch=19,lwd=0.75,cex=point_val,bg="white")
			points(pcoa_bac$points[aa_swe4,1],pcoa_bac$points[aa_swe4,2],col="royalblue4",pch=19,lwd=0.75,cex=point_val,bg="white")

			points(pcoa_bac$points[aa_ita,1],pcoa_bac$points[aa_ita,2],col="red",pch=19,lwd=0.75,cex=point_val,bg="white")
			
			if(ii==3){
				legend(0.5,0.45,c("spa1","spa2","spa3","fr1","fr2","fr3","ger1","ger2","ger3","ger4","ger5","ger6","swe1","swe2","swe3","swe4","ita"),fill=c("tomato4","tomato1","firebrick3","yellow","orange","chocolate3","green","lightcyan","palegreen2","olivedrab1","seagreen","darkslategrey","skyblue","deepskyblue","royalblue1","royalblue4","red"))
				
			}
		
			
		}

		if(i==2){
			aa_1=grep("Y1",row.names(bc_bac),perl=T)
			aa_2=grep("Y2",row.names(bc_bac),perl=T)
			aa_3=grep("Y3",row.names(bc_bac),perl=T)

			points(pcoa_bac$points[aa_1,1],pcoa_bac$points[aa_1,2],col="black",pch=19,lwd=0.75,cex=point_val,bg="white")

			points(pcoa_bac$points[aa_2,1],pcoa_bac$points[aa_2,2],col="grey39",pch=19,lwd=0.75,cex=point_val,bg="white")

			points(pcoa_bac$points[aa_3,1],pcoa_bac$points[aa_3,2],col="grey70",pch=19,lwd=0.75,cex=point_val,bg="white")

			if(ii==3){
				legend(0.5,0.45,c("Year 1","Year 2","Year 3"),fill=c("black","grey39","grey70"))
				
			}
			
		}
		

		if(i==3){
			aa_1=grep("Soil",row.names(bc_bac),perl=T)
			aa_2=grep("RS",row.names(bc_bac),perl=T)
			aa_3=grep("RP",row.names(bc_bac),perl=T)
			aa_4=grep("Root",row.names(bc_bac),perl=T)

			aa_5=grep("swe(1|2|3|4).Soil",row.names(bc_bac),perl=T)
			aa_6=grep("swe(1|2|3|4).RS",row.names(bc_bac),perl=T)
			aa_7=grep("swe(1|2|3|4).RP",row.names(bc_bac),perl=T)
			aa_8=grep("swe(1|2|3|4).Root",row.names(bc_bac),perl=T)

			aa_1_no_swe=setdiff(aa_1,aa_5)
			aa_2_no_swe=setdiff(aa_2,aa_6)
			aa_3_no_swe=setdiff(aa_3,aa_7)
			aa_4_no_swe=setdiff(aa_4,aa_8)

			points(pcoa_bac$points[aa_1_no_swe,1],pcoa_bac$points[aa_1_no_swe,2],bg="brown",pch=21,lwd=0.75,cex=point_val,col="brown")

			points(pcoa_bac$points[aa_2_no_swe,1],pcoa_bac$points[aa_2_no_swe,2],bg="grey",pch=21,lwd=0.75,cex=point_val,col="grey")

			points(pcoa_bac$points[aa_3_no_swe,1],pcoa_bac$points[aa_3_no_swe,2],bg="light green",pch=21,lwd=0.75,cex=point_val,col="light green")

			points(pcoa_bac$points[aa_4_no_swe,1],pcoa_bac$points[aa_4_no_swe,2],bg="dark green",pch=21,lwd=0.75,cex=point_val,col="dark green")

			####

			points(pcoa_bac$points[aa_5,1],pcoa_bac$points[aa_5,2],bg="white",pch=21,lwd=1.5,cex=point_val,col="brown")

			points(pcoa_bac$points[aa_6,1],pcoa_bac$points[aa_6,2],bg="white",pch=21,lwd=1.5,cex=point_val,col="grey")

			points(pcoa_bac$points[aa_7,1],pcoa_bac$points[aa_7,2],bg="white",pch=21,lwd=1.5,cex=point_val,col="light green")

			points(pcoa_bac$points[aa_8,1],pcoa_bac$points[aa_8,2],bg="white",pch=21,lwd=1.5,cex=point_val,col="dark green")
			

			if(ii==3){
				legend(0.5,0.45,c("Soil","RS","RP","Root"),fill=c("brown","grey","light green","dark green"))
				legend(0.5,-0.1,c("Sweden","Non-Sweden"),pch=c(1,19))
			}
		}
	

	}# for ii 1:3 bac fun oo




}# i 1:4 soil, rs, rp ro



