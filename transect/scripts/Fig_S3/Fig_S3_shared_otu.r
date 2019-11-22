###
### need to be called separate for bacteria,fungi and oomycetes 
###

###
### you also need to switch between high and low abundant OTUs
###


#### for Panel b, choose one

tab=read.table("bac_high_abun.txt",header=F)
#tab=read.table("fungi_high_abun.txt",header=F)
#tab=read.table("oomyc_high_abun.txt",header=F)

#### for Panel c, choose one

#tab=read.table("bac_low_abun.txt",header=F)
#tab=read.table("fungi_low_abun.txt",header=F)
#tab=read.table("oomyc_low_abun.txt",header=F)



#####


new_nam=paste(tab[,9],tab[,10])

all_root=new_nam=="RootRoot within"
all_soil=new_nam=="SoilSoil within"
all_rs=new_nam=="RSRS within"
all_rp=new_nam=="RPRP within"

soil1=which(tab[,9]=="SoilSoil" & tab[,10] == "within_.Y1")
soil2=which(tab[,9]=="SoilSoil" & tab[,10] == "within_.Y2")
soil3=which(tab[,9]=="SoilSoil" & tab[,10] == "within_.Y3")
soil4=c(soil1,soil2,soil3)

root1=which(tab[,9]=="RootRoot" & tab[,10] == "within_.Y1")
root2=which(tab[,9]=="RootRoot" & tab[,10] == "within_.Y2")
root3=which(tab[,9]=="RootRoot" & tab[,10] == "within_.Y3")
root4=c(root1,root2,root3)

dev.new(height=5.044,width=8.093)
par(mfrow=c(1,4))

boxplot(rbind(as.matrix(tab[soil1,6]),as.matrix(tab[soil1,7])),xlim=c(0.5,2.5),ylim=c(0,1),las=2,at=1,col="brown",ylab="% shared OTUs per sample",main="Year1 ",width=1.5)
boxplot(rbind(as.matrix(tab[root1,6]),as.matrix(tab[root1,7])),add=T,at=2,col="dark green",yaxt="n",width=1.5)
t_res1=t.test(rbind(as.matrix(tab[soil1,6]),as.matrix(tab[soil1,7])),rbind(as.matrix(tab[root1,6]),as.matrix(tab[root1,7])))
text(1,1,t_res1$p.value,pos=4,cex=0.7,col="red")

boxplot(rbind(as.matrix(tab[soil2,6]),as.matrix(tab[soil2,7])),xlim=c(0.5,2.5),ylim=c(0,1),las=2,at=1,col="brown",ylab="% shared OTUs per sample",main="Year 2",width=1.5)
boxplot(rbind(as.matrix(tab[root2,6]),as.matrix(tab[root2,7])),add=T,at=2,col="dark green",yaxt="n",width=1.5)
t_res2=t.test(rbind(as.matrix(tab[soil2,6]),as.matrix(tab[soil2,7])),rbind(as.matrix(tab[root2,6]),as.matrix(tab[root2,7])))
text(1,1,t_res2$p.value,pos=4,cex=0.7,col="red")


boxplot(rbind(as.matrix(tab[soil3,6]),as.matrix(tab[soil3,7])),xlim=c(0.5,2.5),ylim=c(0,1),las=2,at=1,col="brown",ylab="% shared OTUs per sample",main="Year 3",width=1.5)
boxplot(rbind(as.matrix(tab[root3,6]),as.matrix(tab[root3,7])),add=T,at=2,col="dark green",yaxt="n",width=1.5)
t_res3=t.test(rbind(as.matrix(tab[soil3,6]),as.matrix(tab[soil3,7])),rbind(as.matrix(tab[root3,6]),as.matrix(tab[root3,7])))
text(1,1,t_res3$p.value,pos=4,cex=0.7,col="red")

boxplot(rbind(as.matrix(tab[soil4,6]),as.matrix(tab[soil4,7])),xlim=c(0.5,2.5),ylim=c(0,1),las=2,at=1,col="brown",ylab="% shared OTUs per sample",main="All years",width=1.5)
boxplot(rbind(as.matrix(tab[root4,6]),as.matrix(tab[root4,7])),add=T,at=2,col="dark green",yaxt="n",width=1.5)
t_res4=t.test(rbind(as.matrix(tab[soil4,6]),as.matrix(tab[soil4,7])),rbind(as.matrix(tab[root4,6]),as.matrix(tab[root4,7])))
text(1,1,t_res4$p.value,pos=4,cex=0.7,col="red")




