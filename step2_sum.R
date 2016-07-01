setwd("~/Desktop/drosophila mAC")

df = read.table("dm6_m4C_gene.txt", sep = "\t", header = T)

# READ  reference of CBI from drosophila 
cbi = read.table("dmel_CDS.fasta.indices", sep = "\t", header=T)
library(splitstackshape) #for split strings in data frame
cbi <- cSplit(cbi, "title", sep="_", direction="wide", fixed=T)
cbi <- as.data.frame(cbi)
cbi <- cbi[,c("title_0001", "title_0002", "CBI")]
names(cbi) = c("name", "ID", "CBI")

mrna = read.table("dmel6_gene_stat.txt", sep = "\t", header = T)

ncrna = read.table("dmel6_ncRNA_stat.txt", sep = "\t", header = F)

dfa = merge(df, cbi, by.x = "Gene", by.y="ID")

dfa = merge(dfa, mrna, by.x = "Gene", by.y="transcript_id")

dfa$dens_UTR5 = dfa$UTR5.x/dfa$UTR5.y 
dfa$dens_UTR3 = dfa$UTR3.x/dfa$UTR3.y
dfa$dens_CDS = dfa$CDS.x/dfa$CDS.y  
dfa$dens_Intron = dfa$Intron.x/dfa$Intron.y
dfa$all = dfa$UTR5.x+dfa$CDS.x+dfa$UTR3.x+dfa$Intron.x 
dfa$dens_all = dfa$all/(dfa$UTR5.y+dfa$CDS.y+dfa$UTR3.y+dfa$Intron.y)

#choose the transcript id within the group (gene) with highest density
library(dplyr)
df1 = dfa %>% group_by(gene) %>% top_n(1, dens_all)
df1 = df1[!duplicated(df1$gene),]

par(pch=19, cex=0.25)

par(mfrow=c(2,2))

for (i in c("dens_UTR5", "dens_UTR3", "dens_CDS", "dens_Intron")) {
  cor = cor.test(log2(df1[[i]][df1[[i]]>0]), df1$CBI[df1[[i]]>0])
  plot(log2(df1[[i]][df1[[i]]>0])~df1$CBI[df1[[i]]>0], cex=0.25, 
       col=1, xlab="CBI", ylab= i, main = paste("Person's r=", round(cor$estimate, 2))) }

