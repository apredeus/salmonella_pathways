library(KEGGREST)
library(data.table)

library(fgsea)
library(ggplot2)
library(ggrepel)
library(reshape2)

library(clusterProfiler)
library(enrichplot)

setwd("~/salmonella_pathways/vacuole_cytosol")

############## <- do this once - 

path_list <- keggLink("pathway","sey")
pp        <- as.data.table(cbind(path_list,names(path_list)))
write.table(pp,"sey.gene_to_pathway.tsv",quote = F,sep = "\t",col.names = F,row.names = F)
pvec     <- unique(pp$path_list)
NP       <- length(unique(pp$path_list))

for (i in 1:NP) { 
  p_id       <- pvec[i]
  cat(sprintf("Processing pathway number %d, pathway id %s\n",i,p_id))
  p_rec      <- keggGet(p_id)
  p_name     <- gsub(" - Salmonella enterica subsp. enterica serovar Typhimurium SL1344","",p_rec[[1]]$NAME)
  
  p_string   <- paste(p_id,p_name,sep="\t")
  cat(sprintf("%s\n",p_string))
  write(p_string,"sey.kegg_path_name.tsv",append=T)
} 
################# <---------- 

## after this add 80 names to sey pathways 
## we should now have everything we need 

path2SL    <- read.table("sey.gene_to_pathway2.tsv",header = F)
go2SL      <- read.table("sey_go_terms.tsv",header = F)
go_names   <- read.table("GO_names.tsv",sep = '\t',quote = "")
path_names <- read.table("sey.kegg_path_name.tsv", sep = '\t', quote = "")
SL_uni     <- as.vector(read.table("sey.list",header = F)$V1)
## little bit of this so we can annotate things yeah?
govec <- go_names$V2
names(govec) <- go_names$V1
pathvec <- path_names$V2
names(pathvec) <- path_names$V1


Vacuole <- read.table("Vacuole.tsv",header = T,sep = "\t")
Vac_list <- as.vector(Vacuole[Vacuole$SL1344_fixed != "NONE",]$SL1344_fixed)
#write.table(vac_list,"vac_stm.list",quote = F, row.names = F, col.names = F)

Cytosol <- read.table("Cytosol.tsv",header = T,sep = "\t")
Cyt_list <- as.vector(Cytosol[Cytosol$SL1344_fixed != "NONE",]$SL1344_fixed)
#write.table(cyt_list,"cyt_stm.list",quote = F, row.names = F, col.names = F)

en_vac_path <- enricher(Vac_list, universe = SL_uni, minGSSize = 1, qvalueCutoff  = 0.1, TERM2GENE = path2SL)
en_vac_go   <- enricher(Vac_list, universe = SL_uni, minGSSize = 1, qvalueCutoff  = 0.1, TERM2GENE = go2SL)

en_cyt_path <- enricher(Cyt_list, universe = SL_uni, minGSSize = 1, qvalueCutoff  = 0.1, TERM2GENE = path2SL)
en_cyt_go   <- enricher(Cyt_list, universe = SL_uni, minGSSize = 1, qvalueCutoff  = 0.1, TERM2GENE = go2SL)

## annotate with pathway/GO names
for (i in 1:length(en_vac_go@result$ID)) {
  go_id <- en_vac_go@result$ID[i]
  en_vac_go@result$Description[i] <- as.character(govec[go_id])
  en_vac_go@result$GR[i] <- as.numeric(unlist(strsplit(en_vac_go@result$GeneRatio[i],"/"))[1])/as.numeric(unlist(strsplit(en_vac_go@result$GeneRatio[i],"/"))[2])
}

for (i in 1:length(en_cyt_go@result$ID)) {
  go_id <- en_cyt_go@result$ID[i]
  en_cyt_go@result$Description[i] <- as.character(govec[go_id])
  en_cyt_go@result$GR[i] <- as.numeric(unlist(strsplit(en_cyt_go@result$GeneRatio[i],"/"))[1])/as.numeric(unlist(strsplit(en_cyt_go@result$GeneRatio[i],"/"))[2])
}

for (i in 1:length(en_vac_path@result$ID)) {
  path_id <- en_vac_path@result$ID[i]
  en_vac_path@result$Description[i] <- as.character(pathvec[path_id])
  en_vac_path@result$GR[i] <- as.numeric(unlist(strsplit(en_vac_path@result$GeneRatio[i],"/"))[1])/as.numeric(unlist(strsplit(en_vac_path@result$GeneRatio[i],"/"))[2])
}

for (i in 1:length(en_cyt_path@result$ID)) {
  path_id <- en_cyt_path@result$ID[i]
  en_cyt_path@result$Description[i] <- as.character(pathvec[path_id])
  en_cyt_path@result$GR[i] <- as.numeric(unlist(strsplit(en_cyt_path@result$GeneRatio[i],"/"))[1])/as.numeric(unlist(strsplit(en_cyt_path@result$GeneRatio[i],"/"))[2])
}

## look at the results formatted as a table
DT::datatable(as.data.frame(en_vac_go)) 
DT::datatable(as.data.frame(en_vac_path))
DT::datatable(as.data.frame(en_cyt_go))
DT::datatable(as.data.frame(en_cyt_path))
