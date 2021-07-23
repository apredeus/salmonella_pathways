library(data.table)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(DESeq2)

setwd("~/salmonella_pathways/vacuole_cytosol")

## reference file structure and preparation
## there are three groups of pathways: Custom, GO, and KEGG
## each group has a pathways.tsv file (tab-separated 2 column list, pathway - gene ID), 
## and names.tsv - name of each pathway (e.g. "oxidative phosphorylation")
## the process of custom pathway curation is described in the README and the paper; 
## GO names are taken from GO website, and GO mapping (pathway to gene) is obtained with InterProScan;
## KEGG mapping is obtained using KEGGREST package as follows: 

## path_list <- keggLink("pathway","sey")
## pp        <- as.data.table(cbind(path_list,names(path_list)))
## write.table(pp,"KEGG.pathways.tsv",quote = F,sep = "\t",col.names = F,row.names = F)
## pvec     <- unique(pp$path_list)
## NP       <- length(unique(pp$path_list))

## for (i in 1:NP) { 
##   p_id       <- pvec[i]
##  cat(sprintf("Processing pathway number %d, pathway id %s\n",i,p_id))
##  p_rec      <- keggGet(p_id)
##  p_name     <- gsub(" - Salmonella enterica subsp. enterica serovar Typhimurium SL1344","",p_rec[[1]]$NAME)
##   p_string   <- paste(p_id,p_name,sep="\t")
##  cat(sprintf("%s\n",p_string))
##  write(p_string,"KEGG.names.tsv",append=T)
## } 

## strangely enough, you also need to run dos2unix on the KEGGREST output.. 

################################### PART 1 - Differential expression analysis ######################################

tt            <- read.table("474.annotated.counts.tsv", header=T, check.names=F)
head(tt)
rownames(tt)  <- tt$Locus_tag
tt$Locus_tag  <- NULL
cond          <- read.table("Conditions.txt", header=T, row.names=1, check.names=F,stringsAsFactors = T)

exp           <- tt[,-c(1:3)]

dim(tt)
dim(exp)
barplot(colSums(exp))  ## 9.5-17.5M useful reads per replicate

dds           <- DESeqDataSetFromMatrix(countData = round(exp,0),colData = cond,design = ~ Condition)
deseq         <- DESeq(dds)

## let's also read in TPMs for annotation purposes
tt_tpm           <- read.table("474.annotated.TPM.tsv", header=T, check.names=F)
rownames(tt_tpm) <- tt_tpm$Locus_tag
tt_tpm$Locus_tag <- NULL
tpm              <- tt_tpm[,-c(1:3)]

## our output will have log2FC of EBSS/WTM = cytosole/vacuole

res_orig      <- results(deseq,contrast=c("Condition","WTM","EBSS"))
summary(res_orig,alpha = 0.05)  ## 216 genes up in cytosole, 443 up in vacuole 

ann           <- tt[,c(1:3)]
ann2          <- read.table("474_to_SL.tsv",sep = "\t",header = F,row.names = 1)
ann3          <- read.table("474_to_STM.tsv",sep = "\t",header = F,row.names = 1)

ann           <- merge(ann,ann2,by = 0,all.x = T)
rownames(ann) <- ann$Row.names
ann$Row.names <- NULL
ann           <- merge(ann,ann3,by = 0,all.x = T)
rownames(ann) <- ann$Row.names
ann$Row.names <- NULL
colnames(ann)[4:5] <- c("SL1344_ID","LT2_ID")
ann$SL1344_ID[is.na(ann$SL1344_ID)] <- "NONE"
ann$LT2_ID[is.na(ann$LT2_ID)] <- "NONE"

### this is our annotation, with SL and STM tags 
head(ann)
dim(ann)

## remove unnecessary columns
res           <- as.data.frame(res_orig)
res$lfcSE     <- NULL
res$stat      <- NULL
res$pvalue    <- NULL
colnames(res) <- c("baseMean","log2FC","padj")
res$baseMean  <- round(res$baseMean,digits=0)
res$log2FC    <- round(res$log2FC,digits=2)
res$padj      <- formatC(res$padj, format="e", digits=2)

## also add mean TPMs per condition
tpm1          <- tpm[,grepl("EBSS",colnames(tpm))]
tpm2          <- tpm[,grepl("WTM",colnames(tpm))]
tpm1$tpm1_mean  <- apply(tpm1,1,mean)
tpm2$tpm2_mean  <- apply(tpm2,1,mean)
tpm1$tpm1_mean  <- round(tpm1$tpm1,digits=3)
tpm2$tpm2_mean  <- round(tpm2$tpm2,digits=3)
tpms           <- merge(tpm1,tpm2,by="row.names")
rownames(tpms) <- tpms$Row.names
tpms$Row.names <- NULL

## annotate the output with all these things!
res            <- merge(res,tpms,by="row.names")
row.names(res) <- res$Row.names
res$Row.names  <- NULL
res            <- merge(res,ann,by="row.names")
colnames(res)[1] <- "474_ID"

## reorder, and save filtered and un-filtered tables
res2           <- res[,c(1,14,15,11:13,5:10,2:4)]
res3           <- res2[complete.cases(res2) & res2$padj != " NA",]
res3$padj      <- as.numeric(res3$padj)
res4           <- res3[res3$padj <= 0.05,]

write.table(res2[order(res2$log2FC,decreasing=T),],"EBSS_vs_WTM.deseq2_full.tsv",sep = "\t",quote=F,row.names=F)
write.table(res4[order(res4$log2FC,decreasing=T),],"EBSS_vs_WTM.deseq2_filt.tsv",sep = "\t",quote=F,row.names=F)



### let's make a fancy volcano plot 

res3$logp      <- -log10(res3$padj)
res3$mean_tpm  <- (res3$tpm1_mean + res3$tpm2_mean)/2

res3$significant <- "Not sig"
res3$significant[res3$padj <= 0.05] <- "FDR < 0.05"
## make Volcano plot with dot size proportional to expression

res3$regulation <- "not significant"
res3$regulation[res3$log2FC > 0 & res3$padj <= 0.05] <- "up (cytosol)"
res3$regulation[res3$log2FC < 0 & res3$padj <= 0.05] <- "up (vacuole)"

confirmed <- read.table("confirmed.list")$V1
res_label <- res3[res3$Gene_name %in% confirmed,]

res3$regulation[res3$Gene_name %in% res_label$Gene_name & res3$log2FC > 0] <- "confirmed (cytosol)"
res3$regulation[res3$Gene_name %in% res_label$Gene_name & res3$log2FC < 0] <- "confirmed (vacuole)"

ggplot(res_label, aes(x = log2FC, y = logp)) +
  geom_point(data = res_label, aes(size = log10(mean_tpm+2)))

theme_set(theme_bw(base_size = 14))
ggplot(res3, aes(x = log2FC, y = logp)) +
  geom_point(aes(color = regulation, size = log10(mean_tpm+1)), stroke = 0) + 
  geom_label_repel(data = res_label, aes(label = Gene_name), size = 4) + 
  scale_size_continuous(range = c(0.1,4)) + 
  scale_color_manual(values = c("not significant"="#999999","up (cytosol)"="#E86850","up (vacuole)"="#587498","confirmed (cytosol)"="#00FF00","confirmed (vacuole)"="#00FFFF")) + 
  xlab("log2 fold change") + 
  ylab("-log10 (adjusted p-value)") + 
  ggtitle("Volcano plot, EBSS vs. WTM (vacuole vs. cytosol)") + 
  theme(panel.grid.major = element_line(size = 0.25))


################################# PART 2 - Pathway enrichment analysis ####################################


## write_kegg_gmt("smar")

## write_kegg_gmt is being finicky - seems to fail at keggGet("md:stm_M00001") requests
## so I will use the old GMT file generated for Enteritidis

## it seems like there are more genes (~ 1750 vs 1500) available in gene symbols vs gene names. 
## so as far as KEGG mappings go, it's good to use GMT that use original locus tags.

## GO to names obtained by parsing GO .obo file downloaded from geneontology.org

## get a list of KEGG pathways and modules
## then convert all the IDs to gene names using unique locus_tag to name table from GFF
## and write a combined GMT file with all the names 

############## <- do this once 

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

path2SL    <- read.table("sey.gene_to_pathway3.tsv",header = F)
go2SL      <- read.table("sey_go_terms.tsv",header = F)
go_names   <- read.table("GO_names.tsv",sep = '\t',quote = "")
path_names <- read.table("sey.kegg_path_name.tsv", sep = '\t', quote = "")
SL_uni     <- as.vector(read.table("sey.list",header = F)$V1)
## little bit of this so we can annotate things yeah?
govec <- go_names$V2
names(govec) <- go_names$V1
pathvec <- path_names$V2
names(pathvec) <- path_names$V1


Vacuole <- read.table("vacuole.tsv",header = T,sep = "\t")
Vac_list <- as.vector(Vacuole[Vacuole$SL1344_fixed != "NONE",]$SL1344_fixed)
#write.table(vac_list,"vac_stm.list",quote = F, row.names = F, col.names = F)

Cytosol <- read.table("cytosol.tsv",header = T,sep = "\t")
Cyt_list <- as.vector(Cytosol[Cytosol$SL1344_fixed != "NONE",]$SL1344_fixed)
#write.table(cyt_list,"cyt_stm.list",quote = F, row.names = F, col.names = F)

en_vac_path <- enricher(Vac_list, universe = SL_uni, maxGSSize = 1000, minGSSize = 1, qvalueCutoff  = 0.1, TERM2GENE = path2SL)
en_vac_go   <- enricher(Vac_list, universe = SL_uni, maxGSSize = 1000, minGSSize = 1, qvalueCutoff  = 0.1, TERM2GENE = go2SL)

en_cyt_path <- enricher(Cyt_list, universe = SL_uni, maxGSSize = 1000, minGSSize = 1, qvalueCutoff  = 0.1, TERM2GENE = path2SL)
en_cyt_go   <- enricher(Cyt_list, universe = SL_uni, maxGSSize = 1000, minGSSize = 1, qvalueCutoff  = 0.1, TERM2GENE = go2SL)

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


DT::datatable(as.data.frame(en_vac_go)) ## I likes!!!
DT::datatable(as.data.frame(en_vac_path))
DT::datatable(as.data.frame(en_cyt_go))
DT::datatable(as.data.frame(en_cyt_path))

tt1 <- rbind(as.data.frame(en_vac_path),as.data.frame(en_vac_go))
tt2 <- rbind(as.data.frame(en_cyt_path),as.data.frame(en_cyt_go))

tt1$type <- "vacuole"
tt2$type <- "cytosol"
rownames(tt1) <- 1:nrow(tt1)
rownames(tt2) <- 1:nrow(tt2)
all_res <- rbind(tt1,tt2)
all_res$type <- as.factor(all_res$type)
write.table(all_res,"all_sig_pathways.tsv",quote = F,row.names = F,sep = "\t")
all_res$Description <- factor(all_res$Description)

## oui shall plot
cbPalette     <- c("#0000ff", "#ff0000")
fn            <- read.table("fix_names2.tsv",header = T,sep = "\t")
ar1           <- merge(all_res,fn,by.x = "ID",by.y = "ID")
ar1$Category  <- ifelse(grepl("UP",ar1$Fix_name),"regulon","other")
ar1$Category  <- ifelse(grepl("DN",ar1$Fix_name),"regulon",ar1$Category)
ar1$Category  <- ifelse(grepl("D23",ar1$Fix_name),"other",ar1$Category)
ar1$Category  <- ifelse(grepl("4/74",ar1$Fix_name),"other",ar1$Category)

ar1$type      <- factor(ar1$type,levels = c("vacuole","cytosol"))
ar1$Category  <- factor(ar1$Category,levels = c("regulon","other"))

ar1$Description <- NULL
ar1$Description2 <- NULL
ar1$pvalue    <- NULL
ar1$qvalue    <- NULL
ar1$geneID    <- NULL
ar1$Name      <- NULL

##### change select.list file if you want more or fewer

sel <- read.table("select.list",header = F,sep = "\t")$V1
ar2 <- ar1[ar1$Fix_name %in% sel,]

ggplot(ar2,aes(x = reorder(Fix_name,Count,FUN = sum),y = Count,fill = type)) + geom_bar(stat = "identity") + coord_flip() + 
  scale_fill_manual(values = cbPalette) + ylab("Number of genes in overlap, log2 scale") + xlab("") + 
  facet_grid(cols=vars(type),rows=vars(Category)) + theme_bw() + 
  scale_y_continuous(trans = 'log2',breaks = c(0,5,10,20,40,80,160))

preg <- ar2[ar2$Category == "regulon",]
poth <- ar2[ar2$Category == "other",]

preg$Fix_name <- factor(preg$Fix_name,levels = c("DinvF (Smith) - DN","DssrB - UP","DrtsA (Smith) - DN","DhilA - DN",
                                                 "DhilA (Smith) - DN","DhilD (Smith) - DN","DhilD (MEP) - DN","DhilD (ESP) - DN","DrpoS - UP","DsprB (Smith) - DN",
                                                 "DphoPQ - UP","DrpoS - DN","Dfur - DN","DslyA - UP",
                                                 "DompR DenvZ - DN","DfliZ - DN","DbarA DsirA - DN",
                                                 "Dhfq - DN","Dfur - UP","DslyA - DN","DssrB - DN","DphoPQ - DN"))
poth$Fix_name <- factor(poth$Fix_name)

p1 <- ggplot(preg,aes(x = Fix_name,y = Count,fill = type)) + geom_bar(stat = "identity") + coord_flip() + 
  scale_fill_manual(values = cbPalette) + ylab("Number of genes (log2 scale)") + xlab("") + facet_grid(cols=vars(type)) + theme_bw() + 
  scale_y_continuous(trans = 'log2',breaks = c(0,5,10,20,40,80,160),limits = c(1,160)) + theme(text = element_text(size = 14))
p2 <- ggplot(poth,aes(x = reorder(Fix_name,Count,FUN = sum),y = Count,fill = type)) + geom_bar(stat = "identity") + coord_flip() + 
  scale_fill_manual(values = cbPalette) + ylab("Number of genes (log2 scale)") + xlab("") + facet_grid(cols=vars(type)) + theme_bw() + 
  scale_y_continuous(trans = 'log2',breaks = c(0,5,10,20,40,80,160),limits = c(1,160)) + theme(text = element_text(size = 14))
p1 / p2 + plot_layout(heights = c(1.6,1))

names <- list()
sl_names <- read.table("sl_names.tsv",header = T)

for (i in 1:nrow(sl_names)) { 
  id <- sl_names[i,1]
  gene_name <- sl_names[i,2]
  names[[id]] <- gene_name
}

all_res$Name <- "NA"
for (j in 1:nrow(all_res)) { 
  ids <- all_res[j,8]
  ids_list <- as.list(unlist(strsplit(ids,"/")))
  namestr <- ""
  for (k in 1:length(ids_list)) { 
    id <- ids_list[[k]]
    namestr <- ifelse (namestr == "",names[[id]],paste(namestr,names[[id]],sep = ","))
  } 
  all_res[j,12] <- namestr
}

write.table(all_res,"ann_sig_pathways4.tsv",quote = F,sep = "\t",row.names = F)
dsel <- fn[fn$Fix_name %in% sel,]$Description2
all_resf <- all_res[all_res$Description %in% dsel,]
write.table(all_resf,"select_sig_pathways4.tsv",quote = F,sep = "\t",row.names = F)

