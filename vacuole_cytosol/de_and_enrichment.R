library(data.table)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(DESeq2)

setwd("~/salmonella_pathways/vacuole_cytosol")



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
res_label <- res_label[res_label$significant == "FDR < 0.05",] ## only keep the significant ones
res_label[res_label$Gene_name == "STM1868A",]$Gene_name <- "SL1802"
res_label[res_label$Gene_name == "STM4042",]$Gene_name  <- "SL3990"

res_label2 <- res3[res3$logp > 18,]
res_label2[res_label2$Gene_name == "STM2243",]$Gene_name <- "SL2219"
res_label3 <- res3[res3$logp > 7 & res3$log2FC > 0,]

all_labels <- rbind(res_label,res_label2,res_label3)

res3$confirmed <- "not confirmed"
res3[res3$`474_ID` %in% res_label$`474_ID`,]$confirmed <- "confirmed"

res4 <- res3[res3$confirmed == "not confirmed",]
res5 <- res3[res3$confirmed == "confirmed",]

## this is the best I came up with (and it's honestly not bad)

theme_set(theme_bw(base_size = 14))
ggplot(res3, aes(x = log2FC, y = logp)) +
  geom_point(aes(color = regulation, size = log10(mean_tpm+1)),stroke = 0) + 
  geom_point(aes(color = confirmed, size = log10(mean_tpm+1)),stroke = 0) + 
  geom_label_repel(data = all_labels, aes(label = Gene_name), size = 4, max.overlaps = 20) + 
  scale_size_continuous(range = c(0.1,4)) + 
  scale_color_manual(values = c("not significant"="#999999","up (cytosol)"="#FF0000","up (vacuole)"="#0000FF","not confirmed"="#FFFFFF00","confirmed"="#000000")) + 
  xlab("log2 fold change") + 
  ylab("-log10 (adjusted p-value)") + 
  ggtitle("Volcano plot, EBSS vs. WTM (vacuole vs. cytosol)") + 
  theme(panel.grid.major = element_line(size = 0.25))



################################# PART 2 - Pathway enrichment analysis ####################################

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

path_genes    <- read.table("Custom.pathways.tsv",header = F)
go_genes     <- read.table("GO.pathways.tsv",header = F)
kegg_genes    <- read.table("KEGG.pathways.tsv",header = F)

path_names <- read.table("Custom.names.tsv", sep = '\t', quote = "")
go_names   <- read.table("GO.names.tsv",sep = '\t',quote = "")
kegg_names <- read.table("KEGG.names.tsv",sep = '\t',quote = "")

SL_universe     <- read.table("474_ID.list",header = F)$V1

all_pathways <- rbind(path_genes,go_genes,kegg_genes)
all_names    <- rbind(path_names,go_names,kegg_names)
namevec <- all_names$V2
names(namevec) <- all_names$V1

#### now get lists 
vacuole <- read.table("vacuole.tsv",header = T,sep = "\t",check.names = F)
vac_genes <- vacuole$`474_ID`
cytosol <- read.table("cytosol.tsv",header = T,sep = "\t",check.names = F)
cyt_genes <- cytosol$`474_ID`

### calculate Fisher enrichment score
vac_pathways <- enricher(vac_genes, universe=SL_universe, maxGSSize=1000, minGSSize=3, qvalueCutoff=0.1, TERM2GENE=all_pathways)
cyt_pathways <- enricher(cyt_genes, universe=SL_universe, maxGSSize=1000, minGSSize=3, qvalueCutoff=0.1, TERM2GENE=all_pathways)

## annotate with pathway/GO names
for (i in 1:length(vac_pathways@result$ID)) {
  pathway_id <- vac_pathways@result$ID[i]
  vac_pathways@result$Description[i] <- as.character(namevec[pathway_id])
  vac_pathways@result$GR[i] <- 
    as.numeric(unlist(strsplit(vac_pathways@result$GeneRatio[i],"/"))[1])/as.numeric(unlist(strsplit(vac_pathways@result$GeneRatio[i],"/"))[2])
}

for (i in 1:length(cyt_pathways@result$ID)) {
  pathway_id <- cyt_pathways@result$ID[i]
  cyt_pathways@result$Description[i] <- as.character(namevec[pathway_id])
  cyt_pathways@result$GR[i] <- 
    as.numeric(unlist(strsplit(cyt_pathways@result$GeneRatio[i],"/"))[1])/as.numeric(unlist(strsplit(cyt_pathways@result$GeneRatio[i],"/"))[2])
}

## take a quick look at things
DT::datatable(as.data.frame(vac_pathways))
DT::datatable(as.data.frame(cyt_pathways))

vac_df <- as.data.frame(vac_pathways)
cyt_df <- as.data.frame(cyt_pathways)

vac_df$type <- "vacuole"
cyt_df$type <- "cytosol"
rownames(vac_df) <- 1:nrow(vac_df)
rownames(cyt_df) <- 1:nrow(cyt_df)

## make dataframe of all enriched pathways
all_res      <- rbind(vac_df,cyt_df)
all_res$type <- as.factor(all_res$type)
all_res$Description <- factor(all_res$Description)

## annotate the dataframe - convert overlapping gene IDs to gene names for interpretation
gene_names <- rbind(vacuole,cytosol)[,c(1,4)]
gnamevec <- gene_names$Gene_name
names(gnamevec) <- gene_names$`474_ID`

all_res$Names <- "NA"

for (j in 1:nrow(all_res)) { 
  ids <- all_res[j,8]
  ids_list <- as.list(unlist(strsplit(ids,"/")))
  namestr <- ""
  for (k in 1:length(ids_list)) { 
    id <- ids_list[[k]]
    namestr <- ifelse (namestr == "",gnamevec[[id]],paste(namestr,gnamevec[[id]],sep = ","))
  } 
  all_res[j,12] <- namestr
}

write.table(all_res,"annotated_significant_pathways.tsv",quote = F,sep = "\t",row.names = F)

## we then inspect the pathways and decide which ones to keep and which ones to remove. 

curated           <- read.table("curated_pathways.tsv",header = T,sep = "\t")
curated           <- curated[curated$Retain == "yes",]
curated$Names     <- NULL
curated$Description <- NULL

curated$Condition <- factor(curated$Condition,levels=c("vacuole","cytosol"))
curated$Type      <- factor(curated$Type,levels = c("regulon","other"))
curated$Simple    <- factor(curated$Simple)

## plotting works best like this
cbPalette     <- c("#0000ff", "#ff0000")


preg <- curated[curated$Type == "regulon",]
poth <- curated[curated$Type == "other",]

preg$Simple <- factor(preg$Simple,levels = c("DinvF (Smith) - DN","DssrB - UP","DrtsA (Smith) - DN","DhilA - DN",
                                                 "DhilA (Smith) - DN","DhilD (Smith) - DN","DhilD (MEP) - DN","DhilD (ESP) - DN","DrpoS - UP",
                                                 "DphoPQ - UP","DrpoS - DN","Dfur - DN","DslyA - UP",
                                                 "DompR DenvZ - DN","DfliZ - DN","DbarA DsirA - DN",
                                                 "Dhfq - DN","Dfur - UP","DslyA - DN","DssrB - DN","DphoPQ - DN"))
poth$Simple <- factor(poth$Simple)

p1 <- ggplot(preg,aes(x = Simple,y = Count,fill = Condition)) + geom_bar(stat = "identity") + coord_flip() + 
  scale_fill_manual(values = cbPalette) + ylab("Number of genes (log2 scale)") + xlab("") + facet_grid(cols=vars(Condition)) + theme_bw() + 
  scale_y_continuous(trans = 'log2',breaks = c(0,5,10,20,40,80,160),limits = c(1,200)) + theme(text = element_text(size = 14))
p2 <- ggplot(poth,aes(x = reorder(Simple,Count,FUN = sum),y = Count,fill = Condition)) + geom_bar(stat = "identity") + coord_flip() + 
  scale_fill_manual(values = cbPalette) + ylab("Number of genes (log2 scale)") + xlab("") + facet_grid(cols=vars(Condition)) + theme_bw() + 
  scale_y_continuous(trans = 'log2',breaks = c(0,5,10,20,40,80,160),limits = c(1,200)) + theme(text = element_text(size = 14))
p1 / p2 + plot_layout(heights = c(1.6,1))







