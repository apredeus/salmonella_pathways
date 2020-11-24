# Salmonella pathways

When analyzing gene expression and regulation, custom curated pathways are often used to analyze the observed expression changes using Fisher's exact test or GSEA. However, many organisms lack such pathways, leaving only generic sources like KEGG or GO. This collection of curated Salmonella gene sets was obtained from published RNA-seq experiments of wild-type and mutant *Salmonella* strains. 

Final collection of pathways is given using SL1344 gene IDs. These could be easily converted between strains using ortholog tables given in **orthologs** subdirectory. 

## Sources

The data was extracted from the following published sources: 

* [Colgan et al. 2016](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006258) - used up- and down-regulated genes in regulator mutants were taken from [Table S3](https://doi.org/10.1371/journal.pgen.1006258.s011). The expression values are also available from [SalComRegulon](http://bioinf.gen.tcd.ie/cgi-bin/salcom.pl?db=SalComRegulon_HL). 
* [Canals et al. 2019](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000059) - used expression values for 17 *in vitro* conditions for 4/74 and D23580 strains can be taken from [Table S5](https://doi.org/10.1371/journal.pbio.3000059.s005), or from [SalComD23580](http://bioinf.gen.tcd.ie/cgi-bin/salcom_v2.pl?_HL). 
* [Fink et al. 2007](https://jb.asm.org/content/189/6/2262) - supplementary [Table S1](https://jb.asm.org/highwire/filestream/294368/field_highwire_adjunct_files/0/R_Table_S1__Diff_Expressed_genes_in_FNR_.xls).
* [Troxell et al. 2011](https://bmcmicrobiol.biomedcentral.com/articles/10.1186/1471-2180-11-236) - supplementary [Table S2](https://static-content.springer.com/esm/art%3A10.1186%2F1471-2180-11-236/MediaObjects/12866_2011_1522_MOESM2_ESM.XLS).
* [Smith et al. 2016](https://mbio.asm.org/content/7/5/e01024-16) - supplementary [Table S1](https://mbio.asm.org/content/7/5/e01024-16#DC1). 
* SPI1 plus associated effectors, SPI2 plus associated effectors, iron transporters and sideropores - custom pathways compiled from multiple literature sources. 

## Processing 

Custom curated gene sets (SPI1 plus effectors, SPI2 plus effectors, Iron transporters and siderophores) were compiled using many literature sources that are listed in *Salmonella custom pathways.xlsx*. For data from Colgan 2016, Fink 2007, Troxell 2011, and Smith 2016, supplementary tables with differentially expressed genes were used. Whenever the authors used a different reference strain, gene IDs were converted using provided ortholog table (see **orthologs** directory). 

For single-replicate RNA-seq experiments profiling 4/74 and D23580 under 17 *in vitro* conditions, the following strategy was used. Tables of TPM-normalized expression values (see **supp_tables**/474.tsv,D23.tsv) were used to calculate fold change of a particular condition with regard to 16 others. Genes with log2FC of 2 or more were selected as "marker" genes for a particular *in vitro* condition. All lists can be generated using the following command: 

`./calculate_fc_tpm.pl 474.tsv` 

After all gene lists were generated, the following command was used to make one custom pathway table: 

```bash
for i in *.list; do TAG=${i%%.list}; awk -v v=$TAG '{print v"\t"$1}' $i; done > custom_80_pathways.tsv
```

## Case study: analysis of *Salmonella* genes enriched in vacuole and cytosol

An example of pathway use can be found in *enrichment.R* script. Using our pathway compendium, we have annotated *Salmonella* genes up-regulated either in vacuole or cytosol compartments upon HeLa cell infection.
