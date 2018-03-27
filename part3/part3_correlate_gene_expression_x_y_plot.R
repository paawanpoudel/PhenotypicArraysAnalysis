

rm(list=ls())
options(stringsAsFactors = FALSE)

# sourcing the files
r_fol <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS//Misc/multiorgan_crc/R/"
r_files <- list.files(path = r_fol, pattern = ".R", full.names = TRUE)
lapply(r_files, source)



r_fol <- "/Users/ppoudel/Dropbox (SPM)/Pawan/Scripts/MULTIORGAN_RCODES/2017/functions/source/"
r_files <- list.files(path = r_fol, pattern = ".R", full.names = TRUE)

lapply(r_files, source)

require(ggplot2)
require(reshape2)
require(plyr)
# reading the metabolic genes

output_dir <-  "/Users/ppoudel/Dropbox (SPM)//Pawan/ANALYSIS/Internal/Pancreas/Biolog/analysis/opm_output/carbon_plates/opm_output/PM01/gene_metabolite_correlation/"
#metabolic_genes_file <- "/Users/ppoudel/Dropbox (SPM)//Pawan/ANALYSIS/Internal/Pancreas/Biolog//supporting_materials/metabolic_genes_supplementary table 2 nature10350-s2_v2.txt"
#metabolic_genes <- read.table(metabolic_genes_file, sep = "\t", header = TRUE)

# reading the cell line data
cell_line_file <- "/Users/ppoudel/Dropbox (ICR)/Pawan/DATA/External/Pancreas/Cell lines/KlijnEtAl2014/Rdata/20170509_PP_KlijnEtAl2014_coding_genes_tmm_cpm_log2_plus1_filtered_annotated.Rdata"
load(cell_line_file)

colnames(annotated_pda_tmm_cpm_log2_plus1) <- toupper( colnames(annotated_pda_tmm_cpm_log2_plus1))
colnames(annotated_pda_tmm_cpm_log2_plus1) <- gsub( "PATU", "", colnames(annotated_pda_tmm_cpm_log2_plus1) )
colnames(annotated_pda_tmm_cpm_log2_plus1) <- gsub( "PANC0327", "PANC327", colnames(annotated_pda_tmm_cpm_log2_plus1) )


# now correlating the metabolite with the gene expression data and selecting the gene have the highest and the lowest correlation

# now selecting the metabolic genes
load("/Users/ppoudel/Dropbox (SPM)//Pawan/ANALYSIS/Internal/Pancreas/Biolog/rdata//PM01_Interesting_Metabolites.RData")

costas_metabolic_genes <- c("SLC44A1", "UPP1", "ACBD5", "MAN2A2", "TPI", "GAPDH", "OGFOD2", "CYBRD1", "CHSY1", "SLC25A12", "PCYT2")

plate2analyse<- "PM01"

lapply(costas_metabolic_genes, function(x, gene_expression_data=annotated_pda_tmm_cpm_log2_plus1, metabolite_data=PM.data.normalise.significant.metabolites, plate_num= plate2analyse, od=output_dir){

  setwd(od)
  my_genes <- paste0(x, " gene expression")
  
  #gene_exp_1 <-  sweep( gene_expression_data[ gene_expression_data[,1] %in% x, -c(1)], 1, rowMeans(gene_expression_data[ gene_expression_data[,1] %in% x, -c(1)]), "-" )
  #data2plot_genes <-  t( data.frame(gene_expression_data[ gene_expression_data[,1] %in% x, c(1)]  , gene_exp_1))
  
  data2plot_genes <- t( gene_expression_data[ gene_expression_data[,1] %in% x,  ] )
  
  if( ncol(data2plot_genes) >=1){
    
    data2plot_genes_metabolites <- merge( data2plot_genes, PM.data.normalise.significant.metabolites, by.x="row.names", by.y="row.names")
    
    data2plot <- melt(data2plot_genes_metabolites)
    colnames(data2plot) <- c( "cell_lines", "gene_expression", "metabolite_names", "metabolite_abundance")
    
    data2plot$gene_expression <- as.numeric( data2plot$gene_expression) 
    data2plot$metabolite_abundance <- as.numeric( data2plot$metabolite_abundance) 
    
  
    # Calculate correlation for each group
    
    cors <- ddply(data2plot, .(metabolite_names ), summarise, cor = round(cor(gene_expression, metabolite_abundance, method = "spearman"), 4))
    x_val <- ddply(data2plot, .(metabolite_names ), summarise, x_p = min(gene_expression) +1)
    y_val <- ddply(data2plot, .(metabolite_names ), summarise, y_p = max(metabolite_abundance) -5)
    
    cors_x <- merge(cors, x_val, by.x = colnames(cors)[1], by.y = colnames(x_val)[1] )
    cors_xy <- merge(cors_x, y_val, by.x = colnames(cors_x)[1], by.y = colnames(y_val)[1] )
    
    
    min_x <- min(data2plot$gene_expression) 
    min_y <- mean(data2plot$metabolite_abundance)
    
    p <- ggplot(data2plot, aes(x=gene_expression, y=metabolite_abundance )) +
      geom_point( aes( colour=factor(cell_lines)))+ facet_wrap(~metabolite_names, scales = "free") 
    
    p <- p + geom_text(data=cors_xy, aes(label=paste("r=", cor, sep="")), x=cors_xy$x_p, y=cors_xy$y_p )+ theme_Publication() + labs(x=my_genes,y= "Normalised metabolic uptake") 
    
    plot_file <- paste0( Sys.Date(), "_", plate_num,"_",my_genes, "_spearman_correlation_analysis.pdf")
    plot_file_txt <- paste0( Sys.Date(), "_", plate_num,"_",my_genes, "_spearman_correlation_analysis.txt")
    write.table( data2plot_genes_metabolites, plot_file_txt, sep = "\t", quote = FALSE)
    
    ggsave(p, filename = plot_file)
    
  }
  
  
})






