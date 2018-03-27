rm(list=ls())
options(stringsAsFactors = FALSE)

# sourcing the files
r_fol <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS//Misc/multiorgan_crc/R/"
r_files <- list.files(path = r_fol, pattern = ".R", full.names = TRUE)
lapply(r_files, source)



r_fol <- "/Users/ppoudel/Dropbox (SPM)//Pawan/Scripts/MULTIORGAN_RCODES/2017/functions/source/"
r_files <- list.files(path = r_fol, pattern = ".R", full.names = TRUE)
lapply(r_files, source)

require(ggplot2)
require(reshape2)

# reading the metabolic genes

pm_sig_metabolite_rda <- "/Users/ppoudel/Dropbox (SPM)//Pawan/ANALYSIS/Internal/Pancreas//Biolog/rdata//PM01_Interesting_Metabolites.RData"
plate_num <- "PM01"
output_dir <- "/Users/ppoudel/Dropbox (SPM)//Pawan/ANALYSIS/Internal/Pancreas//Biolog/analysis/opm_output/carbon_plates/opm_output/PM01/gene_metabolite_correlation/"
dir.create(output_dir, showWarnings = FALSE)


# running the script
gene_metabolite_correlation_published_metabolic_genes(plate_num, metabolic_genes, pm_sig_metabolite_rda, output_dir)
gene_metabolite_correlation_all_transcriptome(plate_num, pm_sig_metabolite_rda, output_dir)
    
metabolic_genes_file <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/supporting_materials/metabolic_genes_supplementary table 2 nature10350-s2_v2.txt"
metabolic_genes <- unique( read.table(metabolic_genes_file, sep = "\t", header = TRUE)[,2] )


gene_metabolite_correlation_published_metabolic_genes <- function(plate_num, metabolic_genes, pm_sig_metabolite_rda, output_dir){
  
  load(pm_sig_metabolite_rda)
   
  # reading the cell line data
  cell_line_file <- "/Users/ppoudel/Dropbox (ICR)//Pawan/DATA/External/Pancreas/Cell lines/KlijnEtAl2014/Rdata/20170509_PP_KlijnEtAl2014_coding_genes_tmm_cpm_log2_plus1_filtered_annotated.Rdata"
  load(cell_line_file)
  
  colnames(annotated_pda_tmm_cpm_log2_plus1) <- toupper( colnames(annotated_pda_tmm_cpm_log2_plus1))
  colnames(annotated_pda_tmm_cpm_log2_plus1) <- gsub( "PATU", "", colnames(annotated_pda_tmm_cpm_log2_plus1) )
  colnames(annotated_pda_tmm_cpm_log2_plus1) <- gsub( "PANC0327", "PANC327", colnames(annotated_pda_tmm_cpm_log2_plus1) )
  
  
  # now correlating the metabolite with the gene expression data and selecting the gene have the highest and the lowest correlation
  
  # now selecting the common cell lines between the gene expression and the metabolite abundance datasets
  
  intersecting_cell_lines <- intersect( rownames(PM.data.normalise.significant.metabolites), colnames(annotated_pda_tmm_cpm_log2_plus1))
  
  data2plot <- annotated_pda_tmm_cpm_log2_plus1[ annotated_pda_tmm_cpm_log2_plus1[,1] %in% metabolic_genes,  ]
  data2plot.1 <- apply(data2plot[,-1], 2, as.numeric)
  #rowmed <- apply( data2plot.1, 1, median)
  
  data2plot.2 <- data2plot.1[, colnames(data2plot.1) %in%  intersecting_cell_lines] #- rowmed
  data2plot.2.sd <- apply(data2plot.2, 1, var)
  rownames(data2plot.2) <- data2plot[,1]
  
  names(data2plot.2.sd) <- data2plot[,1]
  var_genes <- sort(data2plot.2.sd, decreasing = TRUE)
  top_1000_var_genes <- names(var_genes)[1:1000]
  
  data2cluster.plot <- data2plot.2[top_1000_var_genes, ]
  #data2cluster.plot <- data2plot.2[top_1000_var_genes, ]
  
  
  # selecting the top variable genes
  d1 <- PM.data.normalise.significant.metabolites[ intersecting_cell_lines, ]
  d2 <- t(data2cluster.plot)
  d3 <- d2[ intersecting_cell_lines, ]
  
  correlation_score <- cor( d1, d3, method = "spearman")
  df.correlation <- melt(correlation_score)
  
  # seleting the to highly correlated genes and metabolites
  ind.sel <- which( df.correlation$value>=0.8 | df.correlation$value <= -0.8)
  
  df.correlation.sel <- df.correlation[ind.sel, ]
  colnames(df.correlation.sel) <- c("metabolites", "genes","spearman_correlation")
  
  df.correlation.sel <- df.correlation.sel[ order( df.correlation.sel$spearman_correlation), ]
  
  my_order <- unique(df.correlation.sel$genes)
  

  df.correlation.sel$order <-factor(df.correlation.sel$genes, levels=my_order)
  
  corr_file_all <- paste0(output_dir,"/", Sys.Date(), plate_num, "_Gene_expression_PDA_cell_lines_metabolic_genes_spearman_correlation_all_values.txt")
  
  corr_file <- paste0(output_dir,"/", Sys.Date(), plate_num, "_Gene_expression_PDA_cell_lines_metabolic_genes_spearman_correlation_cutoff0.8.txt")
  
  
  write.table(df.correlation.sel, file = corr_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(df.correlation, file = corr_file_all, sep = "\t", row.names = FALSE, quote = FALSE)
  
  
  # now plotting the correlation score
  pp <- ggplot(data = df.correlation.sel, aes(y=metabolites, x=order, fill = spearman_correlation))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") + theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed()
  
  
  
  
  ggsave( pp, filename = paste0(output_dir,"/", Sys.Date(), plate_num, "_Gene_expression_PDA_cell_lines_metabolic_genes_spearman_correlation_cutoff0.8.pdf"))
  
  
}



gene_metabolite_correlation_all_transcriptome <- function(plate_num, pm_sig_metabolite_rda, output_dir){
  
  
  ########## Whole transcriptomic analysis ###############################
  
  load(pm_sig_metabolite_rda)
  
  
  # reading the cell line data
  cell_line_file <- "/Users/ppoudel/Dropbox (ICR)//Pawan/DATA/External/Pancreas/Cell lines/KlijnEtAl2014/Rdata/20170509_PP_KlijnEtAl2014_coding_genes_tmm_cpm_log2_plus1_filtered_annotated.Rdata"
  load(cell_line_file)
  
  colnames(annotated_pda_tmm_cpm_log2_plus1) <- toupper( colnames(annotated_pda_tmm_cpm_log2_plus1))
  colnames(annotated_pda_tmm_cpm_log2_plus1) <- gsub( "PATU", "", colnames(annotated_pda_tmm_cpm_log2_plus1) )
  colnames(annotated_pda_tmm_cpm_log2_plus1) <- gsub( "PANC0327", "PANC327", colnames(annotated_pda_tmm_cpm_log2_plus1) )
  
  # now doing this for the whole genes in the transcriptome
  
  whole_gene_exp_data <- apply( annotated_pda_tmm_cpm_log2_plus1[,-1], 2, as.numeric )
  #row.med <- apply( whole_gene_exp_data, 1, median)
  whole_gene_exp_data.f <- whole_gene_exp_data #- row.med
  rownames(whole_gene_exp_data.f) <- annotated_pda_tmm_cpm_log2_plus1[,1]
  
  intersecting_cell_lines <- intersect( rownames(PM.data.normalise.significant.metabolites), colnames(whole_gene_exp_data.f))
  
  d1 <- PM.data.normalise.significant.metabolites[ intersecting_cell_lines, ]
  d2 <- t(whole_gene_exp_data.f)
  d3 <- d2[ intersecting_cell_lines, ]
  
  correlation_score <- cor( d1, d3, method = "spearman")
  df.correlation <- melt(correlation_score)
  
  # seleting the to highly correlated genes and metabolites
  ind.sel <- which( df.correlation$value>=0.85 | df.correlation$value <= -0.85)
  
  df.correlation.sel <- df.correlation[ind.sel, ]
  colnames(df.correlation.sel) <- c("metabolites", "genes","spearman_correlation")
  
  df.correlation.sel <- df.correlation.sel[ order( df.correlation.sel$spearman_correlation), ]
  
  my_order <- unique(df.correlation.sel$genes)
  
  df.correlation.sel$order <-factor(df.correlation.sel$genes, levels=my_order)
  
  
  corr_file_all <- paste0(output_dir,"/", Sys.Date(), plate_num, "_Gene_expression_PDA_cell_lines_entire_transcriptome_spearman_correlation_all_values.txt")
  corr_file <- paste0(output_dir,"/", Sys.Date(), plate_num, "_Gene_expression_PDA_cell_lines__entire_transcriptome_spearman_correlation_cutoff0.85.txt")
  
  write.table(df.correlation.sel, file = corr_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(df.correlation, file = corr_file_all, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # now plotting the correlation score
  
  #ind <- round( length(df.correlation.sel$metabolites)/3, digits = 0)
  
  
  pp <- ggplot(data = df.correlation.sel, aes(y=metabolites, x=order, fill = spearman_correlation))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") + theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed()
  
    
  
  ggsave( pp, filename = paste0(output_dir,"/", Sys.Date(), plate_num, "_Gene_expression_PDA_cell_lines_all_transcriptome_spearman_correlation_cutoff0.85.pdf"))
  
  
}


# creating a x-y plot for the high correlation coefficients





