


rm(list=ls())
options(stringsAsFactors = FALSE)


source("/Users/ppoudel/Dropbox (SPM)/Pawan/Scripts/Phenotypic Microarrays/2017/source/PM_Arrays_source.R")

# for PM1 plate carbon source
data_dir <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/datasets/PM/20171505_PP_PM_Processed/CarbonPlates/PM-M01/"
output_fol <-  "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/analysis/opm_output/carbon_plates/"
plate_num="PM01"
outlierCelllines2Remove <- "CFPAC1"
rda_dir <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/rdata/"
preprocessPMplates(plate_num, data_dir, output_fol,  rda_dir,outlierCelllines2Remove )

# for PM2 plate carbon source
data_dir <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/datasets/PM/20171505_PP_PM_Processed/CarbonPlates/PM-M02/"
output_fol <-  "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/analysis/opm_output/carbon_plates/"
plate_num="PM02_Carbon"
outlierCelllines2Remove <- "CFPAC1"
rda_dir <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/rdata/"
preprocessPMplates(plate_num, data_dir, output_fol,  rda_dir,outlierCelllines2Remove )


# for PM2 plate nitrogen source
data_dir <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/datasets/PM/20171505_PP_PM_Processed/NitrogenPlates/PM-M02/"
output_fol <-  "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/analysis/opm_output/nitrogen_plates/"
plate_num="PM02_Nitrogen"
outlierCelllines2Remove <- ""
rda_dir <- "/Users/ppoudel/Dropbox (SPM)/Pawan/ANALYSIS/Internal/Pancreas/Biolog/rdata/"
preprocessPMplates(plate_num, data_dir, output_fol,  rda_dir,outlierCelllines2Remove )


preprocessPMplates <- function(plate_num, data_dir, output_fol, rda_dir,outlierCelllines2Remove="" ){
  
  # loading the packages
  require(stringr)
  require(opm)
  require(reshape)
  require(ggplot2)
  require(gplots)
  
  #plate_num <- "PM-M02"
  #plate_num_all <- c("PM-M01", "PM-M02", "PM-M03", "PM-M04", "PM-M05", "PM-M06", "PM-M07", "PM-M08", "PM-M11", "PM-M12", "PM-M13", "PM-M14")
  
  
  metadata_dir <- paste0(output_fol, "/opm_output_metadata/", plate_num)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_dir <- paste0(output_fol, "/opm_output/", plate_num)
  dir.create(out_dir, showWarnings = FALSE)
  
  # reading one plate at a time
  csvs <- list.files(path=data_dir, pattern = ".csv", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
  
  # before reading the csvs check for the missing metabolites
  
  check_metabolites <- lapply( csvs, function(x){
    
    # reading the each csv files to check for the missing values
    my_data <- read_single_opm(x)
    my_data_rows <- nrow(my_data@measurements)
    my_data_cols <- ncol(my_data@measurements) - 1  # ignoring the first columns
    
    df <- data.frame(x, my_data_rows, my_data_cols )
    return(df)
    
  } )
  
  qc_check <- do.call("rbind", check_metabolites)
  write.table(qc_check, paste0(metadata_dir, "/", Sys.Date(), "_", plate_num, "checking_data_generated_raw_data.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  setwd(data_dir)
  
  # reading the opm data
  ff <- read_opm( convert = "yes")
  
  
  
  # aggregating the curve parameters
  data_A <- do_aggr(ff, method="splines" )
  
  
  # creatin the metadata
  # this point you can add the subtype information to a metadata file. However, in this case I am just creating the mock cell lines or the real cell lines
  my_temp <- collect_template(data_A, outfile=paste0(metadata_dir, "/", plate_num,"_metadata.csv"))
  
  # reading the metadata
  metadata.example <- to_metadata(paste0(metadata_dir, "/", plate_num,"_metadata.csv"))
  metadata.example$combo <- paste0( metadata.example$cell_lines, "_", metadata.example$Position, "_", metadata.example$replicates )
  metadata.example$cell_lines <- do.call( "rbind", strsplit( metadata.example$File, split = "/"))
  metadata.example$cell_lines <- metadata.example$cell_lines[, ncol(metadata.example$cell_lines)]
  metadata.example$cell_lines <- toupper( gsub( "_PM-.+|\\.", "" , metadata.example$cell_lines ) )
  

  # adding the metadata and combining the replicates
  example.opm <- include_metadata(data_A, md = metadata.example)
  
 
  # extracting the individual curve parameters
  data.mu <- extract(example.opm, dataframe = TRUE, as.labels = list("cell_lines") , subset = "mu")
  data.lambda <- extract(example.opm, dataframe = TRUE, as.labels = list("cell_lines"), subset = "lambda")
  data.A <- extract(example.opm, dataframe = TRUE, as.labels = list("cell_lines"), subset = "A")
  data.AUC <- extract(example.opm, dataframe = TRUE, as.labels = list("cell_lines"), subset = "AUC")
  
  # writing the values
  out_dir_A <- paste0(output_fol, "/opm_output/", plate_num, "/A_OPM/")
  dir.create(out_dir_A, showWarnings = FALSE, recursive = TRUE)
  
  write.table( data.A, file = paste0(out_dir_A ,"/", Sys.Date(), "_OPM_A_values_", plate_num, ".txt"), sep = "\t", quote = FALSE )
  
  
  # based on the paper now we are comparing the results from the A, AUC, lambda and mu
  dir.create( paste0( out_dir,"/plots/"), showWarnings = FALSE, recursive = TRUE)
  plot_dir <- paste0( out_dir,"/plots/")
  
  setwd(plot_dir)
  
  # plotting the heatmap
  dir.create( paste0(plot_dir, "/", "heatmap/"), showWarnings = FALSE)
  setwd( paste0(plot_dir, "/", "heatmap/"))
  
  plot_heatmap(data.mu, num = plate_num)
  plot_heatmap(data.lambda, num = plate_num)
  plot_heatmap(data.A, num = plate_num)
  plot_heatmap(data.AUC, num = plate_num)
  
  
  # plotting the correlation between the curve parameters
  dir.create( paste0(plot_dir, "/", "correlation_param/"), showWarnings = FALSE)
  setwd( paste0(plot_dir, "/", "correlation_param/"))
  
  plot_coor(  data.A,data.lambda, num = plate_num)
  plot_coor(  data.A,data.AUC, num = plate_num)
  plot_coor(  data.A,data.mu, num = plate_num)
  plot_coor(  data.AUC, data.mu, num = plate_num)
  plot_coor(  data.AUC, data.lambda, num = plate_num)
  # again coming back to the plot directory
  
  
  setwd(plot_dir)
  
  # now comapring the variabiltiy between the plates using the PCA
  # before removing the negetive controls
  dat1 <- data.A[, -c(1:2)]
  rownames(dat1) <- data.A[,1]
  
  create.PCA.Plot(mydata = dat1, plot_file = paste0(Sys.Date(), "_PCA_plot_All_data_selecting_A_Cell_lines_", plate_num, ".pdf") )
  create.PCA.Plot(mydata = t(dat1), plot_file = paste0(Sys.Date(), "_PCA_plot_All_data_selecting_A_metabolites_", plate_num, ".pdf") )
  
  # plotting the distribution of negetive controls for all the cell lines
  df.negetive.control <- melt( data.A[, c(1,3:5)])
  colnames(df.negetive.control) <- c("cell_lines", "variable", "value")
  g <- ggplot( df.negetive.control, aes(x=cell_lines, y=value, fill=variable )) +
    geom_point(alpha=0.6,  position = "jitter", aes(colour=variable), size = 3, stroke = 3) + ylim(0,350) + theme_Publication()
  
  ggsave( g,  filename = paste0( plot_dir, "/", Sys.Date(),"_",plate_num, "_negeative_controls_all_cell_lines_points.pdf") )
  
  
  
  # now plotting the distribution of negetive controls
  
  
  # after removing the CFPAC1 and perfroming the analysis again
  if(length(outlierCelllines2Remove) >=1 ){
    
    dat2 <- dat1[ !(rownames(dat1) %in% outlierCelllines2Remove ), ]
    create.PCA.Plot(mydata = dat2, plot_file = paste0(Sys.Date(), "_PCA_plot_All_data_selecting_A_ignoring_metabolites", plate_num, ".pdf") )
    # just removing the CFPAC1 from the analysis
    data.A <- data.A[ !(data.A$cell_lines %in% outlierCelllines2Remove), ]
    
  }else{
    
    dat2=dat1
  
  }
  
  
  # now performing the statistical test to select the high confidence metabolites
  
  dir.create( paste0( plot_dir, "/statistical_test/") , showWarnings = FALSE)
  setwd( paste0( plot_dir, "/statistical_test/") )
  
  data.melt <- melt(data.A)
  my_met <- colnames(data.A)[-c(1:5)]
  
  select_significant_metablites <- lapply(my_met, function(x, df=data.melt){
    
    pf1 <- paste0("Comparing negetive control vs met ", x, " all cell lines.pdf")
    pf2 <- paste0("Comparing negetive control vs met ", x, " differentail A value.pdf")
    
    df[,3] <- as.character(df[,3])
    ind1 <- grep("A01|A02|A03", df[,3])
    x1 <- gsub("\\(.+| ", "", x)
    
    ind2 <- grep(x1, df[,3])
    
    data2plot <- df[c(ind1,ind2), ]
    ind3 <- grep(x1, data2plot[,3])
    
    
    data2plot$value <- as.numeric(data2plot$value)
    #data2plot$variable <- gsub(",", "", data2plot$variable)
    
    f1 <- ggplot(data2plot, aes(x=cell_lines, y=value, colour= factor(cell_lines) ) ) + geom_boxplot() + facet_wrap(~ variable) + theme_bw() + theme(axis.text.y=element_text(hjust=1, angle=0), axis.text.x = element_text(hjust=1, angle=45))
    
    ggsave(f1, filename = pf1)
    
    data2plot$groups <- "Negative"
    data2plot$groups[ind3] <- x1
    
    pv_kruskal <-kruskal.test(value ~ factor(variable), data = data2plot)$p.value
    
    # also performing one tailed t-test
    pv_ttest <-t.test(value ~ factor(groups), data = data2plot, alternative = "greater")$p.value
    
    f2 <- ggplot(data2plot, aes(x=variable, y=value, colour= factor(variable) ) ) + geom_boxplot()
    f2 <- f2 + annotate("text", label = paste0( x," pval: ", pv_kruskal, " pval ttest", pv_ttest ), x=1, y=0) + theme_Publication() #+ theme_bw() + theme(axis.text.y=element_text(hjust=1, angle=0), axis.text.x = element_text(hjust=1, angle=45))
    ggsave(f2, filename = pf2)
    
    kruskall_test <- data.frame(x, pv_kruskal, pv_ttest)
    
    return(kruskall_test)
    
  })
  
  setwd(plot_dir)
  
  # merging the data in list to a data frame
  my_metabolic_signals_kruskall <- do.call("rbind", select_significant_metablites)
  colnames(my_metabolic_signals_kruskall) <- c("metabolites", "kruskall", "one_sided_test")
  
  f2_file <- paste0(out_dir, "/", Sys.Date(), "_", plate_num ,"_Negetive_Controls_vs_metabolites_Kruskall.pval_table.txt")
  
  # now performing the FDR adjustments
  my_metabolic_signals_kruskall$fdr.kruskal.test <- p.adjust(my_metabolic_signals_kruskall$kruskall, method = "fdr")
  my_metabolic_signals_kruskall$fdr.one.sided.ttest <- p.adjust(my_metabolic_signals_kruskall$one_sided_test, method = "fdr")
  # also doing this for the t-test results
  
  my_metabolic_signals_kruskall <- my_metabolic_signals_kruskall[ order(my_metabolic_signals_kruskall$fdr.kruskal.test, my_metabolic_signals_kruskall$fdr.one.sided.ttest), ]
  
  # adding the full names to the metabolites
  write.table(my_metabolic_signals_kruskall, f2_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # just getting the list of the significant metabolites
  significant_metabolites <- my_metabolic_signals_kruskall[ which(my_metabolic_signals_kruskall$fdr.one.sided.ttest <= 0.05), ]
  
  
  # plotting the PCA plot
  
  # normalising with the negetive control and plotting the correlation plot of negetive controls
  data.normalise <- step_before_normalise_data(data.mat = data.A, num = "PM1")
  
  # substracting the intensity of the negetive controls with the other metabolites
  require(gridExtra)
  
  data.normalise.significant.metabolites <- data.normalise[, colnames(data.normalise) %in% significant_metabolites$metabolites]
  PM.data.normalise.significant.metabolites <- data.normalise[, colnames(data.normalise) %in% significant_metabolites$metabolites]
  
  # after substracting the negetive controls
  create.PCA.Plot(mydata = data.normalise, plot_file = paste0(Sys.Date(), "_PCA_plot_All_data_substracting_negetive_probes_", plate_num, ".pdf") )
  
  
  # after selecting the significant metaboloites
  create.PCA.Plot(mydata = data.normalise.significant.metabolites, plot_file = paste0(Sys.Date(), "_PCA_plot_All_data_substracting_negetive_probes_significant_metabolites_Cell_lines_", plate_num, ".pdf") )
  
  
  
  # after selecting the significant metaboloites and plotting the metabolites not the cell lines
  create.PCA.Plot(mydata = t(data.normalise.significant.metabolites), plot_file = paste0(Sys.Date(), "_PCA_plot_All_data_substracting_negetive_probes_significant_metabolites_", plate_num, ".pdf") )
  
  
  # now considering the significant metabolites for the downstream anaysis
  
  # writing to the file
  out_file1 <- paste0( out_dir, "/" , Sys.Date(), "_", plate_num, "_selected_normalised_metabolies.txt")
  write.table(data.normalise.significant.metabolites, out_file1 , sep = "\t", quote = FALSE)
  
  # saving the data as Rdata
  save( PM.data.normalise.significant.metabolites, file = paste0(rda_dir, "/" , plate_num, "_Interesting_Metabolites.RData"))
  
  
  # creating a heatmap after normalising and selecting high confidence metabolites
  require(gplots)
  out_file_heatmap <- paste0( out_dir, "/" , Sys.Date(), "_", plate_num, "_selected_normalised_metabolies.pdf")
  pdf(out_file_heatmap)
  heatmap.2( PM.data.normalise.significant.metabolites, trace = "none", scale = "none", col=bluered)
  dev.off()
  
  # creating the additional plots
  
  require(RColorBrewer)
  n <- length(unique(metadata.example$cell_lines))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  my.col <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
  # now performing the QC of the data
  
  # now producing the x-y plots
  setwd(out_dir)
  
  inf="cell_lines"

  pdf( paste0(Sys.Date(), "_", plate_num, "_selected_metabolites_exploratory_analysis_plot_xy_plot_subtypes.pdf"))
  my_met <- c("A01","A02", "A03" ,gsub(" \\(.+", "", colnames(data.normalise.significant.metabolites) ))
  lapply(my_met, function(x, opm=example.opm, pl_n=plate_num){
    
    cat(x, "\n")
    xy_plot(example.opm[,,as.character(x)], include = inf) 
    
  })
  dev.off() 
  
#  pdf(paste0(out_dir, "/", Sys.Date(), "_", plate_num, "_selected_metabolites_exploratory_analysis_plot_level_plot_subtypes.pdf"))
#  level_plot(example.opm[,,my_met], include=inf)
#  dev.off()
  
  x <- extract(example.opm[,,my_met], as.labels = "cell_lines", dataframe = TRUE)
  
#  # plotting the radial plot
#  pdf(paste0(out_dir, "/", Sys.Date(), "_", plate_num, "_selected_metabolites_exploratory_analysis_radial_plot_subtypes.pdf"))
#  
#  
#  for(f in c(1,10,20,30,40)){
    
    # resetting f
#    if(f >= nrow(x)){
#      f <- nrow(x) - 10
#    }
#    radial_plot(x[f: (f+10), ], as.labels = c("cell_lines"), main = "Test")
    
#  }
  
#  dev.off()
  
  # plotting the ci plot
  #pdf(paste0(out_dir, "/", Sys.Date(), "_", plate_num, "_selected_metabolites_exploratory_analysis_plot_ci_plot_subtypes.pdf"))
  #ci_plot( example.opm[,,my_met], as.labels="cell_lines" )
  #dev.off()
  
  # plotting the heatmap
  pdf(paste0(out_dir, "/", Sys.Date(), "_", plate_num, "_selected_metabolites_exploratory_analysis_plot_heatmap_subtypes.pdf"))
  heat_map(example.opm[,,my_met], as.labels="cell_lines" )
  dev.off()
  
  
  
  
}




















