 # install.packages("conquer", repos='http://cran.us.r-project.org')
 # install.packages("ggpmisc", repos='http://cran.us.r-project.org')
 source("/cs/icore/ekushele/script/r/filippo/functions_atlas.R")
args = commandArgs(trailingOnly=TRUE)
output_dir="/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/new_atlas/results/"

num_selected_cpgs=as.numeric(args[1])*2


tumor_fractions <- data.frame(read.table(text=system("less -S /vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/tumor_fractions/tumorFractionIchorCNA_newRun_091921.txt| awk -F '\t' -v OFS='\t' '{print $1,$8,$16}'",intern=T),header=T))
colnames(tumor_fractions)[3] <- "num_reads"
tumor_fractions$sample <- sub("\\.fragmentomic.hg38","",tumor_fractions$sample)

#file_to_write <- paste0("methylation_in_sample_per_tissue_score_table_all_samples_new_score_",num_selected_cpgs/2,".csv")
file_to_write <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization.csv"
all_table_methylation_in_tissue <- data.frame(fread(file.path(output_dir,"atlas_joshCpGs_vs_healthy",file_to_write)))
print(file.path(output_dir,"atlas_joshCpGs_vs_healthy",file_to_write))
all_table_methylation_in_LUAD <- subset(all_table_methylation_in_tissue,cancer_type=="Lung_cells")
all_table_methylation_in_LUAD$sample <- sub("\\.hg38.*bed","",all_table_methylation_in_LUAD$sample)
all_table_methylation_in_LUAD <- merge(all_table_methylation_in_LUAD, tumor_fractions,by="sample")

illumina_data <- all_table_methylation_in_LUAD[grep("BC08",all_table_methylation_in_LUAD$sample),]
illumina_data$sample="BC08_illumina"
illumina_data$category="Tumor Illumina"
illumina_data$tumorFrac=0.109
illumina_data$num_reads=17076899
all_table_methylation_in_LUAD <- rbind(all_table_methylation_in_LUAD,illumina_data)

title_for_plot=paste("Atlas of",num_selected_cpgs,"selected CpGs")
output_dir_TF <- file.path(output_dir,"atlas_joshCpGs_vs_healthy",paste0("plots_",num_selected_cpgs),"tumor_fractions_vs_LUAD_maxCn3_0")
samples_to_remove_options <- list(c(),"BC09",c("BC08","BC09"))

lapply(samples_to_remove_options[1],function(sample_to_remove) {
	if (length(sample_to_remove) > 0) 
			all_table_methylation_in_LUAD[grep(paste0(sample_to_remove,collapse="|"),all_table_methylation_in_LUAD$sample)[1],c("category")] <- "Low Content Tumor"
	removed_samples <- ifelse(length(sample_to_remove)>0, paste0("_no_",paste0(sample_to_remove,collapse="_")),"")
	output_dir_path <- file.path(output_dir_TF,paste0("TF",removed_samples))
	lapply(c(F,T),function(type_of_line) {
		add_to_suffix=ifelse(isTRUE(type_of_line),"_with_identity","")
		lapply(c(F,T),function(labeled) {
			suffix=paste0(add_to_suffix,ifelse(isTRUE(labeled),"_with_label.jpeg",".jpeg"))
			   
			lapply(grep("nnls_hypo_hyper",colnames(all_table_methylation_in_LUAD),value=T),function(type_of_score) {
			plot_fraction_tumor_vs_data(all_table_methylation_in_tissue=all_table_methylation_in_LUAD,
										main_title="",#"Lung Methylation Score vs. CNA Tumor Fraction",
										sub_title="",#title_for_plot,
										x_axis_data="tumorFrac",y_axis_data=type_of_score,
										xlab="Tumor fractions (ichorCNA estimate)",ylab=expression(bold("Estimated lung fraction"~(beta))),
										sample_to_remove=sample_to_remove,plot_identity_line=type_of_line,
										output_dir_path=output_dir_path,
										save_plot=T,suffix=suffix,add_sample_label=labeled,add_horizonal_line=T,add_arrows=F)
			}) })
	})
})
