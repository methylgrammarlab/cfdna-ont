source("/cs/icore/ekushele/script/r/filippo/functions_atlas.R")
output_dir="/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/new_atlas/results/"


tumor_fractions <- data.frame(read.table(text=system("less -S /vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/tumor_fractions/tumorFractionIchorCNA_newRun_091921.txt| awk -F '\t' -v OFS='\t' '{print $1,$8,$16}'",intern=T),header=T))
colnames(tumor_fractions)[3] <- "num_reads"
tumor_fractions$sample <- sub("\\.fragmentomic.hg38","",tumor_fractions$sample)


file_to_write <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization.csv"
all_table_methylation_in_tissue <- data.frame(fread(file.path(output_dir,"atlas_vs_healthy",file_to_write)))
all_table_methylation_in_LUAD <- subset(all_table_methylation_in_tissue,cancer_type=="LUAD")
all_table_methylation_in_LUAD$sample <- sub("\\.hg38.*bed.*","",all_table_methylation_in_LUAD$sample)
all_table_methylation_in_LUAD <- merge(all_table_methylation_in_LUAD, tumor_fractions,by="sample")


illumina_data <- all_table_methylation_in_LUAD[grep("BC08",all_table_methylation_in_LUAD$sample),]
illumina_data$sample="BC08_illumina"
illumina_data$category="Tumor Illumina"
illumina_data$tumorFrac=0.109
illumina_data$num_reads=17076899
all_table_methylation_in_LUAD <- rbind(all_table_methylation_in_LUAD,illumina_data)

print(all_table_methylation_in_LUAD)

title_for_plot="TCGA Tumor 450k vs. Healthy Plasma 450k"
# output_dir_path <- file.path(output_dir,"atlas_vs_healthy","plots","tumor_fractions_vs_LUAD")
output_dir_TF <- file.path(output_dir,"atlas_vs_healthy","plots","tumor_fractions_vs_LUAD_maxCn3_0")
#output_dir_path <- file.path(output_dir,"atlas_vs_healthy","plots","tumor_fractions_vs_LUAD_maxCn3_0_subtract_mean(healthy)")
samples_to_remove_options <- list(c(),"BC09",c("BC08","BC09"))

lapply(samples_to_remove_options[1],function(sample_to_remove) {
	if (length(sample_to_remove) > 0) 
			all_table_methylation_in_LUAD[grep(paste0(sample_to_remove,collapse="|"),all_table_methylation_in_LUAD$sample),c("category")] <- "Low Content Tumor"
	removed_samples <- ifelse(length(sample_to_remove)>0, paste0("_no_",paste0(sample_to_remove,collapse="_")),"")
	output_dir_path <- file.path(output_dir_TF,paste0("TF",removed_samples))
	lapply(c(F,T),function(labeled) {
	lapply(c(F,T),function(identity_line) {
		suffix_identity=ifelse(isTRUE(identity_line),"_identity_line","")
		suffix=ifelse(isTRUE(labeled),paste0(suffix_identity,"_with_label.jpeg"),paste0(suffix_identity,".jpeg"))
		lapply(grep("hypo_hyper",colnames(all_table_methylation_in_LUAD),value=T)[1:2],function(type_of_score) 		
			plot_fraction_tumor_vs_data(all_table_methylation_in_tissue=all_table_methylation_in_LUAD,
									main_title="",#LUAD Methylation Score vs CNA tumor Fraction",
									sub_title="",#title_for_plot,
									x_axis_data="tumorFrac",y_axis_data=type_of_score,
									xlab="Tumor fractions (ichorCNA estimate)",ylab=expression(bold("Estimated lung fraction"~(beta))),
									sample_to_remove=sample_to_remove,plot_identity_line=identity_line,
									output_dir_path=output_dir_path,
									save_plot=T,suffix=suffix,add_sample_label=labeled,add_horizonal_line=T,add_arrows=F))

	}) 
	})
})

# #predicted tumor fractions vs tumor fractions
# output_dir_predicted_tumor <- file.path(output_dir,"atlas_vs_healthy","plots","predicted_tumor_frac")
# type_of_tumorFrac <- "Cancer.DNA.fraction"
# suffix=""
# # suffix="_with_normal"
# #get regression line from plots/TCGA_score/TCGA_score_vs_Cancer.DNA.fraction/TCGA_data_vs_Cancer.DNA.fraction_relative_weight_hypo_hyper.jpeg
# lapply(samples_to_remove_options,function(sample_to_remove) {
# 	if (length(sample_to_remove) > 0) 
# 			all_table_methylation_in_LUAD[grep(paste0(sample_to_remove,collapse="|"),all_table_methylation_in_LUAD$sample),c("category")] <- "Low Content Tumor"
# 	removed_samples <- ifelse(length(sample_to_remove)>0, paste0("_no_",paste0(sample_to_remove,collapse="_")),"")
# 	output_dir_path <- file.path(output_dir_predicted_tumor,paste0("predicted_tumor_frac_vs_LUAD"),paste0("TF_",removed_samples))
# 	lapply(unique(sub("normalized.","",grep("weight_hypo_hyper",colnames(all_table_methylation_in_LUAD),value=T))),function(x) {
# 		coefficients <- unlist(fread(file.path(output_dir,"atlas_vs_healthy/plots/TCGA_score/",paste0("TCGA_score_vs_",type_of_tumorFrac),paste0("coefficients_",x,suffix))))
# 		all_table_methylation_in_LUAD$predictedTumorFrac <- coefficients[1]+coefficients[2]*all_table_methylation_in_LUAD[,x]
# 		plot_fraction_tumor_vs_LUAD(all_table_methylation_in_LUAD,type_of_score=x,atlas_title=title_for_plot,x_axis_data="predictedTumorFrac",xlab="Tumor Fractions (TCGA Methylation Estimate)",sample_to_remove,output_dir_path=output_dir_path,save_plot=T,suffix=paste0(suffix,".jpeg"))
# 	})
# })

# ###############3

# plot_tumorFractions <- function(all_table_methylation_in_LUAD,type_of_score,suffix) {
# 	all_table_methylation_in_LUAD <- change_name_of_label(all_table_methylation_in_LUAD,"Tumor","Cancers (others)")
# 	all_table_methylation_in_LUAD <- change_name_of_label(all_table_methylation_in_LUAD,"Low Content Tumor","Cancers (BC08, BC09)")
# 	groups_for_plot <- levels(as.factor(all_table_methylation_in_LUAD$category))
# 	colors_for_plot <-c( "lightpink2", "red","blue")

# 	 p <- plot_correlation(data_to_plot=all_table_methylation_in_LUAD,y_axis_data="predictedTumorFrac",x_axis_data="tumorFrac",main_title="Tumor Fractions",sub_title="TCGA Methylation Estimate vs. ichorCNA Estimate",xlab="Tumor Fraction (ichorCNA Estimate)",ylab="Tumor Fraction (TCGA Methylation Estimate)",color_by_category=T)+scale_color_manual(breaks = c(groups_for_plot),values=c(colors_for_plot))
# 	output_dir_path <- file.path(output_dir_predicted_tumor,"predictedTumorFrac_vs_tumorFrac")
# 	create_dir(output_dir_path)
# 	ggsave(file.path(output_dir_path,paste0("predictedTumorFrac_vs_tumorFrac_",type_of_score,suffix)), p,height=6,width=8)
# 	message(paste("saved",file.path(output_dir_path,paste0("predictedTumorFrac_vs_tumorFrac",type_of_score,suffix))))
# }

# type_of_score <- c("relative_weight_hypo_hyper","overlapping_relative_weight_hypo_hyper")
# lapply(c("", "_with_normal"), function(suffix) {
# 	lapply(type_of_score,function(x) {
# 		coefficients <- unlist(fread(file.path(output_dir,"atlas_vs_healthy/plots/TCGA_score/",paste0("TCGA_score_vs_",type_of_tumorFrac),paste0("coefficients_",x,suffix))))
# 		all_table_methylation_in_LUAD$predictedTumorFrac <- coefficients[1]+coefficients[2]*all_table_methylation_in_LUAD[,x]
# 		plot_tumorFractions(all_table_methylation_in_LUAD,x,paste0(suffix,".jpeg"))

# 	})
# })


