source("/sci/home/ekushele/icore-home/script/r/filippo/functions_atlas.R")
output_dir="/sci/labs/bermanb/ekushele/icore-data/filippo/data/new_atlas/results/"
args = commandArgs(trailingOnly=TRUE)
atlas_type <- "atlas_joshCpGs_vs_healthy"

argsLen <- length(args);

cancer_type_to_filter="Lung_cells"
suffix_general <- args[1]

if (argsLen <3) 
{	
	atlas_type <- sub("joshCpGs_","",atlas_type)
	cancer_type_to_filter <- "LUAD"
}
print(atlas_type)
tumor_fractions <- data.frame(read.table(text=system(paste0("less -S ",args[2],"| awk -F '\t' -v OFS='\t' '{print $1,$8,$16}'"),intern=T),header=T))
colnames(tumor_fractions)[3] <- "num_reads"
tumor_fractions$sample <- sub("\\.fragmentomic.hg38","",tumor_fractions$sample)
file_to_write <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization.csv"
all_table_methylation_in_tissue <- data.frame(fread(file.path(output_dir,atlas_type,file_to_write)))
all_table_methylation_in_LUAD <- subset(all_table_methylation_in_tissue,cancer_type==cancer_type_to_filter)
# all_table_methylation_in_LUAD$sample <- sub("\\.hg38.*bed.*","",all_table_methylation_in_LUAD$sample)
all_table_methylation_in_LUAD$sample <- sub("\\.remora1.*","",all_table_methylation_in_LUAD$sample)
all_table_methylation_in_LUAD <- merge(all_table_methylation_in_LUAD, tumor_fractions,by="sample",all.x=T)

all_table_methylation_in_LUAD$tumorFrac[all_table_methylation_in_LUAD$category=="Healthy"] <- 0
all_table_methylation_in_LUAD$num_reads[grep("remora",all_table_methylation_in_LUAD$sample)] <- 2800000

all_table_methylation_in_LUAD[grep("HU005",all_table_methylation_in_LUAD$sample),"category"]="Healthy(HUJI)"
all_table_methylation_in_LUAD[grep("Healthy$",all_table_methylation_in_LUAD$category),"category"]="Healthy(ISPRO)"

print(all_table_methylation_in_LUAD)

title_for_plot="TCGA Tumor 450k vs. Healthy Plasma 450k"
# output_dir_path <- file.path(output_dir,"atlas_vs_healthy","plots","tumor_fractions_vs_LUAD")
output_dir_TF <- file.path(output_dir,atlas_type,"plots","methylation_tumor_fractions")
#output_dir_path <- file.path(output_dir,"atlas_vs_healthy","plots","tumor_fractions_vs_LUAD_maxCn3_0_subtract_mean(healthy)")
samples_to_remove_options <- list(c(),"BC09",c("BC08","BC09"))


lapply(samples_to_remove_options[1],function(sample_to_remove) {
	if (length(sample_to_remove) > 0) 
			all_table_methylation_in_LUAD[grep(paste0(sample_to_remove,collapse="|"),all_table_methylation_in_LUAD$sample),c("category")] <- "Low Content Tumor"
	removed_samples <- ifelse(length(sample_to_remove)>0, paste0("_no_",paste0(sample_to_remove,collapse="_")),"")
	output_dir_path <- file.path(output_dir_TF,paste0("TF",removed_samples))
	lapply(c(F,T),function(labeled) {
		suffix=paste0(suffix_general,ifelse(isTRUE(labeled),"_with_label.pdf",".pdf"))
		lapply(grep("hypo_hyper",colnames(all_table_methylation_in_LUAD),value=T)[2],function(type_of_score) 		
			plot_methylation_tumor_fractions_for_catergory(all_table_methylation_in_tissue=all_table_methylation_in_LUAD,
									main_title="LUAD Methylation Score",
									sub_title=title_for_plot,
									y_axis_data=type_of_score,
									xlab="",ylab=expression(bold("Estimated lung fraction"~(beta))),
									sample_to_remove=sample_to_remove,
									output_dir_path=output_dir_path,
									save_plot=T,suffix=suffix,add_sample_label=labeled,levels_list=c("Healthy(ISPRO)","Healthy(HUJI)","Cancer")))

	}) 
	
})
