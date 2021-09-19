source("/cs/icore/ekushele/script/r/filippo/functions_atlas.R")
output_dir="/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/new_atlas/results/"

#create methylation score table
pvalue=0.001
print(pvalue)
diff_mean=0.3


#get cpgs unique to LUAD-get normal lung CpGs
num_cpg <- 2000
print(paste("Number initial CpGs in atlas:",num_cpg))
#get josh data
best_hypo_hyper_probes_file <- paste0( "/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/josh_moss_probes/creating_new_atlas_scripts/CpGs.100bp-block.",num_cpg,".X1.xls")
best_hypo_hyper_probes <- data.frame(fread(best_hypo_hyper_probes_file))[,c("acc","group","name")]
best_hypo_hyper_probes$name <- gsub(" ","_",best_hypo_hyper_probes$name)
best_probes_tissue <-subset(best_hypo_hyper_probes,name=="Lung_cells")


#####=========score of methylation according to nanopore samples
data_on_files <- function(cur_file,sample_mapped,sample_name,tissue_file) {
	
	#score will be NA if data not exists (mainly LAML)
	score_of_cpg <- nnls_score_of_cpg <- NA
	cur_file <-cur_file[cur_file$FDR_adjust < pvalue & abs(cur_file$mean_diff) > diff_mean & !is.na(cur_file$FDR_adjust), ]
	if (nrow(cur_file)>0) {
		CpGs <- cur_file$CpGs_common
		#CpGs <- setdiff(CpGs,best_hypo_hyper_probes[best_hypo_hyper_probes$name=="Lung_cells","acc"])

		sample_mapped$cpg_index <- sub("^","cg",sample_mapped$cpg_index)

		sample_overlap_meth <- sample_mapped[which(sample_mapped$cpg_index %in% intersect(sample_mapped$cpg_index,CpGs)),c("cpg_index","V5","V6","V7", "methylation")]
		#CpGs includes CpGs that repeats themselved in sample_mapped
		CpGs = sample_overlap_meth$cpg_index
		
		CpGs_numbers_rep <- plyr::count(CpGs)
		CpGs_numbers_rep <- CpGs_numbers_rep[match(unique(CpGs),CpGs_numbers_rep$x),]
		cur_file <-cur_file[which(cur_file$CpGs_common %in% intersect(CpGs,cur_file$CpGs_common)),]

	 	cur_file <- cur_file[match(unique(CpGs),cur_file$CpGs_common),]
	 	cur_file <- cur_file[rep(1:nrow(cur_file), times=CpGs_numbers_rep$freq),]
		# sample_overlap_meth <- sample_overlap_meth[match(cur_file$CpGs_common,sample_overlap_meth$cpg_index),]
		score_of_cpg <- (sample_overlap_meth$methylation - cur_file[,"mean_healthy"])/(cur_file$mean_tumor-cur_file$mean_healthy)
		nnls_score_of_cpg <- nnls::nnls(as.matrix(data.frame(cur_file$mean_tumor,cur_file$mean_healthy)),sample_overlap_meth$methylation)$x
	}
	#write file for ben  
	# if (grepl("LUAD",tissue_file)) 
	if (nrow(cur_file)>0) 

		fwrite(data.frame(chr=sample_overlap_meth$V5,start=sample_overlap_meth$V6,end=sample_overlap_meth$V7,mean_cancer=cur_file$mean_tumor,mean_healthy=cur_file$mean_healthy,methylation=sample_overlap_meth$methylation,p_value=cur_file$p_val,FDR=cur_file$FDR_adjust),file.path(output_dir,"atlas_vs_healthy","check",paste0(sample_name,"_",sub(".csv","",basename(tissue_file)),"_","check.csv")))
	data.frame(mean_data_hypo_hyper=mean(score_of_cpg),nnls_hypo_hyper=nnls_score_of_cpg[[1]]/sum(nnls_score_of_cpg),length_overlapping_meth=length(score_of_cpg))
	
}

#find the hypo-hyper methylation for each sample
all_t_test_files <- dir(file.path(output_dir,"atlas_vs_healthy","t_test"),pattern="csv",full.names = T)
all_corrected_TCGA_files <- dir("/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/new_atlas/TCGA_purity/",pattern="csv",full.names = T)
samples_mapped_to_CpG <- dir(file.path("/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/samples_nanopore/","samples_merged_CpG_450K"),full.names=T,pattern="_shifted_mapped_CpG.bed.gz")

all_table_methylation_in_tissue <- do.call(rbind,lapply(samples_mapped_to_CpG,function(sample_bed) {
	print(sample_bed)
	sample_mapped <- data.frame(fread(sample_bed,colClasses = c('V4'='character')))
	colnames(sample_mapped)[4] <- "cpg_index"
	colnames(sample_mapped)[8] <- "methylation"
	#take only sample_mapped where there are CpGs with reasonable	 FDR exists
	do.call(rbind,parallel::mclapply(all_t_test_files,function(tissue_file) {
		print(tissue_file)
		cur_file <- data.frame(fread(tissue_file))
		if (nrow(cur_file)>0)  {
			cur_corrected_TCGA_file <- grep(basename(tissue_file),all_corrected_TCGA_files,value=T)
			cur_corrected_TCGA_file <- colMeans(data.frame(fread(cur_corrected_TCGA_file)),na.rm=T)
			cur_corrected_TCGA_file <- data.frame(pure_cpg =cur_corrected_TCGA_file) 
			cur_file <- merge(cur_file,cur_corrected_TCGA_file,by.x="CpGs_common",by.y="row.names")
			cur_file$mean_tumor <- cur_file$pure_cpg
		}
			scores_of_meth <- data_on_files(cur_file,sample_mapped,sub("\\..*","",basename(sample_bed)),tissue_file)
		data.frame(cancer_type=sub("\\.csv","",basename(tissue_file)),sample=basename(sample_bed),scores_of_meth)
	},mc.cores=future::availableCores()[[1]]))
}))

####====scoring
samples_healthy <- c("BC02","BC03","BC05","BC04","BC12")
 

all_table_methylation_in_tissue[is.na(all_table_methylation_in_tissue)] <- 0
all_table_methylation_in_tissue$category <- "Tumor"
all_table_methylation_in_tissue$category[grep(paste0(samples_healthy,collapse="|"),all_table_methylation_in_tissue$sample)] <- "Healthy"

types_of_raw_scores <- unique(sub("normalized.","",grep("hypo_hyper",colnames(all_table_methylation_in_tissue),value=T)))
normalized_data <- do.call(cbind,lapply(types_of_raw_scores,function(x) {
	medians <- aggregate(get(x)~category+cancer_type,all_table_methylation_in_tissue,median)[c(T,F),]
	all_table_methylation_in_tissue[,x]- rep(medians[,3],each=length(unique(all_table_methylation_in_tissue$sample)))
}))
colnames(normalized_data) <- paste0("normalized.",types_of_raw_scores)
all_table_methylation_in_tissue <- cbind(all_table_methylation_in_tissue[,grep("normalized",colnames(all_table_methylation_in_tissue),invert=T)],normalized_data)

#file_to_write <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization.csv"
#file_name_to_write <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization_by_beta_mean.csv"
#file_name_to_write <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization_by_median_healthy.csv"
# file_to_write <- "methylation_in_sample_per_tissue_score_table_all_samples_new_score.csv"
#file_to_write <- "methylation_in_sample_LUAD_score_table_all_samples_new_score_tumor_only.csv"
# file_to_write <- "methylation_in_sample_LUAD_score_table_all_samples_new_score_tumor_only_unique_to_LUAD.csv"
file_to_write <- "methylation_in_sample_per_tissue_table_all_samples_new_score_tumor_only.csv"
fwrite(all_table_methylation_in_tissue,file.path(output_dir,"atlas_vs_healthy",file_to_write))
file_name_to_link <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization.csv"

system(paste("ln -sf",file.path(output_dir,"atlas_vs_healthy",file_to_write),file.path(output_dir,"atlas_vs_healthy",file_name_to_link)))
