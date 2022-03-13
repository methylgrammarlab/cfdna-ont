source("/sci/labs/bermanb/ekushele/icore-home/script/r/filippo/functions_atlas.R")
output_dir="/sci/labs/bermanb/ekushele/icore-data/filippo/data/new_atlas/results/"

#create methylation score table
pvalue=0.001
print(pvalue)
diff_mean=0.3
args = commandArgs(trailingOnly=TRUE)

 
#get cpgs unique to LUAD-get normal lung CpGs
num_cpg <- as.numeric(args[4]) 
print(paste("Number initial CpGs in atlas:",num_cpg))
#get josh data
best_hypo_hyper_probes_file <- paste0( "/sci/labs/bermanb/ekushele/icore-data/filippo/data/josh_moss_probes/creating_new_atlas_scripts/CpGs.100bp-block.",num_cpg,".X1.xls")
best_hypo_hyper_probes <- data.frame(fread(best_hypo_hyper_probes_file))[,c("acc","group","name")]
best_hypo_hyper_probes$name <- gsub(" ","_",best_hypo_hyper_probes$name)
best_probes_tissue <-subset(best_hypo_hyper_probes,name=="Lung_cells")


#####=========score of methylation according to nanopore samples
data_on_files <- function(cur_file,sample_mapped,sample_name,tissue_file) {
	
	#score will be NA if data not exists (mainly LAML)
	score_of_cpg <- nnls_score_of_cpg <- NA
	
	#this is to remove any probes with only zero or one values
	cur_file <-cur_file[cur_file$FDR_adjust < pvalue & abs(cur_file$mean_diff) > diff_mean & !is.na(cur_file$FDR_adjust), ]
	if (nrow(cur_file)>0) {
		CpGs <- cur_file$CpGs_common
		CpGs <- setdiff(CpGs,best_hypo_hyper_probes[best_hypo_hyper_probes$name=="Lung_cells","acc"])

		sample_mapped$cpg_index <- sub("^","cg",sample_mapped$cpg_index)

		# sample_overlap_meth <- sample_mapped[which(sample_mapped$cpg_index %in% intersect(sample_mapped$cpg_index,CpGs)),c("cpg_index","V5","V6","V7", "methylation")]
		sample_overlap_meth_coordinates <- sample_mapped[which(sample_mapped$cpg_index %in% intersect(sample_mapped$cpg_index,CpGs)),c("cpg_index","V5","V6","V7", "methylation")]
		sample_overlap_meth <- aggregate(methylation~cpg_index,data=sample_overlap_meth_coordinates,FUN=mean, na.rm=TRUE, na.action=na.pass)
	
		#CpGs includes CpGs that repeats themselved in sample_mapped
		#sample_overlap_meth <- aggregate(.~cpg_index,data=sample_overlap_meth[,-c(2,3,4)],FUN=mean, na.rm=TRUE, na.action=na.pass) 

		CpGs = sample_overlap_meth$cpg_index
		
		# CpGs_numbers_rep <- plyr::count(CpGs)
		# CpGs_numbers_rep <- CpGs_numbers_rep[match(unique(CpGs),CpGs_numbers_rep$x),]
		cur_file <- cur_file[which(cur_file$CpGs_common %in% intersect(CpGs,cur_file$CpGs_common)),]
	 	cur_file <- cur_file[match(unique(CpGs),cur_file$CpGs_common),]
	 	#cur_file <- cur_file[rep(1:nrow(cur_file), times=CpGs_numbers_rep$freq),]
		# sample_overlap_meth <- sample_overlap_meth[match(cur_file$CpGs_common,sample_overlap_meth$cpg_index),]
		sample_overlap_meth <- merge(sample_overlap_meth,aggregate(.~cpg_index,sample_overlap_meth_coordinates[,-c(5)],head,1),by="cpg_index")

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
all_corrected_TCGA_files <- dir("/sci/labs/bermanb/ekushele/icore-data/filippo/data/new_atlas/TCGA_purity/",pattern=".csv",full.names = T)
#samples_mapped_to_CpG <- dir(file.path("/sci/labs/bermanb/ekushele/icore-data/filippo/data/samples_nanopore/","samples_merged_CpG_450K"),full.names=T,pattern="_shifted_mapped_CpG.bed.gz")
# samples_mapped_to_CpG <- dir(file.path("/sci/labs/bermanb/ekushele/icore-data/filippo/data/samples_nanopore_020222/","samples_merged_CpG_450K_links_060222"),full.names=T,pattern="_shifted_mapped_CpG.bed.gz")
samples_mapped_to_CpG <- dir(args[1],full.names=T,pattern="_shifted_mapped_CpG.bed.gz")
 
all_table_methylation_in_tissue <- do.call(rbind,lapply(samples_mapped_to_CpG,function(sample_bed) {
	print(sample_bed)
	sample_mapped <- data.frame(fread(sample_bed,colClasses = c('V4'='character')))
	colnames(sample_mapped)[4] <- "cpg_index"
	colnames(sample_mapped)[8] <- "methylation"
	#take only sample_mapped where there are CpGs with reasonable	 FDR exists
	do.call(rbind,parallel::mclapply(all_t_test_files[17],function(tissue_file) {
		print(tissue_file)
		cur_file <- data.frame(fread(tissue_file))
		if (nrow(cur_file)>0)  {
			cur_corrected_TCGA_file <- grep(basename(tissue_file),all_corrected_TCGA_files,value=T)
			cur_corrected_TCGA_file <- colMeans(data.frame(fread(cur_corrected_TCGA_file)),na.rm=T)
			cur_corrected_TCGA_file <- data.frame(pure_cpg =cur_corrected_TCGA_file) 
###TODO

			# cur_corrected_TCGA_file <- data.frame(fread(grep(sub(".csv",basename(tissue_file),all_corrected_TCGA_files,value=T)))
			
			cur_file <- merge(cur_file,cur_corrected_TCGA_file,by.x="CpGs_common",by.y="row.names")
			cur_file$mean_tumor <- cur_file$pure_cpg
		}
			scores_of_meth <- data_on_files(cur_file,sample_mapped,sub("\\..*","",basename(sample_bed)),tissue_file)
		data.frame(cancer_type=sub("\\.csv","",basename(tissue_file)),sample=sub("\\.hg38.*","",basename(sample_bed)),scores_of_meth)
	},mc.cores=future::availableCores()[[1]]))
}))

####====scoring
# samples_healthy <- c("BC02","BC03","BC05","BC04","BC12","HU005.10","HU005.11","HU005.12")
samples_healthy <- unlist(strsplit(args[3],split=",")) 

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
# file_to_write <- "methylation_in_sample_LUAD_score_table_all_samples_new_healthy_new_score_tumor_only_unique_to_LUAD.csv" 
# file_to_write <- "methylation_in_sample_LUAD_score_table_all_samples_remora1_060222_new_score_tumor_only_unique_to_LUAD.csv"  
# file_to_write <- "methylation_in_sample_LUAD_score_table_all_samples_remora1_080222_new_score_tumor_only_unique_to_LUAD.csv"  
file_to_write <- paste0("methylation_in_sample_LUAD_new_score_tumor_only_unique_to_LUAD",args[2],".csv" )
# file_to_write <- "methylation_in_sample_per_tissue_table_all_samples_new_score_tumor_only.csv"

file_to_write_path <- file.path(output_dir,"atlas_vs_healthy",file_to_write)
fwrite(all_table_methylation_in_tissue,file_to_write_path)
file_name_to_link <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization.csv"
system(paste("ln -sf",file.path(output_dir,"atlas_vs_healthy",file_to_write),file.path(output_dir,"atlas_vs_healthy",file_name_to_link)))
message(paste("The file",file_to_write_path,"was written"))
