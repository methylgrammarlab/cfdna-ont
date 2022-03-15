source("/sci/labs/bermanb/ekushele/icore-home/script/r/filippo/functions_atlas.R")
args = commandArgs(trailingOnly=TRUE)

cores=future::availableCores()[[1]]


output_dir="/sci/labs/bermanb/ekushele/icore-data/filippo/data/new_atlas/results/"

num_selected_cpgs=as.numeric(args[1])*2
best_hypo_hyper_probes_file <- paste0( "/sci/labs/bermanb/ekushele/icore-data/filippo/data/josh_moss_probes/creating_new_atlas_scripts/CpGs.100bp-block.",num_selected_cpgs/2,".X1.xls")
best_hypo_hyper_probes <- data.frame(fread(best_hypo_hyper_probes_file))
best_hypo_hyper_probes$name <- gsub(" ","_",best_hypo_hyper_probes$name)
best_hypo_hyper_probes$name <- sub("-","\\.",best_hypo_hyper_probes$name)




#####=========score of methylation
# data_on_files <- function(cur_tissue_file,sample_mapped,cell_type,sample_name) {
# 	CpGs <- unique(cur_tissue_file$acc)
# 	#methylation data
# 	sample_overlap_meth <- sample_mapped[which(sample_mapped$cpg_index %in% sub("cg","",CpGs)),c("cpg_index","methylation")]
# 	# sample_overlap_meth <- sample_mapped[intersect(sample_mapped$cpg_index, sub("cg","",CpGs)),c("cpg_index","methylation")]
# 	sample_overlap_meth$cpg_index <- sub("^","cg",sample_overlap_meth$cpg_index)
	
# 	#include only CpGs that are in the nanopore methylation
# 	CpGs <- sample_overlap_meth$cpg_index
# 	sample_overlap_meth <- sample_overlap_meth[match(CpGs,sample_overlap_meth$cpg_index),]
  
# 	#healthy data
# 	cur_all_healthy_mean <- all_healthy_mean[CpGs]
# 	cur_all_healthy_mean <- as.numeric(cur_all_healthy_mean[match(CpGs,names(cur_all_healthy_mean))])
# 	#cell-type data
# 	cur_all_data_meth_mean <- all_cell_type_mean[[cell_type]]
# 	cur_all_data_meth_mean <- cur_all_data_meth_mean[CpGs]
# 	cur_all_data_meth_mean <- as.numeric(cur_all_data_meth_mean[match(CpGs,names(cur_all_data_meth_mean))])
	
# 	index_not_na <- which(!is.na(cur_all_healthy_mean) & !is.na(cur_all_data_meth_mean))

# 	cur_all_healthy_mean <- cur_all_healthy_mean[index_not_na]
# 	cur_all_data_meth_mean <- cur_all_data_meth_mean[index_not_na]
# 	sample_overlap_meth <- sample_overlap_meth[index_not_na,]
	
# 	print(length(cur_all_healthy_mean))
# 	score_of_cpg <- (sample_overlap_meth$methylation - cur_all_healthy_mean)/(cur_all_data_meth_mean-cur_all_healthy_mean)
# 	# score_of_cpg <- (sample_overlap_meth$methylation - cur_all_data_meth_mean)/-(cur_all_data_meth_mean-cur_all_healthy_mean)
# 	nnls_score_of_cpg <- nnls::nnls(as.matrix(data.frame(cur_all_data_meth_mean,cur_all_healthy_mean)),sample_overlap_meth$methylation)$x
# 	if (cell_type=="Lung_cells")
# 		fwrite(data.frame(mean_Lung=cur_all_data_meth_mean,mean_healthy=cur_all_healthy_mean,mehtylation=sample_overlap_meth$methylation),file.path(output_dir,"atlas_joshCpGs_vs_healthy","check",paste0(sample_name,"_",num_selected_cpgs,"_check.csv")))
# 	data.frame(mean_data_hypo_hyper=mean(score_of_cpg),nnls_hypo_hyper=nnls_score_of_cpg[[1]]/sum(nnls_score_of_cpg))

# }
pvalue=1 #as.numeric(args[2])# 0.01
print(pvalue)
diff_mean=0 #as.numeric(args[3]) #0.3
print(diff_mean)
data_on_files <- function(cur_file_atlas,sample_mapped,cell_type,sample_name) {
	cur_file <- data.frame(fread(file.path(output_dir,"atlas_joshCpGs_vs_healthy","t_test",paste0(cell_type,".csv"))))
	
	#this is to remove any probes with only zero or one values
	cur_file <-cur_file[cur_file$FDR_adjust < pvalue & abs(cur_file$mean_diff) > diff_mean & !is.na(cur_file$FDR_adjust), ]
	
	CpGs <- unique(cur_file_atlas$acc) 
	sample_mapped$cpg_index <- sub("^","cg",sample_mapped$cpg_index)

	sample_overlap_meth_coordinates <- sample_mapped[which(sample_mapped$cpg_index %in% intersect(sample_mapped$cpg_index,CpGs)),c("cpg_index","V5","V6","V7", "methylation")]
	sample_overlap_meth <- aggregate(methylation~cpg_index,data=sample_overlap_meth_coordinates,FUN=mean, na.rm=TRUE, na.action=na.pass)
	CpGs <- sample_overlap_meth$cpg_index

	# CpGs_numbers_rep <- plyr::count(CpGs)
	# CpGs_numbers_rep <- CpGs_numbers_rep[match(unique(CpGs),CpGs_numbers_rep$x),] 
	cur_file <- cur_file[which(cur_file$CpGs_common %in% intersect(CpGs,cur_file$CpGs_common)),]
	cur_file <- cur_file[match(unique(CpGs),cur_file$CpGs_common),]
	# cur_file <- cur_file[rep(1:nrow(cur_file), times=CpGs_numbers_rep$freq),]
	
	#remove NA data, which is probably because met_healthy had only one or zero value so t-test on that CpG couldn't be preformed
	index_not_na <- which(!is.na(cur_file$CpGs_common))
	sample_overlap_meth <- sample_overlap_meth[index_not_na,]
	cur_file <- cur_file[index_not_na,]
	sample_overlap_meth <- merge(sample_overlap_meth,aggregate(.~cpg_index,sample_overlap_meth_coordinates[,-c(5)],head,1),by="cpg_index")
	score_of_cpg <- (sample_overlap_meth$methylation - cur_file$mean_healthy)/(cur_file$mean_tumor-cur_file$mean_healthy)
	nnls_score_of_cpg <- nnls::nnls(as.matrix(data.frame(cur_file$mean_tumor,cur_file$mean_healthy)),sample_overlap_meth$methylation)$x
	
	#write file for ben 
	# if (grepl("Lung_cells",cell_type)) 
		fwrite(data.frame(chr=sample_overlap_meth$V5,start=sample_overlap_meth$V6,end=sample_overlap_meth$V7,mean_cell=cur_file$mean_tumor,mean_healthy=cur_file$mean_healthy,methylation=sample_overlap_meth$methylation,p_value=cur_file$p_val,FDR=cur_file$FDR_adjust),file.path(output_dir,"atlas_joshCpGs_vs_healthy",paste0("check_FDR",pvalue,"_diffmean",0),paste0(cell_type,"_",sample_name,"_",num_selected_cpgs,"_check.csv")))
	data.frame(mean_data_hypo_hyper=mean(score_of_cpg),nnls_hypo_hyper=nnls_score_of_cpg[[1]]/sum(nnls_score_of_cpg),length_overlapping_meth=length(score_of_cpg))
	
}

#find the hypo-hyper methylation for each sample
samples_mapped_to_CpG <- dir(args[2],full.names=T,pattern="_shifted_mapped_CpG.bed.gz")
# samples_mapped_to_CpG <- dir(file.path("/sci/labs/bermanb/ekushele/icore-data/filippo/data/samples_nanopore_020222/","samples_merged_CpG_450K_links_060222"),full.names=T,pattern="_shifted_mapped_CpG.bed.gz")
# samples_mapped_to_CpG <- dir(file.path("/sci/labs/bermanb/ekushele/icore-data/filippo/data/samples_nanopore/","samples_merged_CpG_450K"),full.names=T,pattern="_shifted_mapped_CpG.bed.gz")
#samples_mapped_to_CpG <- dir(file.path("/sci/labs/bermanb/ekushele/icore-data/filippo/data/samples_nanopore/","samples_merged_CpG_450K"),full.names=T,pattern=".hg38_shifted_check.bed")
#samples_mapped_to_CpG <- dir(file.path("/sci/labs/bermanb/ekushele/icore-data/filippo/data/samples_nanopore/","samples_merged_CpG_450K"),full.names=T,pattern="_unique_for_new_file.bed")


all_table_methylation_in_tissue <- do.call(rbind,lapply(samples_mapped_to_CpG,function(sample_bed) {
	print(sample_bed)
	sample_mapped <- data.frame(fread(sample_bed,colClasses = c('V4'='character')))
	colnames(sample_mapped)[4] <- "cpg_index"
	colnames(sample_mapped)[8] <- "methylation"
	#take only sample_mapped where there are CpGs for specific cell types
	list_of_names <- unique(best_hypo_hyper_probes$name)
	do.call(rbind,lapply(list_of_names,function(cell_type) {
		print(cell_type)
		cur_tissue_file <- subset(best_hypo_hyper_probes,name==cell_type)
		scores_of_meth <- data_on_files(cur_tissue_file,sample_mapped,cell_type,sub("\\..*","",basename(sample_bed)))

		data.frame(cancer_type=gsub(" ","_",cell_type),sample=basename(sample_bed),scores_of_meth)
	}))
}))


#scoring
#samples_healthy <- c("BC02","BC03","BC05","BC04","BC12","HU005.10","HU005.11","HU005.12")
samples_healthy <- unlist(strsplit(args[4],split=","))
all_table_methylation_in_tissue[is.na(all_table_methylation_in_tissue)] <- 0
all_table_methylation_in_tissue$category <- "Tumor"
all_table_methylation_in_tissue$category[grep(paste0(samples_healthy,collapse="|"),all_table_methylation_in_tissue$sample)] <- "Healthy"
all_table_methylation_in_tissue <- all_table_methylation_in_tissue[order(all_table_methylation_in_tissue$cancer_type),]


#file_to_write <- paste0("methylation_new_score_check_num_cpgs_",num_selected_cpgs,"_new_mapping.csv")
# file_to_write <- paste0("methylation_new_score_check_num_cpgs_",num_selected_cpgs,"_old_mapping.csv")
# file_to_write <- paste0("methylation_new_score_check_num_cpgs_",num_selected_cpgs,"_unique_to_new_mapping.csv")
# file_to_write_path <- file.path(output_dir,"atlas_joshCpGs_vs_healthy",file_to_write)
# fwrite(all_table_methylation_in_tissue,file_to_write_path)
# message(paste("The file",file_to_write_path),"was written")


# types_of_raw_scores <- unique(sub("normalized.","",grep("hypo_hyper",colnames(all_table_methylation_in_tissue),value=T)))
# normalized_data <- do.call(cbind,lapply(types_of_raw_scores,function(x) {
	# medians <- aggregate(get(x)~category+cancer_type,all_table_methylation_in_tissue,median)[c(T,F),]
	# all_table_methylation_in_tissue[,x]- rep(medians[,3],each=length(unique(all_table_methylation_in_tissue$sample)))
# }))
# colnames(normalized_data) <- paste0("normalized.",types_of_raw_scores)
# all_table_methylation_in_tissue <- cbind(all_table_methylation_in_tissue[,grep("normalized",colnames(all_table_methylation_in_tissue),invert=T)],normalized_data)

all_table_methylation_in_tissue$sample <- sub("\\.hg38.*","",all_table_methylation_in_tissue$sample)

######for Ben to get correlation with RMSE
# tumor_fractions <- data.frame(read.table(text=system("less -S /sci/labs/bermanb/ekushele/icore-data/filippo/data/tumor_fractions/tumorFracTable_subDT_percCNA_015_maxCN3setto0_with_num_reads.txt| awk -F '\t' -v OFS='\t' '{print $1,$8,$20}'",intern=T),header=T))
# tumor_fractions$sample <- sub("\\.fragmentomic.hg38","",tumor_fractions$sample)
# all_table_methylation_in_tissue <- merge(all_table_methylation_in_tissue, tumor_fractions,by="sample")


# num_samples=length(unique(all_table_methylation_in_tissue$sample))
# all_table_methylation_in_tissue$RMSE_mean_data_hypo_hyper <- sqrt(((all_table_methylation_in_tissue$tumorFrac - all_table_methylation_in_tissue$mean_data_hypo_hyper)^2)/num_samples)
# all_table_methylation_in_tissue$RMSE_nnls_hypo_hyper <- sqrt(((all_table_methylation_in_tissue$tumorFrac - all_table_methylation_in_tissue$nnls_hypo_hyper)^2)/num_samples)
# all_table_methylation_in_tissue$tumorFrac <-c()
# all_table_methylation_in_tissue$num_reads <-c()
# file_to_write <- "methylation_new_score_2000_FDR1_meandiff0_with_RMSE.csv"
# file_to_write <- paste0("methylation_new_score_",args[1],"_FDR",pvalue,"_meandiff",diff_mean,".csv")
# file_to_write <- paste0("methylation_new_score_all_lung_cells",args[1],"_FDR",pvalue,"_meandiff",diff_mean,".csv")
# file_to_write <- paste0("methylation_new_score_all_lung_cells_new_healthy",args[1],"_FDR",pvalue,"_meandiff",diff_mean,".csv")
#file_to_write <- paste0("methylation_new_score_all_lung_cells_new_healthy_020222",args[1],"_FDR",pvalue,"_meandiff",diff_mean,".csv")
# file_to_write <- paste0("methylation_new_score_all_lung_cells_remora1_060222",args[1],"_FDR",pvalue,"_meandiff",diff_mean,".csv")
#file_to_write <- paste0("methylation_new_score_all_lung_cells_remora1_080222",args[1],"_FDR",pvalue,"_meandiff",diff_mean,".csv")

output_dir="/sci/labs/bermanb/ekushele/icore-data/filippo/data/new_atlas/results/"
file_name_to_link <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization.csv"
#all_table_methylation_in_tissue <- data.frame(fread(file.path(output_dir,"atlas_joshCpGs_vs_healthy",file_name_to_write)))
# file_to_write <- paste0("methylation_new_score_all_lung_cells",args[3]	"_",args[1],"_FDR",pvalue,"_meandiff",diff_mean,".csv")
file_to_write <- paste0("methylation_in_sample_new_score_tumor",args[3],".csv" )

file_to_write_path <- file.path(output_dir,"atlas_joshCpGs_vs_healthy",file_to_write)
fwrite(all_table_methylation_in_tissue,file_to_write_path)
message(paste("The file",file_to_write_path,"was written"))
system(paste("ln -sf",file.path(output_dir,"atlas_joshCpGs_vs_healthy",file_to_write),file.path(output_dir,"atlas_joshCpGs_vs_healthy",file_name_to_link)))
