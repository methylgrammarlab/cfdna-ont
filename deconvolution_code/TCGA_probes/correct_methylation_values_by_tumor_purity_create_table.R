library(ELMER)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(data.table)
source("/cs/icore/ekushele/script/r/filippo/functions_atlas.R")

cores=future::availableCores()[[1]]

#########Get all the data: TCGA, normal cfDNA 
files_TCGA <- dir("/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/new_atlas/DNA_met_sesame/",full.names = TRUE,pattern = "rda")

TCGA_master_call_table <- data.frame(fread(file.path("/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data","tumor_fractions","TCGA_mastercalls.abs_tables_JSedit.fixed.txt")))
TCGA_master_call_table <- subset(TCGA_master_call_table,solution=="new")


all_white_cells_methylation <- data.frame(fread("cell_type file
B.cells_EPIC /cs/cbio/josh/data/array/deconvolution/reference/samples/B_cells_EPIC.csv.gz
CD4T.cells_EPIC /cs/cbio/josh/data/array/deconvolution/reference/samples/CD4T_cells_EPIC.csv.gz
CD8T.cells_EPIC /cs/cbio/josh/data/array/deconvolution/reference/samples/CD8T_cells_EPIC.csv.gz
Monocytes_EPIC /cs/cbio/josh/data/array/deconvolution/reference/samples/Monocytes_EPIC.csv.gz
Neutrophils_EPIC /cs/cbio/josh/data/array/deconvolution/reference/samples/Neutrophils_EPIC.csv.gz
NK.cells_EPIC /cs/cbio/josh/data/array/deconvolution/reference/samples/NK_cells_EPIC.csv.gz
"))

# rbind.fill.DT <- function(ll) {
    # # changed sapply to lapply to return a list always
    # all.names <- lapply(ll, names)
    # unq.names <- unique(unlist(all.names))
    # ll.m <- rbindlist(lapply(seq_along(ll), function(x) {
        # tt <- ll[[x]]
        # setattr(tt, 'class', c('data.table', 'data.frame'))
        # data.table:::settruelength(tt, 0L)
        # invisible(alloc.col(tt))
        # tt[, c(unq.names[!unq.names %chin% all.names[[x]]]) := NA_character_]
        # setcolorder(tt, unq.names)
    # }))
# }

methylation_white_blood_cell <- lapply(all_white_cells_methylation$file,function(cell_type) {
	met_tissue <- data.frame(fread(cell_type))
	rownames(met_tissue) <- met_tissue$V1
	met_tissue$V1 <- c()
	#get cpg values, remove rows with all na
	met_tissue <- met_tissue[rowSums(is.na(met_tissue)) != ncol(met_tissue), ] # result
	colnames(met_tissue) <- c()	
	mean_cpg <- data.frame(sample_name=rowMeans(met_tissue,na.rm=T))
	colnames(mean_cpg)[1] <- sub("_EPIC.*","",basename(cell_type))
	mean_cpg
})

for (i in seq(1:length(methylation_white_blood_cell))) {
	for (j in i+1:length(methylation_white_blood_cell)-i) {
		# print(paste(i,j))
		# print(setdiff(rownames(methylation_white_blood_cell[[i]]),rownames(methylation_white_blood_cell[[j]])))
		# print(setdiff(rownames(methylation_white_blood_cell[[j]]),rownames(methylation_white_blood_cell[[i]])))
		
		setdiff2 <- setdiff(rownames(methylation_white_blood_cell[[i]]),rownames(methylation_white_blood_cell[[j]]))
		if (length(setdiff2) >0 ) {
			dataframe2 <- data.frame(cpg=setdiff2,value=NA)
			rownames(dataframe2) <- dataframe2$cpg
			dataframe2$cpg <- c()
			colnames(dataframe2) <- colnames( methylation_white_blood_cell[[j]])
			methylation_white_blood_cell[[j]] <- rbind(methylation_white_blood_cell[[j]],dataframe2)
		}
		setdiff1 <- setdiff(rownames(methylation_white_blood_cell[[j]]),rownames(methylation_white_blood_cell[[i]]))

		if (length(setdiff1) > 0) {
			dataframe1 <- data.frame(cpg=setdiff1,value=NA)
			rownames(dataframe1) <- dataframe1$cpg
			dataframe1$cpg <- c()
			colnames(dataframe1) <- colnames( methylation_white_blood_cell[[i]])
			methylation_white_blood_cell[[i]] <- rbind(methylation_white_blood_cell[[i]],dataframe1)
			}
		}
}
for (i in 2:length(methylation_white_blood_cell)) {
	methylation_white_blood_cell[[i]][,1] <- methylation_white_blood_cell[[i]][order(match(rownames(methylation_white_blood_cell[[1]]),rownames(methylation_white_blood_cell[[i]]))),1]
	rownames(methylation_white_blood_cell[[i]]) <- rownames(methylation_white_blood_cell[[1]])
}
methylation_white_blood_cell <- do.call(cbind,methylation_white_blood_cell)

average_beta_white_blood_cell <- rowMeans(methylation_white_blood_cell,na.rm=T)

parallel::mclapply(files_TCGA, function(tcga_tissue) {
		filename <- sub(".rda","",basename(tcga_tissue))
		filename_to_write <- file.path("/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/new_atlas/TCGA_purity/",paste0(filename,".csv"))
		print(filename_to_write)
		if (!file.exists(filename_to_write)) {
			data <- get(load(tcga_tissue))
			print(paste("@@@@@@@@@@@@@@@@@",filename))
			met <- ELMER:::makeSummarizedExperimentFromDNAMethylation(
			data, genome = "hg38", met.platform = "450K")
			# print(all(met$samples == rownames(TCGAbiolinks:::colDataPrepareTCGA(met$samples))))
			colData(met) <- TCGAbiolinks:::colDataPrepareTCGA(met$samples)
			# print(met)

			#remove un-needed groups: take only solid tumor
			group = "Primary solid Tumor"
			met <- met[,colData(met)[, "definition", drop = TRUE] %in% c(group)]

			#get cpg values, remove rows with all na
			met_tissue <- assay(met)
			met_tissue <- met_tissue[rowSums(is.na(met_tissue)) != ncol(met_tissue), ] # result
			print(dim(met_tissue))
			if (nrow(met_tissue) >0) {
				met_tissue <- met_tissue[intersect(names(average_beta_white_blood_cell),rownames(met_tissue)),]
				cur_average_beta_white_blood_cell <- average_beta_white_blood_cell[intersect(rownames(met_tissue),names(average_beta_white_blood_cell))]
				# met_tissue <- bind(met_tissue,average_beta_white_blood_cell)

				# average_meth_tissue <-  data.frame(t(colMeans(met_tissue,na.rm=T)))
				met_tissue <- data.frame(t(met_tissue))
				met_tissue$sample <- substr(rownames(met_tissue) ,0,15)
				met_tissue <- merge(met_tissue,TCGA_master_call_table[,c("array","purity")],by.x="sample",by.y="array")
				cur_average_beta_white_blood_cell <- c(sample=NA,cur_average_beta_white_blood_cell,purity=NA)
				# met_tissue <- data.frame(t(met_tissue))
				# colnames(met_tissue) <- met_tissue[1,]
				# met_tissue <- met_tissue[-c(1),]
				met_tissue <- rbind(met_tissue,cur_average_beta_white_blood_cell)

				cg_columns <- grep("^cg",colnames(met_tissue))
				index_average_WBC <- nrow(met_tissue)
				# met_tissue[-c(index_average_WBC),c(cg_columns)] <- 
				pure_methylation <- simplify2array(lapply(colnames(met_tissue[c(cg_columns)]),function(x)  {

					average_WBC <- met_tissue[[x]][index_average_WBC]
					methylation_data <- met_tissue[[x]][-c(index_average_WBC)]
				#	(met_tissue[[x]]-(1-met_tissue$purity)*cur_average_beta_white_blood_cell[[x]])/met_tissue$purity
					(methylation_data - (1-met_tissue$purity[-c(index_average_WBC)])*average_WBC)/met_tissue$purity[-c(index_average_WBC)]
				})) 	  
				colnames(pure_methylation) <- colnames(met_tissue[c(cg_columns)])
				print(dim(pure_methylation))
				fwrite(pure_methylation,filename_to_write)
			}
			
		}
		else {
			print(paste("The file:",filename_to_write,"already exiss"))   
		}
})



