library(data.table)
library(argparse)

parser <- ArgumentParser(formatter_class= 'argparse.RawTextHelpFormatter')

parser$add_argument("-i","--inputdir", action="store", help="Input directory of files to be merged")#,required=T)
parser$add_argument("-o","--output_dir", action="store", help="Output dir path. Default: input dir")
parser$add_argument("-p","--pattern", action="store", help="pattern of files to be searched for")
parser$add_argument("-s","--suffix", action="store", help="suffix file name,start with '_'",default="")
parser$add_argument("-g","--genome",action="store",help="genome assembly for list of cpgs from 450k",default="hg38")
parser$add_argument("-e", "--empty_samples", action="store_false", help="Add empty samples with 'NA's to table")

args <- parser$parse_args()  
print(args$genome)
print(paste0("/sci/labs/bermanb/bermanb/icore-data/projects/genomic-data-misc/Infinium/illumina-methyl-450k-manifest.cgs.0based.",args$genome,".bed"))

list_of_cpgs <- read.csv(paste0("/sci/labs/bermanb/bermanb/icore-data/projects/genomic-data-misc/Infinium/illumina-methyl-450k-manifest.cgs.0based.",args$genome,".bed"), header=F,sep="\t")[,4]
if (!grepl("cg",list_of_cpgs[1])) #column doesn't starts with cg
	list_of_cpgs <- sub("^","cg",list_of_cpgs)
# args$inputdir<-"/sci/labs/bermanb/ekushele/icore-data/filippo/data/samples_nanopore_020222/samples_merged_CpG_450K_links"
# args$pattern<-".hg38_shifted_mapped_CpG.bed.gz" 
# args$suffix<-"_020226_check"

files_to_loop=dir(args$inputdir, pattern=args$pattern,full.names=T)
print(files_to_loop)

files_to_loop_not_empty<-lapply(files_to_loop,function(x) file.size(x)>28)
files_to_loop_empty <-lapply(files_to_loop,function(x) file.size(x)==28)

all_data_frames <- lapply(files_to_loop[unlist(files_to_loop_not_empty)], function(x) { 
	# cur_sample <- fread(x)[,c(4,5,6,7,8)]
	cur_sample <-fread(x,colClasses = c('V4'='character'))[,c(4,5,6,7,8)]
	colnames(cur_sample) <- c("CpGs", "chr","start","end",sub(args$pattern,"",basename(x)))
	if (!grepl("cg",cur_sample$CpGs[1])) #column doesn't starts with cg
		cur_sample$CpGs <- sub("^","cg",cur_sample$CpGs)	
	cur_sample
})

output_dir_path <- ifelse(is.null(args$output_dir),args$inputdir,args$output_dir)
final_data_frame <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by=c("CpGs","chr","start","end"),allow.cartesian=TRUE, all = TRUE),all_data_frames) 
final_data_frame <- aggregate(.~CpGs,data=final_data_frame[,-c(2,3,4)],FUN=mean, na.rm=TRUE, na.action=na.pass)

#add empty columns
empty_column <- unlist(lapply(files_to_loop[unlist(files_to_loop_empty)],function(x) {
		sub(args$pattern,"",basename(x))
}))

#manage CpGs to contain exactly the number of CpGs from "list of cpgs")
all_cpgs_not_in_data_frame <- data.frame(CpGs=setdiff(list_of_cpgs,final_data_frame$CpGs))
final_data_frame <- merge(final_data_frame,data.frame(all_cpgs_not_in_data_frame),by="CpGs",all=T)
if (isTRUE(args$empty_column))
	final_data_frame[,empty_column] <- NA
dim(final_data_frame)

file_name=file.path(output_dir_path,paste0("all_samples_for_deconvolusion",args$suffix,".csv"))
fwrite(final_data_frame, file_name)
message(paste0("File name: ",file_name," was created"))
	
