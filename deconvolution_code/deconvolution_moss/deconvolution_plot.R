source("/sci/labs/bermanb/ekushele/icore-home/script/r/filippo/functions_atlas.R")
library(argparse)

parser <- ArgumentParser(formatter_class= 'argparse.RawTextHelpFormatter')

parser$add_argument("-i","--input_deconv", action="store", help="Input deconvulotion results from deconvolate.py script, for plot")#,required=T)
parser$add_argument("-o","--output_dir", action="store", help="Output dir path. Default: ../../input_matrix")
parser$add_argument("--suffix", action="store", help="suffix file name,start with '_'",default="")
parser$add_argument("-g","--group_file", action="store", help="CSV file with 'source' group name and 'destination' group name. All groups not shown will be grouped as 'others'. Optional:group file with 'color' column (3rd column) that will have the color for each group")
parser$add_argument("-s","--sample_order", action="store", help="sample list in the needed order, devided by commas (without whitespace)")
parser$add_argument("-t","--trim_string", action="store", help="string to trim from sample name")
parser$add_argument("-t2","--trim_string2", action="store", help="string to trim from sample name  #2, usually trim from end")
parser$add_argument("-r","--replace", action="store", help="string to replace after trimming",default="")
parser$add_argument("-ht","--hiding_threshold", action="store", help="Threshold to show results in the 'Other' category. Default: 0.01",default=0.01)
parser$add_argument("-c","--category", action="store", help="category for facet wrap, in the same order of sample order")

### #A8A7A5
args <- parser$parse_args()  
suffix_file_name=".pdf"

##debug
# args$input_deconv="/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/deconvolusion/CpGs.100bp-block.2000/all_samples_for_deconvolusion_CpGs.100bp-block.2000_deconv_output.csv"
# output_dir_path="/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/deconvolusion/CpGs.100bp-block.2000"
# suffix="2000"
# args$sample_order="BC02,BC03,BC04,BC05,BC09,BC08,BC01,19_326,BC10,BC11"
# args$group_file<-"/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/deconvolusion/group_file_for_plot.csv"
# args$category <- "Healthy,Healthy,Healthy,Healthy,Cancer,Cancer,Cancer,Cancer,Cancer,Cancer"

#orgenize input data 

deconv_result <- read.csv(args$input_deconv,check.names=F) #use read.csv because of leading x in column name
names(deconv_result)[1] <- "cell_type"
if (!is.null(args$trim_string)) {
	colnames(deconv_result)[-c(1)] <- sub(args$trim_string,args$replace,colnames(deconv_result)[-c(1)])
}
if (!is.null(args$trim_string2)) {
        colnames(deconv_result)[-c(1)] <- sub(args$trim_string2,"",colnames(deconv_result)[-c(1)])
}

if (!is.null(args$sample_order)) {
	sample_order <- unlist(strsplit(args$sample_order,split=","))
	print(sample_order)
	print(colnames(deconv_result))
	deconv_result[,-c(1)] <- deconv_result[,sample_order]
	 colnames(deconv_result)[-c(1)] <-  sample_order
}

indexes_colors=c(2,1,3:8)

colors_for_plot <- NA
if (!is.null(args$group_file)) {
	group_data <- data.frame(fread(args$group_file))
	# cell_types_levels <- rev(unique(group_data[,2]))
	cell_types_levels <- c("Others",rev(c("Monocytes","Lymphocytes","Granulocytes","Megakaryocytes","Hepatocytes","Vascular endothelial cells","Epithelial cells")))
	group_data[group_data[,2]=="",3] <- "#A7A7A7"
	group_data[group_data[,2]=="",2] <- "Others" #group_data[group_data[,2]=="",1]
	deconv_result$cell_type <- group_data[,2]
	# cell_types_levels <- c("Others",cell_types_levels)
	deconv_result$cell_type <- factor(deconv_result$cell_type,levels=cell_types_levels[indexes_colors])
	colors_for_plot <- setNames(unique(group_data[,3]), unique(group_data[,2]))
	colors_for_plot <- colors_for_plot[cell_types_levels]
}
deconv_result <- reshape2::melt(deconv_result, id.vars = "cell_type")

if (!is.null(args$sample_order)) {
	sample_order <- unlist(strsplit(args$sample_order,split=","))
	deconv_result$variable <- factor(deconv_result$variable,levels=sample_order)
}


if (!is.null(args$category)) {
	category_for_plot <- unlist(strsplit(args$category,split=",")) 
	category_for_plot <- factor(category_for_plot,levels=unique(category_for_plot))
	#category_for_plot <- category_for_plot[levels(deconv_result$variable)]
	deconv_result$category <- rep(category_for_plot,each=length(which(deconv_result$variable==deconv_result$variable[1])))
}

library(grid)
colors_for_plot <-colors_for_plot[indexes_colors]

sizeit <- function(p, panel.size = 2, default.ar=1){

	gb <- ggplot_build(p)
	# first check if theme sets an aspect ratio
	ar <- gb$plot$coordinates$ratio

	# second possibility: aspect ratio is set by the coordinates, which results in 
	# the use of 'null' units for the gtable layout. let's find out
	g <- ggplot_gtable(gb)
	nullw <- sapply(g$widths, attr, "unit")
	nullh <- sapply(g$heights, attr, "unit")

	# ugly hack to extract the aspect ratio from these weird units
	if(any(nullw == "null"))
	ar <- unlist(g$widths[nullw == "null"]) / unlist(g$heights[nullh == "null"])

	if(is.null(ar)) # if the aspect ratio wasn't specified by the plot
	   ar <- default.ar

	# ensure that panel.size is always the larger dimension
	if(ar <= 1 ) panel.size <- panel.size / ar

	g$fullwidth <- convertWidth(sum(g$widths), "in", valueOnly=TRUE) + 
	panel.size
	g$fullheight <- convertHeight(sum(g$heights), "in", valueOnly=TRUE) + 
	panel.size / ar
	print(paste("ranges","width=",g$fullwidth, "height=",g$fullheight))

	class(g) <- c("sizedgrob", class(g))
	g
}


print.sizedgrob <- function(x){
	print(paste("ranges","width=",x$fullwidth, "height=",x$fullheight))
	dev.new(width=x$fullwidth, height=x$fullheight)
	grid.draw(x)
}

plot_barplot <- function(deconv_result) {
	font_size=21

	p <- ggplot(deconv_result,aes(x=variable,y=value,fill=cell_type))+
		geom_bar(stat = "identity", position = "fill")+
		ylab("Tissue contribution (Normalized)")+
		xlab("Samples")+
		guides(fill = guide_legend(ncol = 1))+#override.aes = list(fill = colors_for_plot))+
		ggtitle(paste("Deconvolution using",args$suffix,"CpGs"))+
		theme_classic()+ 
			 theme(legend.title = element_blank(),
					axis.text.x = element_text(face = "bold", size = font_size,angle=90,hjust=1,vjust=0.5), 
					axis.text = element_text(face = "bold", size = font_size),
					# axis.title = element_text(face = "bold", size = font_size),
					axis.title = element_text(face = "bold", size = font_size),
					legend.text=element_text(size=15,face="bold"),
					panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
					#		 panel.border = element_rect(colour = "black", fill=NA, size=2),
	if (!is.na(colors_for_plot[1])) 
		p <- p + scale_fill_manual(labels=names(colors_for_plot),values=colors_for_plot)
	if (!is.null(args$category))
		p <- p + facet_grid(.~category,scales="free_x",space = "free_x")

	return(p)
}

deconv_result<- aggregate(value~.,deconv_result,sum)
p <- plot_barplot(deconv_result)
num_samples <- length(unique(deconv_result$value))


# p <- sizeit(p, num_samples*4/9) 
# print(p$fullheight)
# print(p$fullwidth)
output_dir_path <- ifelse(is.null(args$output_dir),file.path(dirname(args$input_deconv)),args$output_dir)
create_dir(output_dir_path)  
file_name <- file.path(output_dir_path,paste0("deconvolution_plot",args$suffix,suffix_file_name))
ggsave(file_name,p,width=10,height=7.5)
message(paste("File was saved:",file_name))
