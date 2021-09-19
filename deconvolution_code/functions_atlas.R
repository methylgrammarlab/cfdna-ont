suppressPackageStartupMessages(library(ComplexHeatmap))

library(data.table)
library(ComplexHeatmap)
library(patchwork)

plot_methylation <- function(all_meth_data,names_annotate_rows,names_annotate_type,p_val_annotate,mean_meth_diff,num_cpg=500,title_all="",cluster_rows=F,remove_TGCT=F,order_by_mean=T,order_by_pval=F) {
	num_cpg_to_plot=ifelse(num_cpg>ncol(all_meth_data),ncol(all_meth_data),num_cpg)

	all_meth_data <- data.frame(all_meth_data)[,1:num_cpg_to_plot]
	
	column_order <- if(isTRUE(order_by_pval)) p_val_annotate else abs(mean_meth_diff)
	
	mean_annotate_order_rows <- ave(rowMeans(all_meth_data,na.rm=T), names_annotate_rows, FUN = mean)
	order_rows <- if (isTRUE(order_by_mean)) mean_annotate_order_rows else names_annotate_rows

	####plot
	ht <-draw(Heatmap(all_meth_data,
		show_column_dend=F,show_column_names =F, show_row_names=F,show_row_dend=F, use_raster=T,name="methylation",column_title=title_all,
		left_annotation = rowAnnotation(Tissue = names_annotate_rows,Type=names_annotate_type,show_legend = c(T,F),,annotation_name_side="top",
			col=list(Type=c("LUAD"="darkorchid","Non_LUAD"="plum2"))),
		bottom_annotation = HeatmapAnnotation( "-log10 p_val" = anno_points(-log10(p_val_annotate))  ,
	   			"mean_meth_diff"= anno_points(mean_meth_diff)),
	    cluster_rows =cluster_rows,row_order=order(order_rows),
	    cluster_columns=F,column_order=order(column_order,decreasing=F),
		annotation_legend_list=list(#Legend(title="Source",labels=c("TCGA","Non TCGA"),legend_gp=gpar(fill=c("royalblue","seagreen1"))),
		Legend(title="Type",labels=c("LUAD","Non LUAD"),legend_gp=gpar(fill=c("darkorchid","plum2")))),
	merge_legend = TRUE
	))
	return(ht)

}

#===============meth data according to names cutoff 
get_meth_data <- function(all_csv_files,names_cutoff_cpg) {
	meth_data <- parallel::mclapply(all_csv_files,function(name) {
		base_name=basename(name)
		#name=all_csv_files[1]
		full_meth_data <- data.frame(fread(name)[,-c(1)])
		full_meth_data <- t(full_meth_data[which(full_meth_data$CpG %in% names_cutoff_cpg),])
		indexes <- grep("values",row.names(full_meth_data))
		row.names(full_meth_data)[indexes] <- paste0(full_meth_data["tissue",1],"_",row.names(full_meth_data)[indexes])
		colnames(full_meth_data) <- full_meth_data[1,]
		full_meth_data[-c(1,2),]
	},mc.cores=future::availableCores()[[1]])
	meth_data <- do.call(rbind,meth_data)
	# meth_data <- meth_data[,colSums(is.na(meth_data)) != nrow(meth_data),]
	row_names <- row.names(meth_data)
	meth_data <- data.frame(as.matrix(apply(meth_data,2,as.numeric)))
	row.names(meth_data) <- row_names
	meth_data
}
 
#===============get meth data from TCGA ELMER library
get_met_for_group <- function(met,group) {
	cur_met <- met[,colData(met)[, "definition", drop = TRUE] %in% c(group)]
	full_meth_data <- assay(cur_met)
	full_meth_data <- full_meth_data[rowSums(is.na(full_meth_data)) != ncol(full_meth_data), ] # result
	colnames(full_meth_data) <- gsub(" ","_",paste(cur_met$definition,colnames(full_meth_data)))
	full_meth_data
}

#============rbind rows with different number of rows, get a list of dataframes with same column, different number of rows
rbind_different_number_of_rows <- function(all_data) {
	all_names <- unique(unlist(lapply(all_data,colnames)))
	all_test_rbind <- do.call(rbind,
        lapply(all_data,
               function(x) data.frame(c(x, sapply(setdiff(all_names, names(x)),
                                                  function(y) NA)))))
}

#================create dir if not exists
create_dir <- function(output_dir_path) {
	if(!dir.exists(output_dir_path)) {
			dir.create(output_dir_path,recursive=T)
			message(paste("output dir:",output_dir_path,"was created"))
		}
		
}

#===============get diff between tomur and normal samples: mean or median
diff_in_means <- function(all_table_methylation_in_tissue,type_of_score,samples_to_remove=c(),type_of_analysis="mean") {
	invert_type <- ifelse(length(samples_to_remove)>0,T,F)
	all_means_score <- aggregate(get(type_of_score)~category+cancer_type,
		data=all_table_methylation_in_tissue[grep(paste(samples_to_remove,collapse="|"),all_table_methylation_in_tissue$sample,invert=invert_type),],get(type_of_analysis))
	all_means_score <- all_means_score[grep("Tumor|Healthy",all_means_score$category),]
	#all_means_score <- data.frame(cancer_type=unique(all_table_methylation_in_tissue$cancer_type),diff_in_means=all_means_score[ c(F,T),3 ]-all_means_score[ c(T,F),3 ],category="diff_means_tumor_normal")
	all_means_score <- data.frame(cancer_type=unique(all_table_methylation_in_tissue$cancer_type),diff_in_means=all_means_score[ c(F,T),3 ],category="diff_means_tumor_normal")
	
	colnames(all_means_score)[2] <- type_of_score
	all_means_score
}

#change names of label according to what ben wanted
change_name_of_label <- function(all_table_methylation_in_tissue,label_name,wanted_label_name)
{	
	all_table_methylation_in_tissue[all_table_methylation_in_tissue$category==label_name,"category"] <- wanted_label_name
	all_table_methylation_in_tissue
}

#===============plot of the relative and equal score for tissues
library(ggplot2)
plot_score <- function(all_table_methylation_in_tissue,type_of_score,samples_to_remove,ylab,atlas_title,type_of_analysis,output_dir_path="",suffix=".jpeg",save_plot=T) {
	levels_order <- na.omit(all_table_methylation_in_tissue[all_table_methylation_in_tissue$category=="diff_means_tumor_normal" ,c("cancer_type",type_of_score)])
	all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"diff_means_tumor_normal",paste0(type_of_analysis,"[cancer]"))#paste0(type_of_analysis,"[cancer] - ",type_of_analysis,"[healthy]"))
	all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Tumor","Cancers (others)")
	all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Low Content Tumor",paste0("Cancers (",paste(samples_to_remove,collapse=","),")"))
	groups_for_plot <- levels(as.factor(all_table_methylation_in_tissue$category))
	colors_for_plot <- if (length(groups_for_plot)==3) c("red", "blue", "green") else c("lightpink2","red", "blue", "green")
	p <- ggplot(all_table_methylation_in_tissue, 
		aes(x=factor(cancer_type,levels=levels_order[order(levels_order[,2],decreasing=T),"cancer_type"]), y=get(type_of_score))) + 
		geom_point(aes(color=category))+
		ggtitle("Cancer Signature Score",subtitle=atlas_title)+
		#ylab("beta[sample] - beta[mean]")+
		ylab(ylab)+
		xlab("Cancer Types")+
		labs(color="Nanopore Plasmas")+
		theme_bw()+
		scale_color_manual(breaks = c(groups_for_plot),
                       values=c(colors_for_plot))+
		geom_hline(yintercept=0, linetype="dashed")+
		theme(axis.text.x = element_text(face = "bold", size = 12, angle = 90,vjust = 0.5, hjust=1),
			axis.text.y = element_text(face = "bold", size = 12),
			axis.title.y = element_text(face = "bold", size = 12),
			legend.title=element_text(size=10,face="bold"),
			legend.text=element_text(size=10,face="bold"),
			plot.title = element_text(hjust = 0.5,face = "bold", size = 15),
			plot.subtitle=element_text(hjust = 0.5,face = "bold.italic", size = 12))
  
  	if (isTRUE(save_plot)) {
  		create_dir(output_dir_path)
  		ggsave(file.path(output_dir_path,paste0("score_plot_",type_of_score,suffix)), p) #,height=5,width=8)
		message(paste("saved",file.path(output_dir_path,paste0("score_plot_",type_of_score,suffix))))
	}
	else {
		return(p)
	}
} 

#========================plot LUAD vs. LUSC and identity line y=xlab
plot_LUAD_vs_LUSC <- function(all_table_methylation_in_tissue,type_of_score,output_dir_path="",suffix=".jpeg",save_plot=T,atlas_title) {
	levels_order <- na.omit(all_table_methylation_in_tissue[all_table_methylation_in_tissue$category=="diff_means_tumor_normal" ,c("cancer_type",type_of_score)])
	groups_for_plot <- levels(as.factor(all_table_methylation_in_tissue$category))
	colors_for_plot <- if (length(groups_for_plot)==3) c("green", "blue", "red") else c("green","blue", "lightpink2", "red")
	all_table_methylation_in_tissue_LUAD_LUSC <- all_table_methylation_in_tissue[which(all_table_methylation_in_tissue$cancer_type %in% c("LUAD","LUSC") & all_table_methylation_in_tissue$category != "diff_means_tumor_normal"),]
	all_table_methylation_in_tissue_LUAD_LUSC_wide <- reshape(all_table_methylation_in_tissue_LUAD_LUSC, idvar=c("sample","category"),timevar = "cancer_type", direction = "wide")
	limits_to_show <- c(min(all_table_methylation_in_tissue_LUAD_LUSC[,c(type_of_score)]),max(all_table_methylation_in_tissue_LUAD_LUSC[,c(type_of_score)]))
	p <- ggplot(all_table_methylation_in_tissue_LUAD_LUSC_wide, aes(x=get(paste0(type_of_score,".LUAD")),
		y=get(paste0(type_of_score,".LUSC")))) + 
		geom_point(aes(color=category))+
  	
		ggtitle("Cancer Signature Score",subtitle=atlas_title)+
		ylab(paste("LUSC",sub("\\."," ",sub("_.*","",type_of_score)),"Weight of Hypo/Hyper-meth CpGs"))+
		xlab(paste("LUAD",sub("\\."," ",sub("_.*","",type_of_score)),"Weight of Hypo/Hyper-meth CpGs"))+
		theme_bw()+
		scale_color_manual(breaks = c(groups_for_plot),
						   values=c(colors_for_plot))+
		theme(axis.text.x = element_text(face = "bold", size = 12),
			axis.text = element_text(face = "bold", size = 12),
			axis.title = element_text(face = "bold", size = 12),
			legend.title=element_blank(),
			legend.text=element_text(size=10,face="bold"),
			plot.title = element_text(hjust = 0.5,face = "bold", size = 15),
			plot.subtitle=element_text(hjust = 0.5,face = "bold.italic", size = 12))
	p <- p+geom_abline()+expand_limits(x=limits_to_show,y=limits_to_show)  
  	if (isTRUE(save_plot)) {
  		create_dir(output_dir_path)
  		ggsave(file.path(output_dir_path,paste0("LUAD_vs_LUSC_",type_of_score,suffix)), p,height=6,width=8)
		message(paste("saved",file.path(output_dir_path,paste0(type_of_score,suffix))))
	}
	else {
		return(p)
	}
} 

#========================run the whole analysis
plot_all_type <- function(atlas_type,title_for_plot,plot_LUAD=T,type_of_analysis="mean",samples_to_remove=c()) {
	output_dir="/vol/sci/bio/data/benjamin.berman/ekushele/filippo/data/new_atlas/results/"
	file_to_write <- "methylation_in_sample_per_tissue_score_table_all_samples_all_normalization.csv"
	#fwrite(all_table_methylation_in_tissue,file.path(output_dir,atlas_type,file_to_write))

	all_table_methylation_in_tissue <- data.frame(fread(file.path(output_dir,atlas_type,file_to_write)))
	if (length(samples_to_remove) > 0) 
		all_table_methylation_in_tissue[grep(paste0(samples_to_remove,collapse="|"),all_table_methylation_in_tissue$sample),c("category")] <- "Low Content Tumor"
	 
	type_of_scores_to_plot <- grep("hypo_hyper",colnames(all_table_methylation_in_tissue),value=T)
	all_table_methylation_in_tissue <- plyr::rbind.fill(all_table_methylation_in_tissue,do.call(cbind,lapply(type_of_scores_to_plot, function(type_of_score) diff_in_means(all_table_methylation_in_tissue=all_table_methylation_in_tissue,type_of_score=type_of_score,type_of_analysis=type_of_analysis,samples_to_remove=samples_to_remove))))
	
	#plot scores
	removed_samples <- ifelse(length(samples_to_remove)>0, paste0("_no_",paste0(samples_to_remove,collapse="_")),"")
	output_dir_path=file.path(output_dir,atlas_type,"plots","score_plots",paste0("score_plots_",type_of_analysis,removed_samples))
	atlas_title=title_for_plot
	lapply(type_of_scores_to_plot,function(type_of_score) 
		plot_score(all_table_methylation_in_tissue,type_of_score,samples_to_remove,
		#	ifelse(any(grepl("normalized",type_of_score)), "beta[sample] - beta[median(healthy)]","beta[sample]"),
			"Methylation Beta\nmedian(cancer)-median(healthy)",
			suffix=".jpeg",output_dir_path=output_dir_path,atlas_title=atlas_title,type_of_analysis=type_of_analysis))
	if(isTRUE(plot_LUAD)) lapply(type_of_scores_to_plot,function(type_of_score) plot_LUAD_vs_LUSC(all_table_methylation_in_tissue,type_of_score,suffix="_all_samples.jpeg",output_dir_path=output_dir_path,atlas_title=atlas_title))

}


#================plot correlation line for given data
plot_correlation <- function (data_to_plot,y_axis_data,x_axis_data,main_title,sub_title,xlab,ylab,color_by_category=T,plot_identity_line=F) {
	text_size=28
	p <- ggplot(data_to_plot, aes(y=get(y_axis_data),
		x=get(x_axis_data)))+   
		#ggtitle(main_title,subtitle=sub_title)+
		xlab(xlab)+
		ylab(ylab)+
		theme_bw()+
		theme(axis.text.x = element_text(face = "bold", size = text_size),
			axis.text = element_text(face = "bold", size = text_size),
			axis.title = element_text(face = "bold", size = text_size-5),
			legend.title=element_text(size=text_size,face="bold"),
			legend.text=element_text(size=text_size,face="bold"),
			plot.title = element_text(hjust = 0.5,face = "bold", size = 12),
			plot.subtitle=element_text(hjust = 0.5,face = "bold.italic", size = 12))
	if (plot_identity_line) {
		p <- p+ geom_abline(intercept=0, slope=1)
	}
	else {
		p <- p + geom_smooth(method = "lm", formula = y ~ x,fullrange=F) +
		ggpmisc::stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x.npc = "right",
               vstep = 0.05) # sets vertical spacing
	}
	if (isTRUE(color_by_category))
		p <- p+geom_point(aes(color=category,size = num_reads),alpha=0.5)+  scale_size_continuous(range = c(5,10))+labs(color="",size="Num reads (M)") +  guides(colour = guide_legend(override.aes = list(size=8)))

			#,breaks=c(3,6,9,12,15), labels=c("3-5M","6M","9M","12M","15M"))

	else
		p <- p+geom_point(aes(size=num_reads),alpha=0.5)+  scale_size_continuous(range = c(1, 7))


	return(p) 
}

#====================== plot tumor fractions vs LUAD score
plot_fraction_tumor_vs_LUAD <- function(all_table_methylation_in_tissue,type_of_score,atlas_title,x_axis_data,xlab,sample_to_remove,output_dir_path="",save_plot=T,suffix=".jpeg",add_sample_label=F) {
	
	all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Tumor","Cancers (others)")
	all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Low Content Tumor",paste0("Cancers (",paste(sample_to_remove,collapse=","),")"))
	
	groups_for_plot <- levels(as.factor(all_table_methylation_in_tissue$category))
	colors_for_plot <-c( "lightpink2", "red","blue")
	if(length(samples_to_remove)==0)
		colors_for_plot <- colors_for_plot[2:3]

	p <- plot_correlation(data_to_plot=all_table_methylation_in_tissue,y_axis_data=type_of_score,x_axis_data=x_axis_data,main_title="LUAD Methylation Score vs CNA tumor fraction",sub_title=atlas_title,xlab=xlab,ylab="LUAD Specific Methylation",color_by_category=T)
	p <- p + scale_color_manual(breaks = c(groups_for_plot),
		values=c(colors_for_plot))

	if(isTRUE(add_sample_label))
		p <-p+ggrepel::geom_text_repel(aes(label = num_reads),box.padding = 0.5,xlim = c(NA, NA),ylim = c(NA, NA))#, position = position_dodge(width=0.9), box.padding = 0.5,  size=6)
		
  	if (isTRUE(save_plot)) {
  		create_dir(output_dir_path)
  		ggsave(file.path(output_dir_path,paste0("LUAD_vs_tumorFrac_",type_of_score,suffix)), p,height=6,width=8)
		message(paste("saved",file.path(output_dir_path,paste0(type_of_score,suffix))))
	}
	else {
		return(p)
	}
} 
library(ggpubr)
# install.packages("ggsignif", repos='http://cran.us.r-project.org')
library(ggsignif)
#===================== plot methylation tumor fractions value against category
plot_methylation_tumor_fractions_for_catergory <- function(all_table_methylation_in_tissue,main_title,sub_title,y_axis_data,xlab,ylab,sample_to_remove,output_dir_path="",save_plot=T,suffix=".pdf",add_sample_label=F) {
	text_size=28
	colors_for_plot <-c( "lightpink2", "red","blue")
	if (length(sample_to_remove) > 0) {
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Tumor","Cancers (others)")
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Low Content Tumor",paste0("Cancers (",paste(sample_to_remove,collapse=","),")"))
	} else {
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Tumor","Cancer")
		colors_for_plot <-colors_for_plot[2:3]
	}
	subset(all_table_methylation_in_tissue,category=="Healthy")[,"nnls_hypo_hyper"]
		groups_for_plot <- levels(as.factor(all_table_methylation_in_tissue$category))
	print(y_axis_data)
	p <- ggplot(all_table_methylation_in_tissue, aes(y=get(y_axis_data),
		x=factor(category,levels=c("Healthy","Cancer")),color=category,size = num_reads))+
		scale_size_continuous(range = c(5,10))+labs(color=colors_for_plot,size="Num reads (M)") +
	 	ggbeeswarm::geom_quasirandom(dodge.width=0.3, alpha=0.5, bandwidth = 2, varwidth = TRUE) +
		#geom_point(aes(color=category,size = num_reads),alpha=0.5)+  scale_size_continuous(range = c(5,10))+labs(color="",size="Num reads (M)") +    
		#ggtitle(main_title,subtitle=sub_title)+
		xlab(xlab)+
		ylab(ylab)+
		theme_bw()+  
		scale_color_manual(breaks = c(groups_for_plot),values=c(colors_for_plot))+
			scale_y_continuous(breaks=c(seq(0,0.6,0.1)),limits=c(0,0.6))+#values=c(0,0.1,0.2,0.3,0.4,0.5))#
 
		#stat_compare_means(comparisons =list(c("Healthy","Cancer")),method="t.test",paired=F,aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))+ 
		geom_signif(comparisons=list(c("Healthy","Cancer")),test="t.test",test.args="less",color="black",textsize=12,size=1,
				map_signif_level = function(p) sprintf("p = %.2g", p))+
		
		theme(axis.text.x = element_text(face = "bold", size = text_size-2),
			axis.text = element_text(face = "bold", size = text_size),
			axis.title = element_text(face = "bold", size = text_size),
			legend.title=element_text(size=12,face="bold"),
			legend.text=element_text(size=12,face="bold"),
			legend.position = 'none', 
			plot.title = element_text(hjust = 0.5,face = "bold", size = 12),
			plot.subtitle=element_text(hjust = 0.5,face = "bold.italic", size = 12))
			
		 y_range=ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
	 print(y_range)


	if (isTRUE(save_plot)) {
  		create_dir(output_dir_path)
		file_name=file.path(output_dir_path,paste0(gsub(" ","_",y_axis_data),suffix))
  		ggsave(file_name, p,height=8,width=4,dpi=600)
		message(paste("saved",file_name))
	}
	else {
		return(p)
	}


}

#====================== plot tumor fractions vs any data: elaboration of plot_fraction_tumor_vs_LUAD
plot_fraction_tumor_vs_data <- function(all_table_methylation_in_tissue,main_title,sub_title,x_axis_data,y_axis_data,xlab,ylab,sample_to_remove,plot_identity_line,output_dir_path="",save_plot=T,suffix=".jpeg",add_sample_label=F,add_horizonal_line=F,add_arrows=F) {
	colors_for_plot <-c( "lightpink2", "red","blue")
	if (length(sample_to_remove) > 0) {
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Tumor","Cancers (others)")
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Low Content Tumor",paste0("Cancers (",paste(sample_to_remove,collapse=","),")"))
	} else {
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Tumor","Cancers")
		colors_for_plot <-colors_for_plot[2:3]
	}
	groups_for_plot <- levels(as.factor(all_table_methylation_in_tissue$category))
	p <- plot_correlation(data_to_plot=all_table_methylation_in_tissue,y_axis_data=y_axis_data,x_axis_data=x_axis_data,
		main_title=main_title,sub_title=sub_title,xlab=xlab,ylab=ifelse(isTRUE(add_arrows),"\n\n\n\n\n",ylab),color_by_category=T,
		plot_identity_line=plot_identity_line)+
		scale_y_continuous(breaks=c(seq(0,0.6,0.1)),limits=c(0,0.6))+#values=c(0,0.1,0.2,0.3,0.4,0.5))#
		scale_x_continuous(breaks=c(seq(0,0.6,0.1)),limits=c(0,0.6))#values=c(0,0.1,0.2,0.3,0.4,0.5))#
	p <- p + scale_color_manual(breaks = c(groups_for_plot),
		values=c(colors_for_plot))
	if(isTRUE(add_sample_label))
		p <-p+ggrepel::geom_text_repel(aes(label = num_reads),box.padding = 0.5,xlim = c(0, NA),ylim = c(NA, NA))#, position = position_dodge(width=0.9), box.padding = 0.5,  size=6)
		
	if(isTRUE(add_horizonal_line))
		p <- p + geom_hline(yintercept=0, linetype="dashed")

	max_healthy <- max(subset(all_table_methylation_in_tissue,category=="Healthy")[,y_axis_data])
	min_healthy <- min(subset(all_table_methylation_in_tissue,category=="Healthy")[,y_axis_data])


	max_cancer <- max(subset(all_table_methylation_in_tissue,category!="Healthy")[,y_axis_data])
	min_cancer <- min(subset(all_table_methylation_in_tissue,category!="Healthy")[,y_axis_data])

	space_healthy_cancer <- min_cancer-max_healthy
	 y_range=ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
	 print(y_range)
		# x_range=ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
	
	ranges=data.frame(y_range,y_end=c(min_healthy,max_cancer),y_data=c(max_healthy+space_healthy_cancer*0.7,min_cancer),labels=c("Healthy","Cancer"))
	text_high <- grid::textGrob("Healthy", gp=gpar(fontsize=13, fontface="bold"),rot = 90) 

	font_size=20
	if(isTRUE(add_arrows)) {
		arrow_plot <- ggplot(ranges)+ geom_segment(aes(x=0, xend = 0,
			yend = y_end, y = y_data), size=1,inherit.aes=F,
               arrow = arrow(length = unit(0.6,"cm")))+
			# scale_y_continuous(limits=y_range) + 
			scale_x_continuous(limits=c(0,0))+
			geom_text(data=ranges,aes(x=0,y=(y_data+y_end)/2,label=labels,angle=90),vjust=-1.5,size=font_size/4)+#,   ylim = c(-Inf, Inf),xlim=c(NA,Inf))+
			coord_cartesian(clip="off",ylim=y_range) + 
			theme_classic()+ylab(ylab)+ggtitle("","")+xlab("")+
			theme(aspect.ratio = 20/1,
				plot.title = element_text(hjust = 0.5,face = "bold", size = font_size),
		
		        axis.title = element_text(face = "bold", size = font_size),
		         axis.text=element_blank(),
		         axis.ticks=element_blank(),
		         axis.line=element_blank())#		ggplot(ranges)+geom_line(aes(x=0,y=y_range[1]),arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), size = 3)

		
		# p <- p + scale_x_continuous(expand = c(0, 0))+coord_cartesian(xlim = c(0, NA),clip="off")+ geom_segment(aes(x=x_range[1], xend = x_range[1] , y=y_range[1], yend = positions[2]), size=1,inherit.aes=F,
               # arrow = arrow(length = unit(0.6,"cm"))) #+
  			# coord_cartesian(clip = "off")
		p <- p + inset_element(arrow_plot, left = 0, bottom = 0.1, right = 0.3, top = 1,align_to = 'full')


		#annotation_custom(text_high,xmin=x_range[1]-0.02,xmax=x_range[1]-0.02,ymin=-0.01,ymax=positions[2]) 
	}
  	if (isTRUE(save_plot)) {
  		create_dir(output_dir_path)
		suffix <-sub("jpeg","pdf",suffix)
		file_name=file.path(output_dir_path,paste0(gsub(" ","_",y_axis_data),"_vs_tumorFrac",suffix))
  		ggsave(file_name, p,height=7,width=10,dpi=600)
		message(paste("saved",file_name))
	}
	else {
		return(p)
	}
} 



#====================== plot tumor fractions vs any data: elaboration of plot_fraction_tumor_vs_LUAD
plot_fraction_tumor_vs_data_all_samples <- function(all_sum_RMSE,all_table_methylation_in_tissue,main_title,sub_title,x_axis_data,y_axis_data,xlab,ylab,sample_to_remove,plot_identity_line,add_sample_label=F,add_horizonal_line=F) {
	all_sum_RMSE$value <- round(all_sum_RMSE$value,4)
	all_sum_RMSE<-  reshape(all_sum_RMSE,idvar="cancer_type",timevar ="variable",direction="wide")
	 annotations <- data.frame(
        xpos = c(0),
        ypos =  c(1),
        annotateText = paste("RMSE(deep):",all_sum_RMSE$value.RMSE_nnls_hypo_hyper_19326_BC01,"\nRMSE(all):",all_sum_RMSE$value.RMSE_nnls_hypo_hyper,"\npearson(all):",all_sum_RMSE$value.pearson_data),
         hjustvar = 0,
         vjustvar = 1) #<- adjust
		
	colors_for_plot <-c( "lightpink2", "red","blue")
	if (length(sample_to_remove) > 0) {
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Tumor","Cancers (others)")
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Low Content Tumor",paste0("Cancers (",paste(sample_to_remove,collapse=","),")"))
	} else {
		all_table_methylation_in_tissue <- change_name_of_label(all_table_methylation_in_tissue,"Tumor","Cancers")
		colors_for_plot <-colors_for_plot[2:3]
	}
	groups_for_plot <- levels(as.factor(all_table_methylation_in_tissue$category))
	p <- plot_correlation(data_to_plot=all_table_methylation_in_tissue,y_axis_data=y_axis_data,x_axis_data=x_axis_data,
		main_title=main_title,sub_title=sub_title,xlab=xlab,ylab=ylab,color_by_category=T,
		plot_identity_line=plot_identity_line)+ggtitle(main_title)+
		theme(plot.title = element_text(hjust = 0.5,face = "bold", size = 20),
				legend.position='none',
				axis.text.x=element_text(angle=90))+
		  
		scale_y_continuous(breaks=c(seq(0,1,0.1)),limits=c(0,1))+#values=c(0,0.1,0.2,0.3,0.4,0.5))#
		scale_x_continuous(breaks=c(seq(0,1,0.1)),limits=c(0,1))#values=c(0,0.1,0.2,0.3,0.4,0.5))#
	p <- p + scale_color_manual(breaks = c(groups_for_plot),
		values=c(colors_for_plot))
	if(isTRUE(add_sample_label))
		p <-p+ggrepel::geom_text_repel(aes(label = length_overlapping_meth),size=6,box.padding = 0.5,xlim = c(0, NA),ylim = c(NA, NA))#, position = position_dodge(width=0.9), box.padding = 0.5,  size=6)
		
	if(isTRUE(add_horizonal_line))
		p <- p + geom_hline(yintercept=0, linetype="dashed") 
	p <-p+ geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText,size=8,fontface="bold"))




		return(p)
	
} 


#======================bar plot of all r^2 from correlation

bar_plot_for_final_figure <- function(data_to_plot,y_axis_category,xlab,ylab,main_title,sub_title,decreasing=F,output_dir_path="",suffix=".jpeg",save_plot=T) {
print(decreasing)
	levels_order_x_axis <-data_to_plot$cell_type[order(data_to_plot[,y_axis_category],decreasing=decreasing)]
	p<-ggplot(data=data_to_plot, aes(y=get(y_axis_category), x=factor(cell_type,levels=levels_order_x_axis))) +
		geom_bar(stat="identity")+
		ggtitle(main_title,subtitle=sub_title)+
		ylab(ylab)+
		xlab(xlab)+
		theme_bw()+
		theme(axis.text.x = element_text(face = "bold", size = 12, angle = 90,hjust = 1),
			axis.text.y = element_text(face = "bold", size = 12),
			axis.title.y = element_text(face = "bold", size = 12),
			legend.title=element_text(size=10,face="bold"),
			legend.text=element_text(size=10,face="bold"),
			plot.title = element_text(hjust = 0.5,face = "bold", size = 15),
			plot.subtitle=element_text(hjust = 0.5,face = "bold.italic", size = 12))
  
	if (isTRUE(save_plot)) {
  		create_dir(output_dir_path)
		suffix=sub("jpeg","pdf",suffix)
		file_name=file.path(output_dir_path,paste0("bar_plots_",y_axis_category,suffix))
  		ggsave(file_name, p,height=5,width=8,dpi=600)
		message(paste("saved",file_name))
	}
	else {
		return(p)
	}
}
