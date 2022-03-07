source("/Users/shario/Google Drive/My Drive/BenBermanLab/nonvolition/newDeeptools/functions_atlas.R")
setwd("/Users/shario/Google Drive/My Drive/BenBermanLab/proj-filippo/shari_analysis/DataFromFilo_2022/rerunNew700filt_270222/outDeeptools/")

library(argparse)
require(pals)
library(circlize)
library(dplyr)
library(ComplexHeatmap)
library(gsubfn)

#seqType<-"nanopore"
#subType<-"original"
#nanoType<-"sub2M"
#illType<-"sub4M"

######funtions#####
plotTFTable<-function(nanoType,illType,TF="CTCF"){
  tables.list.nano<-list.files(path = paste0("nanopore/",nanoType,"/","tables/"), pattern = ".txt", full.names = T)
  short=ifelse(nanoType=="original","Orig",ifelse(nanoType=="sub2M","sub2M",NA))
  tumorFracFileNano<-list.files(path = paste0(dirname(getwd())), pattern = ".txt", full.names = T)
  #tumorFracFileNano<-grep("percCNA_015_maxCN3setto0",tumorFracFileNano, value = T)
  tumorFracFileNano<-grep("percCNA_LT015setto0showingMaxCN7",tumorFracFileNano, value = T)
  tumorFracFileNano<-grep(paste0("nanopore",short),tumorFracFileNano, value = T,ignore.case = T)
  
  tumorFracTabNano<-read.delim(tumorFracFileNano,header=T,sep="\t")
  tumorFracTabNano$sample<-gsub(".sub2M|.sub4M|.subNANO","",tumorFracTabNano$sample)
  
  tables.list.ill<-list.files(path = paste0("illumina/",illType,"/","tables/"), pattern = ".txt", full.names = T)
  short=ifelse(illType=="original","Orig",ifelse(illType=="sub4M","Sub4M",ifelse(illType=="subNANO","subNANO",NA)))
  tumorFracFileIll<-list.files(path = paste0(dirname(getwd())), pattern = ".txt", full.names = T)
  tumorFracFileIll<-grep("percCNA_LT015setto0showingMaxCN7",tumorFracFileIll, value = T)
  tumorFracFileIll<-grep(paste0("illumina",short),tumorFracFileIll, value = T,ignore.case = T)
  
  tumorFracTabIll<-read.delim(tumorFracFileIll,header=T,sep="\t")
  tumorFracTabIll$sample<-gsub(".sub2M|.sub4M|.subNANO","",tumorFracTabIll$sample)
  
  tumorFracTab<-rbind(tumorFracTabNano,tumorFracTabIll)
  tumorFracTab$type<-ifelse(grepl("ILL",tumorFracTab$sample),"Illumina WGS","Nanopore WGS")
  tumorFracTab$library<-ifelse(grepl("ILL",tumorFracTab$sample),illType,nanoType)
  
  tumorFracTab$disease<-ifelse(grepl("19_326|BC01|BC08|BC09|BC10|BC11",tumorFracTab$sample),"Cancer",
                               ifelse(grepl("BC02|BC03|BC04|BC05|BC12|HU002",tumorFracTab$sample),"Healthy",NA))
  
    CTCFtable<-lapply(c(tables.list.nano, tables.list.ill), function(tab){
     print(tab)
      sampleName<-gsub(paste0(".Norm-RPGC-RemoveNotPrimary-min130max155.hg38.",TF,".profile1.txt"),"",basename(tab))
      sampleName<-gsub(".sub2M|.sub4M|.subNANO","",sampleName)
      #type<-ifelse(grepl("nanopore",tab),"nanopore",ifelse(grepl("illumina",tab),"illumina",NA))
      df<-read.delim(tab,sep="\t", check.names = F)
      colnames(df)<-df[grep("bins",df$`bin labels`),]
      df<-df[grep("min130max155",df$bins),-c(2)]
      #df<-df[rowSums(is.na(df))<ncol(df),]
      df <- df[,colSums(is.na(df))<nrow(df)]
      df[1,1]<-sampleName
      return(df)
      })
  
  CTCFtable<-do.call(rbind.data.frame, CTCFtable)
  
  CTCFtable<-merge(CTCFtable,tumorFracTab, by.x = c("bins"), by.y = c("sample"))
  CTCFtable <- CTCFtable[,colSums(is.na(CTCFtable))<nrow(CTCFtable)]
  rownames(CTCFtable)<-CTCFtable$bins
  CTCFtable<-CTCFtable[,-c(1)]
  CTCFtable$sampleType<-gsub("_ILL|.sub4M|.subNANO|.sub2M","",CTCFtable$sampleType)
  CTCFtable<-CTCFtable[order(CTCFtable$type,CTCFtable$tumorFrac),]
  CTCFtable<-CTCFtable[-grep("19_326|BC12",rownames(CTCFtable)),]
  
  rownames(CTCFtable)<-gsub(".HAC","",rownames(CTCFtable))
  rownames(CTCFtable)<-gsub("HU002_","HU005.",rownames(CTCFtable))
  rownames(CTCFtable)<-gsub("BC","ISPRO.bc",rownames(CTCFtable))
  
 # newNames<-c("HU005.10","HU005.11","HU005.12","ISPRO.bc01","ISPRO.bc02","ISPRO.bc03","ISPRO.bc04","ISPRO.bc05","ISPRO.bc08","ISPRO.bc09","ISPRO.bc10","ISPRO.bc11","ISPRO.S1")
#  oldNames<-c("HU002_10","HU002_11","HU002_12","BC01","BC02","BC03","BC04","BC05","BC08","BC09","BC10","BC11","19_326")
  return(CTCFtable)
}

makeTFHeatmap<-function(TFtable=CTCFtableForPlotSub,TF="CTCF",column101=100,xlimMin=100,xlimMax=400, binSize=10,font=25, 
                        save=F,w=5, h=10, name, title,gap=0.5){
  
  TFtable<-subset(TFtable, select=-c(maxCN,n_initialNormalFrac,ploidy,percentGenomeCNA,BamReadNum,OrigTumorFrac))
  labCol <- c(seq(0, ncol(TFtable)-5, 1))
  labCol[labCol %% 500 != 0] <- ""
  labCol<-c(labCol[1:500],"sampleType","type","library","tumorFrac","disease")
  #column101=column101
  colnames(TFtable)<-labCol
  colnames(TFtable)[1]<-"-2500"
  colnames(TFtable)[100]<-"-1500"
  colnames(TFtable)[150]<-"-1000"
  colnames(TFtable)[200]<-"-500"
  colnames(TFtable)[250]<-TF
  colnames(TFtable)[300]<-"500"
  colnames(TFtable)[350]<-"1000"
  colnames(TFtable)[400]<-"1500"
  colnames(TFtable)[500]<-"2500"
  maxColumnName<-colnames(TFtable)[xlimMax]
  TFtable<-TFtable[,c(xlimMin:xlimMax,grep("sampleType|type|library|tumorFrac|disease",colnames(TFtable)))]
  colnames(TFtable)<-ifelse(grepl("\\.",colnames(TFtable)),"",colnames(TFtable))
  TFtable$disease<-factor(TFtable$disease, levels = c("Healthy","Cancer"))
  TFtable<-TFtable[order(TFtable$disease,TFtable$tumorFrac),]
  
  hmNanoporeDF<-subset(TFtable,TFtable$type=="Nanopore WGS")
  hmNanoporeDF$group<-ifelse(grepl("ISPRO",rownames(hmNanoporeDF)),"ISPRO","HU")
  hmNanoporeDF$group<-factor(hmNanoporeDF$group, levels = c("HU","ISPRO"))
  hmNanoporeDF<-hmNanoporeDF[order(hmNanoporeDF$disease,hmNanoporeDF$group,hmNanoporeDF$tumorFrac),]
  
  #hmNanoporeDF$disease<-factor(hmNanoporeDF$disease, levels = c("Cancer","Healthy"))
  #hmNanoporeDF<-hmNanoporeDF[order(colnames(hmNanoporeDF)),]
  
  #hmNanoporeDF$disease<-factor(hmNanoporeDF$disease, levels = c("Healthy","Cancer"))
  #hmNanoporeDF<-hmNanoporeDF[rowSums(is.na(hmNanoporeDF)) != ncol(hmNanoporeDF), ]
  
  hmIlluminaDF<-subset(TFtable,TFtable$type=="Illumina WGS")
  
  dfTMP<-hmNanoporeDF
  rightAnnNano=HeatmapAnnotation(empty = anno_empty(border = FALSE, width = unit(1,"mm")),
                             annotation_name_gp= gpar(fontsize = font, col="white"),
                             show_annotation_name = FALSE,
                             simple_anno_size = unit(0.5, "cm"),
                             disease = dfTMP$disease,
                             #library = list_df$hmNanopore$library,
                             tumorFrac = dfTMP$tumorFrac,
                             col = list(disease = c(Healthy =  cols25()[14], Cancer = cols25()[15]),
                                        #library = c(`Nanopore plasma WGS` ="#f8f0d2", `Illumina plasma WGS` ="#d1f2d2",`Tissue WGBS` =cols25()[11]),
                                        tumorFrac = colorRamp2(c(min(dfTMP$tumorFrac, na.rm = T),
                                                                 max(dfTMP$tumorFrac,na.rm = T)), c("white", "black"))),
                             which="row",
                             show_legend = FALSE,
                             annotation_legend_param = list(title_gp = gpar(fontsize = font, color="white"),
                                                            labels_gp = gpar(fontsize = font, color="white")),
                             gap = unit(1.65, "mm"))
  dfTMP<-hmIlluminaDF
  rightAnnILL=HeatmapAnnotation(empty = anno_empty(border = FALSE, width = unit(1,"mm")),
                                 annotation_name_gp= gpar(fontsize = font, col="black"),
                                 show_annotation_name = TRUE,
                                 simple_anno_size = unit(0.5, "cm"),
                                 disease = dfTMP$disease,
                                 #library = list_df$hmNanopore$library,
                                 tumorFrac = dfTMP$tumorFrac,
                                 col = list(disease = c(Healthy =  cols25()[14], Cancer = cols25()[15]),
                                            #library = c(`Nanopore plasma WGS` ="#f8f0d2", `Illumina plasma WGS` ="#d1f2d2",`Tissue WGBS` =cols25()[11]),
                                            tumorFrac = colorRamp2(c(min(dfTMP$tumorFrac, na.rm = T), 
                                                                     max(dfTMP$tumorFrac,na.rm = T)), c("white", "black"))),
                                 which="row",
                                 show_legend = FALSE,
                                 annotation_legend_param = list(title_gp = gpar(fontsize = font, color="white"),
                                                                labels_gp = gpar(fontsize = font, color="white")),
                                 gap = unit(1.65, "mm"))
  
 
  hmNanoporeMat=as.matrix(hmNanoporeDF[,c(1:grep(paste0("^",maxColumnName,"$"),colnames(hmNanoporeDF)))])
  colnames(hmNanoporeMat)<-ifelse(grepl("\\.",colnames(hmNanoporeMat)),"",colnames(hmNanoporeMat))
  
  htNanopore = Heatmap(hmNanoporeMat,
                       name = "Nanopore\n WGS",
                       #split = matTF$library,
                       col = colorRamp2(c(min(hmNanoporeMat,na.rm = T),max(hmNanoporeMat,na.rm = T)), c("white", "red")),
                       right_annotation = rightAnnNano,
                       cluster_rows = F, cluster_columns = F, 
                       show_row_names = T,row_names_side = "left",row_names_gp = gpar(fontsize=font),
                       row_title = "Nanopore WGS",
                       na_col = "white",  show_heatmap_legend = FALSE,
                       column_names_gp = gpar(fontsize=font),row_title_gp = gpar(fontsize=font, fontface="bold"),
                       use_raster = T,
                       border_gp = gpar(col = "black", lty = 1, lwd = 3),
                       width = unit(8, "cm"),
                       height = unit(10, "cm")
  )
  
  hmIlluminaMat=as.matrix(hmIlluminaDF[,c(1:grep(paste0("^",maxColumnName,"$"),colnames(hmIlluminaDF)))])
  colnames(hmIlluminaMat)<-ifelse(grepl("\\.",colnames(hmIlluminaMat)),"",colnames(hmIlluminaMat))
  rownames(hmIlluminaMat)<-gsub("_ILL","",rownames(hmIlluminaMat))
  
  htIllumina = Heatmap(hmIlluminaMat, 
                        name = "Illumina\n WGS", 
                        col = colorRamp2(c(min(hmNanoporeMat,na.rm = T), max(hmNanoporeMat,na.rm = T)), c("white", "red")),
                        right_annotation = rightAnnILL, 
                        row_names_side = "left",row_names_gp = gpar(fontsize=font),
                        cluster_rows = F, cluster_columns = F, 
                        show_row_names = T, row_title = "Illumina WGS",
                        na_col = "white", show_heatmap_legend = FALSE,
                        column_names_gp = gpar(fontsize=font),row_title_gp = gpar(fontsize=font, fontface="bold"),
                        use_raster = T,
                        border_gp = gpar(col = "black", lty = 1, lwd = 3),
                        width = unit(8, "cm"),
                        height = unit(4, "cm")
                        #column_names_rot = 45
  )
  
  htlist<-htNanopore %v% htIllumina

  if(save==T){
    jpeg(paste0(name,".jpeg"),width = w, height = h, units = "in",res = 2000)
    print(paste0("heatmap will be output to", paste0(name,".jpeg"), "in ", getwd()))
    }
    
    draw(htlist,
         column_title= title,
         #show_heatmap_legend = FALSE,
         legend_labels_gp=gpar(fontsize=font),
         legend_title_gp=gpar(fontsize=font),  
         #merge_legend = TRUE,
         heatmap_row_names_gp=gpar(fontsize=font),
         heatmap_column_names_gp=gpar(fontsize=font),
         column_title_gp = gpar(fontsize = font+5, fontface="bold"),
         row_title_gp = gpar(fontsize = font+5, fontface="bold"),
         # legend_grouping = "original",
         # annotation_legend_list = lgd_list,
         # heatmap_legend_side = "right", annotation_legend_side = "right",
         #padding = unit(c(15, 20, 2, 2), "mm"),
         ht_gap = unit(gap, "cm"),
         #auto_adjust = FALSE,
         use_raster = T) 
    
    if(save==T)
    {
      dev.off()
    }
  
  
}



#nano orignal and illumina original
CTCFtableForPlotSub<-plotTFTable(nanoType = "original", illType = "original", TF="CTCF")
tmp<-CTCFtableForPlotSub[,c(501:504,506:510)]
tmp$sampleName<-rownames(tmp)
tmp<-tmp[,c("sampleName","sampleType","type","library","maxCN","tumorFrac","ploidy","percentGenomeCNA","BamReadNum")]
write.table(tmp,"nanoOrigIlluminaOrig.TFtable.txt", sep = "\t",quote = F, row.names = F)
makeTFHeatmap(TFtable=CTCFtableForPlotSub,column101=100,xlimMin=100,xlimMax=400, 
              binSize=10,font=25,gap=0.5, TF="CTCF",
              save=F,w=7, h=10, name="nanoOrigIlluminaOrig", 
              title="               Fragment Coverage")


#nano sub2M and illumina sub4M
CTCFtableForPlotSub<-plotTFTable(nanoType = "sub2M", illType = "sub4M", TF="CTCF")
tmp<-CTCFtableForPlotSub[,c(501:504,506:510)]
tmp$sampleName<-rownames(tmp)
tmp<-tmp[,c("sampleName","sampleType","type","library","maxCN","tumorFrac","ploidy","percentGenomeCNA","BamReadNum")]
write.table(tmp,"nanoSub2MIlluminasub4M.TFtable.txt", sep = "\t",quote = F, row.names = F)
makeTFHeatmap(TFtable=CTCFtableForPlotSub,column101=100,xlimMin=100,xlimMax=400, 
              binSize=10,font=25,gap=0.5,
              save=T,w=7, h=10, name="nanoSub2MIlluminasub4M", 
              title="               Fragment Coverage")
#CTCFtableForPlotSub[,c(501:504,506:510)]

#nano original and illumina subNANO
CTCFtableForPlotSub<-plotTFTable(nanoType = "original", illType = "subNANO", TF="CTCF")
tmp<-CTCFtableForPlotSub[,c(501:504,506:510)]
tmp$sampleName<-rownames(tmp)
tmp<-tmp[,c("sampleName","sampleType","type","library","maxCN","tumorFrac","ploidy","percentGenomeCNA","BamReadNum")]
write.table(tmp,"nanoOrigIlluminasubNANO.TFtable.txt", sep = "\t",quote = F, row.names = F)
makeTFHeatmap(TFtable=CTCFtableForPlotSub,column101=100,xlimMin=100,xlimMax=400, 
              binSize=10,font=25,gap=0.5,
              save=T,w=7, h=10, name="nanoOrigIlluminasubNANO", 
              title="               Fragment Coverage")




###legend
lh=7
lw=7
lgd_disease = Legend(at = 1:3, title = "Disease       ", labels=c("Healthy","Cancer", "Adjacent Lung"), 
                     legend_gp = gpar(fill=c(cols25()[14],  cols25()[15],cols25()[22])), labels_gp = gpar(fontsize = font), 
                     title_gp = gpar(fontsize=font, fontface="bold"), title_gap = unit(4, "mm"),
                     legend_height = unit(lh, "cm"),
                     legend_width = unit(lw, "cm"),
                     grid_width = unit(0.5, "cm"), grid_height = unit(1, "cm"),
                     border = "black", direction = "horizontal", title_position = "lefttop",nrow = 1)

CTCFtableForPlotSub<-plotTFTable(nanoType = "original", illType = "original", TF="CTCF")
minCTCF<-min(na.omit(CTCFtableForPlotSub[,1:500]))
maxCTCF<-max(na.omit(CTCFtableForPlotSub[,1:500]))
colCTCF = colorRamp2(c(minCTCF,maxCTCF), c("white", "red"))

col_fun = colCTCF
lgd_hp = Legend(at = seq(0,1.5,0.5), title = "Fragment\nCoverage     ", col_fun = col_fun, labels_gp = gpar(fontsize = font), 
                title_gp = gpar(fontsize=font, fontface="bold"), title_gap = unit(4, "mm"),
                legend_height = unit(lh, "cm"),
                legend_width = unit(lw, "cm"),grid_height = unit(1, "cm"),
                grid_width = unit(1, "cm"),border = "black",direction = "horizontal", title_position = "lefttop")#,  title_position = "topcenter")

col_fun_meth = colorRamp2(c(0,1), c("yellow", "blue4"))
lgd_meth = Legend(at = seq(0,1), title = "DNA\nMethylation ", col_fun = col_fun_meth, labels_gp = gpar(fontsize = font), 
                  title_gp = gpar(fontsize=font, fontface="bold"), title_gap = unit(4, "mm"),
                  legend_height = unit(lh, "cm"),
                  legend_width = unit(lw, "cm"),grid_height = unit(1, "cm"),
                  grid_width = unit(1, "cm"),border = "black",direction = "horizontal", title_position = "lefttop")#, title_position = "topcenter")

tmp<-CTCFtableForPlotSub[,c(501:504,506:510)]
tf<-tmp$tumorFrac
#min(tf)
#max(tf)
col_fun_tf = colorRamp2(c(min(tf),0.4), c("white", "black"))
lgd_tf = Legend(at = seq(min(tf),0.4,0.2), title = "Tumor\nFraction       ", col_fun = col_fun_tf, labels_gp = gpar(fontsize = font),
                title_gp = gpar(fontsize=font, fontface="bold"), #title_gap = unit(4, "mm"), 
                legend_height = unit(lh, "cm"),
                legend_width = unit(lw, "cm"),grid_height = unit(1, "cm"),
                grid_width = unit(1, "cm"),border = "black",direction = "horizontal", title_position = "lefttop")

#lgd_list=c(lgd_hp, lgd_meth, lgd_tf, lgd_disease, lgd_lib)

dev.off()

pd=packLegend( 
               
               
               lgd_hp,
               lgd_meth,
               lgd_tf,
               lgd_disease,
               
               max_height = unit(10, "cm"), max_width = unit(20, "cm"), gap = unit(1, "cm"),
               column_gap = unit(1, "cm"), row_gap = unit(1, "cm"), direction = "horizontal")

draw(pd)
dev.off()

png(paste0("legend_hm.png"),width = 10, height = 8, units = "in",res = 300)
draw(pd)
dev.off()
