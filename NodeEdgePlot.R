# node edge feature interaction model

# read in interaction df

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)

# read in lead signal document for conversion
leads<-read.table('/Users/tjsears/Code/T1D/interactionPlots/LeadSignals.txt',sep=" ",header=F)
colnames(leads)<-c('Gene','Loci')

# Feature importance
featureImportance<-read.table("/Users/tjsears/Code/T1D/featureImportances/ALL_featureImportanceMatrix_noPCs_jun4_frac_0.3_199_features.txt",sep="\t",header=T,row.names = 1)
#featureImportance<-read.table("/Users/tjsears/Code/T1D/data/catboost_all_NoPCS_feature_interactions.csv",sep=",",header=T,row.names = 1)
#featureImportance<-featureImportance*10000

colnames(featureImportance)<-gsub("SNPS_","",colnames(featureImportance))
rownames(featureImportance)<-gsub("SNPS_","",rownames(featureImportance))
colnames(featureImportance)<-gsub("AA_A_","HLA_",colnames(featureImportance))
rownames(featureImportance)<-gsub("AA_A_","HLA_",rownames(featureImportance))
#colnames(featureImportance)<-gsub("^([^_]*_[^_]*)_", "\\1\n",colnames(featureImportance))
#rownames(featureImportance)<-gsub("^([^_]*_[^_]*)_", "\\1\n",rownames(featureImportance))
colnames(featureImportance)<-gsub("\\.",":",colnames(featureImportance))
rownames(featureImportance)<-gsub("\\.",":",rownames(featureImportance))


for (i in 1:nrow(featureImportance)){
  for (j in 1:nrow(featureImportance)){
    if (i<j){
      featureImportance[i,j]<-0
    }
  }
}

cols=brewer.pal(3,"Set2")
cols=cols[c(3,1,2)]
library(reshape2)

feet_long<-melt(featureImportance)
feet_long$AltVar<-colnames(featureImportance)
feet_long$FinalVar<-paste(feet_long$variable,"/",feet_long$AltVar)
feet_long<-feet_long[feet_long$value>0,]
feet_long<-feet_long[order(feet_long$value,decreasing = T),]

#feet_long<-feet_long[1:400,]

# subset feature importance into entries only containing variables of interest
featureImportance_smaller<-featureImportance[unique(feet_long$variable),unique(feet_long$variable)]

library(scales)
feet_long$FinalVar<-factor(feet_long$FinalVar,levels=feet_long$FinalVar)

featureImportance_smaller$names<-rownames(featureImportance_smaller)

links<-pivot_longer(featureImportance_smaller,cols=colnames(featureImportance_smaller)[1:(ncol(featureImportance_smaller)-1)])
links$value<-round(links$value,digits = 5)

links_table<-links
links_table$finalVar<-paste(links_table$name,links_table$names)
links_table_final<-aggregate(links_table$value,by = list(links_table$name),FUN = sum)
colnames(links_table_final)<-c('name','value')

links<-links[!duplicated(links$value),] #get rid of zeroes
topHits<-links_table_final$name[order(links_table_final$value,decreasing = T)][1:20]
topHits<-c("rs1064173"                  ,   "chr11:2160994:A:T"            , "rs9276235"       ,'rs9268652','chr1:113834946:A:G','chr14:68792124:C:G', 'chr10:8064496:C:T',            
            "rs1391373"                  ,   "chr7:20346715:T:G"          ,   "chr14:98019683:T:C",'chr8:58959618:A:G','chr18:12777326:T:C','chr2:241345437:GACTTT:G','chr1:192541882:G:C','chr11:122727387:G:A','')
  
links<-links[ (links$value>750 | #keep if high mag
               links$names%in%topHits |  # keep if related to a top hit
               links$name%in%topHits) & links$value>350,]  # keep if related to a top hit

#links<-links[links$value>400|links$value<100,]

nodes<-as.data.frame(unique(c(links$names,links$name)))
colnames(nodes)<-c('name')
nodes$type<-ifelse(!grepl('chr',nodes$name),"MHC","NonMHC")
  
# youd have to order by mhc / non mhc
# color by non mhc / mhc 

# chord diagram
library(circlize)
library(RColorBrewer)

cols=brewer.pal(3,"Set2")
cols=cols[c(3,1,2)]

final_var<-links
final_var$FinalVar<-paste(links$names,links$name)

#66C2A5
col_list<-rep('#FC8D62',nrow(links))
col_list<-ifelse(!grepl('chr',final_var$FinalVar),"#8DA0CB",col_list)
col_list[grepl('chr',final_var$FinalVar)]<-"#66C2A5"
col_list[grepl('chr',final_var$FinalVar) & !(grepl('intron',final_var$FinalVar)|grepl('rs',final_var$FinalVar)|grepl('HLA',final_var$FinalVar)|grepl('exon',final_var$FinalVar))]<-"#FC8D62"

unique_names<-c(unique(c(links$name,links$names)))
unique_names_order<-order(unique_names)
gap_determine<-unique_names[order(unique_names)]

gap_determine[grepl('chr',gap_determine)]<-'chr'
gap_determine[!grepl('chr',gap_determine)]<-'mhc'

gaps<-c(1)
for (i in 2:length(gap_determine)){
  if (gap_determine[i]==gap_determine[i-1]){
    gaps[i]<-1
  } else {
    gaps[i-1]<-10
    gaps[i]<-1
  }
}
gaps[length(gaps)]<-10
gaps[1]<-1

# alter links to account for lead signals
match_vec<-match(links$names,leads$Loci,nomatch = 0)
links$names[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

match_vec<-match(links$name,leads$Loci,nomatch = 0)
links$name[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

match_vec<-match(nodes$name,leads$Loci,nomatch = 0)
nodes$name[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

grid.col = ifelse(nodes$type=='MHC','#6457A6','#AB2346')
names(grid.col)<-nodes$name

# respecify names
unique_names<-c(unique(c(links$name,links$names)))

pdf("/Users/tjsears/Code/T1D/figs_may/fig4/interaction.pdf",width=10.9,height=8.23) 

circos.par(gap.after=gaps,canvas.xlim=c(-1.05,1),canvas.ylim=c(-1.1,1))

chordDiagram(annotationTrack = c("grid",'axis'),links,order = unique_names[unique_names_order],grid.col = grid.col,col=col_list,transparency = 0.2,
             preAllocateTracks = list(track.height = 0.2),link.lwd = 0.8, link.lty = 1, link.border = "gray9")

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  if(abs(xplot[2] - xplot[1]) < 25) {
    circos.text(cex=0.8,mean(xlim), ylim[1]+0.25, sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = "gray9")
  } else {
    circos.text(cex=0.95,mean(xlim), ylim[1]+0.25, sector.name, facing = "inside", 
                niceFacing = TRUE, adj = c(0.5, 0), col= "gray9")
  }
}, bg.border = NA)

library(ComplexHeatmap)
# Origin Legend
lgd_origin = Legend(at = c("MHC", "Non-MHC"), type = "points", 
                    legend_gp = gpar(col=c("#6457A6","#AB2346"),cex=3), title_position = "topleft", 
                    title = "Origin",background=NULL,pch=15,size = unit(4, "mm"))

# Interaction Legend
cols_1<-brewer.pal(3,"Set2")
lgd_int = Legend(at = c("MHC / MHC",'Non-MHC / Non-MHC',"MHC / Non-MHC"), type ='points',
                    legend_gp = gpar(col=c("#8DA0CB","#FC8D62","#66C2A5"),cex=3), title_position = "topleft", 
                    title = "Interaction",background=NULL,pch=15,size = unit(4, "mm"))

lgd_list_vertical = packLegend(lgd_origin, lgd_int)

draw(lgd_list_vertical, x = unit(35, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))
circos.clear()


background='white'  
dev.off()


if (FALSE) {


# MHC only

# Feature importance
featureImportance<-read.table("/Users/tjsears/Code/T1D/featureImportances/HLA_featureImportanceMatrix_noPCs_jun4_frac_0.3_199_features.txt",sep="\t",header=T,row.names = 1)
colnames(featureImportance)<-gsub("SNPS_","",colnames(featureImportance))
rownames(featureImportance)<-gsub("SNPS_","",rownames(featureImportance))
colnames(featureImportance)<-gsub("^([^_]*_[^_]*)_", "\\1\n",colnames(featureImportance))
rownames(featureImportance)<-gsub("^([^_]*_[^_]*)_", "\\1\n",rownames(featureImportance))
rownames(featureImportance)<-gsub("\\*",".",rownames(featureImportance))
colnames(featureImportance)<-gsub("\\.",":",colnames(featureImportance))
rownames(featureImportance)<-gsub("\\.",":",rownames(featureImportance))


for (i in 1:nrow(featureImportance)){
  for (j in 1:nrow(featureImportance)){
    if (i<j){
      featureImportance[i,j]<-0
    }
  }
}

cols=brewer.pal(3,"Set2")
cols=cols[c(3,1,2)]
library(reshape2)

feet_long<-melt(featureImportance)
feet_long$AltVar<-colnames(featureImportance)
feet_long$FinalVar<-paste(feet_long$variable,"/",feet_long$AltVar)
feet_long<-feet_long[feet_long$value>0,]
feet_long<-feet_long[order(feet_long$value,decreasing = T),]
feet_long<-feet_long[seq(from=1,to=100,by=1),]

# subset feature importance into entries only containing variables of interest
featureImportance_smaller<-featureImportance[unique(feet_long$variable),unique(feet_long$variable)]

library(scales)
feet_long$FinalVar<-factor(feet_long$FinalVar,levels=feet_long$FinalVar)

featureImportance_smaller$names<-rownames(featureImportance_smaller)

links<-pivot_longer(featureImportance_smaller,cols=colnames(featureImportance_smaller)[1:(ncol(featureImportance_smaller)-1)])

links$value<-round(links$value)

links<-links[!duplicated(links$value),]

links<-links[links$value>1000,]

nodes<-as.data.frame(unique(c(links$names,links$name)))
colnames(nodes)<-c('name')
nodes$type<-ifelse(!grepl('chr',nodes$name),"MHC","NonMHC")

# youd have to order by mhc / non mhc
# color by non mhc / mhc 

# chord diagram
library(circlize)
library(RColorBrewer)

cols=brewer.pal(3,"Set2")
cols=cols[c(3,1,2)]

grid.col = ifelse(nodes$type=='MHC','#6457A6','#AB2346')
names(grid.col)<-nodes$name
final_var<-links
final_var$FinalVar<-paste(links$names,links$name)

#66C2A5
col_list<-rep('#FC8D62',nrow(links))
col_list<-ifelse(!grepl('chr',final_var$FinalVar),"#8DA0CB",col_list)
col_list[grepl('chr',final_var$FinalVar)]<-"#66C2A5"
col_list[grepl('chr',final_var$FinalVar) & !(grepl('intron',final_var$FinalVar)|grepl('rs',final_var$FinalVar)|grepl('HLA',final_var$FinalVar))]<-"#FC8D62"

unique_names<-c(unique(c(links$name,links$names)))
gap_determine<-unique_names[order(unique_names)]

gap_determine[grepl('chr',gap_determine)]<-'chr'
gap_determine[!grepl('chr',gap_determine)]<-'mhc'

gaps<-c(1)
for (i in 2:length(gap_determine)){
  if (gap_determine[i]==gap_determine[i-1]){
    gaps[i]<-1
  } else {
    gaps[i-1]<-10
    gaps[i]<-1
  }
}
gaps[length(gaps)]<-10

pdf("/Users/tjsears/Code/T1D/interactionPlots/MHC_SampleInteraction.pdf",width=9.69,height=6.23) 

circos.par(gap.after=gaps,canvas.xlim=c(-1,1),canvas.ylim=c(-1.3,1.1))

chordDiagram(annotationTrack = c("grid",'axis'),links,order = unique_names[order(unique_names)],grid.col = grid.col,col=col_list,transparency = 0.2,
             preAllocateTracks = list(track.height = 0.2),link.lwd = 0.8, link.lty = 1, link.border = "gray9")

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  if(abs(xplot[2] - xplot[1]) < 30) {
    circos.text(cex=0.65,mean(xlim), ylim[1]+0.3, sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = "gray9")
  } else {
    circos.text(cex=0.95,mean(xlim), ylim[1]+0.3, sector.name, facing = "inside", 
                niceFacing = TRUE, adj = c(0.5, 0), col= "gray9")
  }
}, bg.border = NA)

library(ComplexHeatmap)
# Origin Legend
lgd_origin = Legend(at = c("MHC", "Non-MHC"), type = "points", 
                    legend_gp = gpar(col=c("#6457A6","#AB2346"),cex=3), title_position = "topleft", 
                    title = "Origin",background=NULL,pch=15,size = unit(4, "mm"))

# Interaction Legend
cols_1<-brewer.pal(3,"Set2")
lgd_int = Legend(at = c("MHC / MHC",'Non-MHC / Non-MHC',"MHC / Non-MHC"), type ='points',
                 legend_gp = gpar(col=c("#8DA0CB","#FC8D62","#66C2A5"),cex=3), title_position = "topleft", 
                 title = "Interaction",background=NULL,pch=15,size = unit(4, "mm"))

lgd_list_vertical = packLegend(lgd_origin, lgd_int)

draw(lgd_list_vertical, x = unit(35, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))
circos.clear()


background='white'  
dev.off()







# non MHC

# Feature importance
featureImportance<-read.table("/Users/tjsears/Code/T1D/featureImportances/nonHLA_featureImportanceMatrix_noPCs_jun4_frac_0.3_199_features.txt",sep="\t",header=T,row.names = 1)
colnames(featureImportance)<-gsub("SNPS_","",colnames(featureImportance))
rownames(featureImportance)<-gsub("SNPS_","",rownames(featureImportance))
colnames(featureImportance)<-gsub("^([^_]*_[^_]*)_", "\\1\n",colnames(featureImportance))
rownames(featureImportance)<-gsub("^([^_]*_[^_]*)_", "\\1\n",rownames(featureImportance))
colnames(featureImportance)<-gsub("\\.",":",colnames(featureImportance))
rownames(featureImportance)<-gsub("\\.",":",rownames(featureImportance))


for (i in 1:nrow(featureImportance)){
  for (j in 1:nrow(featureImportance)){
    if (i<j){
      featureImportance[i,j]<-0
    }
  }
}

cols=brewer.pal(3,"Set2")
cols=cols[c(3,1,2)]
library(reshape2)

feet_long<-melt(featureImportance)
feet_long$AltVar<-colnames(featureImportance)
feet_long$FinalVar<-paste(feet_long$variable,"/",feet_long$AltVar)
feet_long<-feet_long[feet_long$value>0,]
feet_long<-feet_long[order(feet_long$value,decreasing = T),]
feet_long<-feet_long[seq(from=1,to=100,by=1),]

# subset feature importance into entries only containing variables of interest
featureImportance_smaller<-featureImportance[unique(feet_long$variable),unique(feet_long$variable)]

library(scales)
feet_long$FinalVar<-factor(feet_long$FinalVar,levels=feet_long$FinalVar)

featureImportance_smaller$names<-rownames(featureImportance_smaller)

links<-pivot_longer(featureImportance_smaller,cols=colnames(featureImportance_smaller)[1:(ncol(featureImportance_smaller)-1)])

links$value<-round(links$value)

links<-links[!duplicated(links$value),]

links<-links[links$value>600,]

nodes<-as.data.frame(unique(c(links$names,links$name)))
colnames(nodes)<-c('name')
nodes$type<-ifelse(!grepl('chr',nodes$name),"MHC","NonMHC")

# youd have to order by mhc / non mhc
# color by non mhc / mhc 

# chord diagram
library(circlize)
library(RColorBrewer)

cols=brewer.pal(3,"Set2")
cols=cols[c(3,1,2)]

grid.col = ifelse(nodes$type=='MHC','#6457A6','#AB2346')
names(grid.col)<-nodes$name
final_var<-links
final_var$FinalVar<-paste(links$names,links$name)

#66C2A5
col_list<-rep('#FC8D62',nrow(links))
col_list<-ifelse(!grepl('chr',final_var$FinalVar),"#8DA0CB",col_list)
col_list[grepl('chr',final_var$FinalVar)]<-"#66C2A5"
col_list[grepl('chr',final_var$FinalVar) & !(grepl('intron',final_var$FinalVar)|grepl('rs',final_var$FinalVar)|grepl('HLA',final_var$FinalVar))]<-"#FC8D62"

unique_names<-c(unique(c(links$name,links$names)))
unique_names_order<-order(unique_names)
gap_determine<-unique_names[order(unique_names)]

gap_determine[grepl('chr',gap_determine)]<-'chr'
gap_determine[!grepl('chr',gap_determine)]<-'mhc'

gaps<-c(1)
for (i in 2:length(gap_determine)){
  if (gap_determine[i]==gap_determine[i-1]){
    gaps[i]<-1
  } else {
    gaps[i-1]<-10
    gaps[i]<-1
  }
}
gaps[length(gaps)]<-10

# alter links to account for lead signals
match_vec<-match(links$names,leads$Loci,nomatch = 0)
links$names[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

match_vec<-match(links$name,leads$Loci,nomatch = 0)
links$name[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

match_vec<-match(nodes$name,leads$Loci,nomatch = 0)
nodes$name[match_vec>0]<-leads$Gene[match_vec[match_vec>0]]

grid.col = ifelse(nodes$type=='MHC','#6457A6','#AB2346')
names(grid.col)<-nodes$name

# respecify names
unique_names<-c(unique(c(links$name,links$names)))

pdf("/Users/tjsears/Code/T1D/interactionPlots/nonMHC_SampleInteraction.pdf",width=9.69,height=6.23) 

circos.par(gap.after=gaps,canvas.xlim=c(-1,1),canvas.ylim=c(-1,1))

chordDiagram(annotationTrack = c("grid",'axis'),links,order = unique_names[order(unique_names)],grid.col = grid.col,col=col_list,transparency = 0.2,
             preAllocateTracks = list(track.height = 0.2),link.lwd = 0.8, link.lty = 1, link.border = "gray9")

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  if(abs(xplot[2] - xplot[1]) < 100) {
    circos.text(cex=0.7,mean(xlim), ylim[1]+0.3, sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = "gray9")
  } else {
    circos.text(cex=0.85,mean(xlim), ylim[1]+0.3, sector.name, facing = "inside", 
                niceFacing = TRUE, adj = c(0.5, 0), col= "gray9")
  }
}, bg.border = NA)

library(ComplexHeatmap)
# Origin Legend
lgd_origin = Legend(at = c("MHC", "Non-MHC"), type = "points", 
                    legend_gp = gpar(col=c("#6457A6","#AB2346"),cex=3), title_position = "topleft", 
                    title = "Origin",background=NULL,pch=15,size = unit(4, "mm"))

# Interaction Legend
cols_1<-brewer.pal(3,"Set2")
lgd_int = Legend(at = c("MHC / MHC",'Non-MHC / Non-MHC',"MHC / Non-MHC"), type ='points',
                 legend_gp = gpar(col=c("#8DA0CB","#FC8D62","#66C2A5"),cex=3), title_position = "topleft", 
                 title = "Interaction",background=NULL,pch=15,size = unit(4, "mm"))

lgd_list_vertical = packLegend(lgd_origin, lgd_int)

draw(lgd_list_vertical, x = unit(35, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))
circos.clear()


background='white'  
dev.off()

}
