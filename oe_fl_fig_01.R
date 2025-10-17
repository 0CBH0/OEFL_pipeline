library(dplyr)
library(future)
library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(ggbeeswarm)
library(ggpointdensity)
library(RColorBrewer)
library(clusterProfiler)
library(patchwork)
library(viridis)
library(scCustomize)
library(philentropy)
library(monocle3)
library(ggrepel)
library(SeuratWrappers)
library(extrafont)
library(cowplot)
library(magick)
options(stringsAsFactors=FALSE)

col_list <- c(c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344"), brewer.pal(12,"Set3")[c(6,5,4,3,1)])
text_size <- 13
title_size <- 15
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, colour="black"), plot.margin=margin(-3,-3,-3,-3), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", colour=NA),  plot.background=element_rect(fill="transparent", colour=NA), 
	legend.box.spacing=unit(0, "pt"))

types <- c("Macrophage", "Monocyte", "Neutrophil", "Tcell", "Bcell1", "Bcell2", "Fibroblast", "Pericyte", "Microvillar", "Ensheathing", "Sustentacular", "Erythrocyte", "Bowman", 
	"mOSN", "iOSN", "INP", "GBC", "HBC")

cell_col <- col_list[1:length(types)]
names(cell_col) <- types
types_osn <- c("mOSN", "iOSN", "INP", "GBC", "HBC")
sce_rna <- readRDS("OER_fix.rds")
sce_nano <- readRDS("oe_nano_fix.rds")
cells <- intersect(colnames(sce_rna), colnames(sce_nano))
sce_rna <- subset(sce_rna, cells=cells)
sce_rna$cell.subtype_fix <- factor(sce_rna$cell.subtype_fix, levels=types)
sce_nano <- subset(sce_nano, cells=cells)
sce_nano$cell.subtype_fix <- factor(sce_nano$cell.subtype_fix, levels=types)
sce_nano$cell.subtype_rna <- sce_rna$cell.subtype_fix[match(colnames(sce_nano), colnames(sce_rna))]

rec_a <- data.frame(X=sce_rna$nFeature_RNA, Y=sce_nano$nFeature_RNA)
#rec_a <- data.frame(X=log10(sce_rna$nFeature_RNA), Y=log10(sce_nano$nFeature_RNA))
pa <- wrap_elements(ggplot(rec_a, aes(x=X, y=Y))+
	geom_pointdensity()+scale_color_viridis()+
	labs(title=NULL, x="Genes of Illumina (10^3)", y="Genes of Nanopore (10^3)", color="Density")+
	stat_smooth(method=lm, se=F, colour="#e13344", linetype="dashed", linewidth=2)+
	scale_x_continuous(limits=c(min(rec_a$X), max(rec_a$X)+0.1), breaks=c(1000, 2000, 3000, 4000), labels=c(1, 2, 3, 4), expand=c(0, 0))+
	scale_y_continuous(limits=c(min(rec_a$Y), max(rec_a$Y)+0.1), breaks=c(1000, 2000, 3000), labels=c(1, 2, 3), expand=c(0, 0))+
	geom_text(data=data.frame(x=4000, y=min(rec_a$Y)+50), aes(x=x,y=y), 
	label=paste0("R=", round(cor(rec_a$X, rec_a$Y), 2)), size=4.5, vjust=0, hjust=0.5)+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.box.spacing=unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))+tag_thm)+tag_thm
rec_b <- data.frame(X=sce_rna$nCount_RNA, Y=sce_nano$nCount_RNA)
#rec_b <- data.frame(X=log10(sce_rna$nCount_RNA), Y=log10(sce_nano$nCount_RNA))
pb <- wrap_elements(ggplot(rec_b, aes(x=X, y=Y))+
	geom_pointdensity()+scale_color_viridis()+
	labs(title=NULL, x="UMI of Illumina (10^3)", y="UMI of Nanopore (10^3)", color="Density")+
	stat_smooth(method=lm, se=F, colour="#e13344", linetype="dashed", linewidth=2)+
	scale_x_continuous(limits=c(min(rec_b$X), max(rec_b$X)+0.1), breaks=c(10000, 20000, 30000), labels=c(10, 20, 30), expand=c(0, 0))+
	scale_y_continuous(limits=c(min(rec_b$Y), max(rec_b$Y)+0.1), breaks=seq(2000, 8000, 2000), labels=c(2, 4, 6, 8), expand=c(0, 0))+
	geom_text(data=data.frame(x=30000, y=min(rec_b$Y)+100), aes(x=x,y=y), 
	label=paste0("R=", round(cor(rec_b$X, rec_b$Y), 2)), size=4.5, vjust=0, hjust=0.5)+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.box.spacing=unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))+tag_thm)+tag_thm

#length(which(sce_rna$cell.subtype_fix == sce_nano$cell.subtype_fix))/ncol(sce_rna)
rec_c <- data.frame()
for (i in 1:length(types))
{
	for (j in 1:length(types))
	{
		cells_rna <- colnames(sce_rna)[which(sce_rna$cell.subtype_fix == types[i])]
		cells_nano <- colnames(sce_nano)[which(sce_nano$cell.subtype_fix == types[j])]
		if (length(cells_rna) == 0 | length(cells_nano) == 0) next
		rec_c <- rbind(rec_c, data.frame(X=types[i], Y=types[j], 
			Val=length(intersect(cells_rna, cells_nano))*100/length(union(cells_rna, cells_nano))))
	}
}
rec_c$X <- factor(rec_c$X, levels=types)
rec_c$Y <- factor(rec_c$Y, levels=types)
pca <- ggplot(rec_c, aes(x=X, y=Y, fill=Val))+geom_tile()+labs(title=NULL, x=NULL, y=NULL, fill="%")+
	scale_fill_viridis()+scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.ticks=element_blank(), 
	strip.text=element_blank(), strip.background=element_blank(), axis.text=element_blank(), 
	plot.title=element_text(size=title_size, colour="black", face="bold", hjust=0.5), 
	legend.box.spacing=unit(0, "pt"), legend.margin=margin(0, 0, 0, 0), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"))+tag_thm
type_x <- intersect(types, sce_rna$cell.subtype_fix)
clsList_x <- data.frame(Term=type_x, Info="Line")
clsList_x$Term <- factor(clsList_x$Term, levels=type_x)
pcb <- ggplot(clsList_x, aes(x=Term, y="", fill=Term))+geom_tile()+theme_minimal()+
	labs(title=NULL, x="Identified with illumina", y=NULL, fill="Line")+scale_fill_manual(values=cell_col[type_x])+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.ticks=element_blank(), 
	strip.text=element_blank(), strip.background=element_blank(), axis.text=element_blank(), 
	axis.title.x=element_text(size=title_size, colour="black"), legend.position="none")+tag_thm
type_y <- intersect(types, sce_nano$cell.subtype_fix)
clsList_y <- data.frame(Term=type_y, Info="Line")
clsList_y$Term <- factor(clsList_y$Term, levels=type_y)
pcc <- ggplot(clsList_y, aes(y=Term, x="", fill=Term))+geom_tile()+theme_minimal()+
	labs(title=NULL, y="Identified with nanopore", x=NULL, fill="Line")+scale_fill_manual(values=cell_col[type_y])+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.ticks=element_blank(), 
	strip.text=element_blank(), strip.background=element_blank(), axis.text=element_blank(), 
	axis.title.y=element_text(size=title_size, colour="black"), legend.position="none")+tag_thm
pc <- wrap_elements(wrap_plots(A=pca, B=pcb, C=pcc, design="CA\n#B", widths=c(1, 30), heights=c(30, 1))+
	plot_layout(guides="collect"))+tag_thm

sce_rna$cell.subtype_fix <- factor(sce_rna$cell.subtype_fix, levels=types)
rec_da <- data.frame(X=sce_rna@reductions$umap@cell.embeddings[,2], 
	Y=-sce_rna@reductions$umap@cell.embeddings[,1], Type=sce_rna$cell.subtype_fix)
pda <- ggplot(rec_da, aes(x=X, y=Y, color=Type))+geom_point(size=1)+
	labs(title="Identified with illumina", x="UMAP1", y="UMAP2", colour="Type")+
	scale_color_manual(values=cell_col, drop=F)+
	guides(colour=guide_legend(ncol=2, override.aes=list(size=4), reverse=T))+
	theme(axis.line=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=1), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black", angle=90, hjust=0.5, vjust=0.5), 
	axis.title=element_text(size=title_size, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, colour="black"), 
	legend.box.spacing=unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	panel.background=element_rect(0, linetype=0))+tag_thm
#sce_rna$cell.subtype_nano <- sce_nano$cell.subtype_fix[match(colnames(sce_rna), colnames(sce_nano))]
#sce_rna$cell.subtype_nano <- factor(sce_rna$cell.subtype_nano, levels=types)
#rec_db <- data.frame(X=sce_rna@reductions$umap@cell.embeddings[,2], 
#	Y=-sce_rna@reductions$umap@cell.embeddings[,1], Type=sce_rna$cell.subtype_nano)
rec_db <- data.frame(Y=sce_nano@reductions$umap@cell.embeddings[,2], 
	X=sce_nano@reductions$umap@cell.embeddings[,1], Type=sce_nano$cell.subtype_fix)
pdb <- ggplot(rec_db, aes(x=X, y=Y, color=Type))+geom_point(size=1)+
	labs(title="Identified with nanopore", x="UMAP1", y="UMAP2", colour="Type")+
	scale_color_manual(values=cell_col, drop=F)+
	guides(colour=guide_legend(ncol=2, override.aes=list(size=4), reverse=T))+
	theme(legend.position="none", panel.border=element_rect(color="black", fill=NA, linewidth=1), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black", angle=90, hjust=0.5, vjust=0.5), 
	axis.title=element_text(size=title_size, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, colour="black"), 
	legend.box.spacing=unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))+tag_thm
rec_dc <- data.frame(table(sce_rna$cell.subtype_fix))
rec_dc$Var1 <- factor(rec_dc$Var1, levels=types)
#rec_dc$Freq <- rec_dc$Freq*100/sum(rec_dc$Freq)
pdc <- ggplot(rec_dc, aes(x=log10(Freq), y=Var1, fill=Var1))+
	geom_bar(stat="identity", position=position_dodge(0.8))+
	scale_fill_manual(values=cell_col, drop=F)+labs(title=NULL, x="#Cell(log10)", y=NULL)+
	scale_x_continuous(breaks=c(0, 3), expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.line.x=element_line(linetype=1,colour="black"), 
	axis.ticks.x=element_line(color="black"), axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title.x=element_text(size=title_size, colour="black"), 
	axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
	axis.title.y=element_blank(), legend.position="none")+tag_thm
rec_dd <- data.frame(table(sce_nano$cell.subtype_fix))
rec_dd$Var1 <- factor(rec_dd$Var1, levels=types)
#rec_dd$Freq <- rec_dd$Freq*100/sum(rec_dd$Freq)
pdd <- ggplot(rec_dd, aes(x=log10(Freq), y=Var1, fill=Var1))+
	geom_bar(stat="identity", position=position_dodge(0.8))+
	scale_fill_manual(values=cell_col, drop=F)+labs(title=NULL, x="#Cell(log10)", y=NULL)+
	scale_x_continuous(breaks=c(0, 3), expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.line.x=element_line(linetype=1,colour="black"), 
	axis.ticks.x=element_line(color="black"), axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title.x=element_text(size=title_size, colour="black"), 
	axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
	axis.title.y=element_blank(), legend.position="none")+tag_thm
pd <- wrap_elements(wrap_plots(list(pda, pdc, pdb, pdd), nrow=1, widths=c(10, 1, 10, 1))+
	plot_layout(guides="collect"))+tag_thm

cell_marker <- read.delim("cellmarker_sel3.txt")
cell_marker <- cell_marker[match(types, cell_marker$Type),]
marker_list <- unique(unlist(strsplit(cell_marker$Marker, split=",")))
marker_list <- rev(intersect(marker_list, rownames(sce_rna)))
rec_ea <- data.frame()
for (type in types)
{
	data.use <- data.frame(t(sce_rna[["SCT"]]@data[marker_list, which(sce_rna$cell.subtype_fix == type)]))
	avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
	pct.exp <- apply(data.use, 2, PercentAbove, threshold=0)
	rec_ea <- rbind(rec_ea, data.frame(Avg.exp=avg.exp, Pct.exp=pct.exp, Gene=marker_list, Type=type))
}
rec_ea$Avg.exp <- as.vector(t(sapply(marker_list, function(x) MinMax(
	scale(log1p(rec_ea[rec_ea$Gene == x, "Avg.exp"])), min=0, max=2.5))))
rec_ea$Pct.exp <- rec_ea$Pct.exp * 100
rec_ea$Pct.exp[rec_ea$Pct.exp < 0.1] <- NA
rec_ea$Group <- "Others"
for (type in types_osn) rec_ea$Group[which(rec_ea$Type == type)] <- "OSN lineage"
rec_ea$Group <- factor(rec_ea$Group, levels=c("OSN lineage", "Others"))
rec_ea$Gene <- factor(rec_ea$Gene, levels=marker_list)
rec_ea$Type <- factor(rec_ea$Type, levels=types)
rec_eb <- data.frame()
for (type in types)
{
	data.use <- data.frame(t(sce_nano[["SCT"]]@data[marker_list, which(sce_nano$cell.subtype_fix == type)]))
	avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
	pct.exp <- apply(data.use, 2, PercentAbove, threshold=0)
	rec_eb <- rbind(rec_eb, data.frame(Avg.exp=avg.exp, Pct.exp=pct.exp, Gene=marker_list, Type=type))
}
rec_eb$Avg.exp <- as.vector(t(sapply(marker_list, function(x) MinMax(
	scale(log1p(rec_eb[rec_eb$Gene == x, "Avg.exp"])), min=0, max=2.5))))
rec_eb$Pct.exp <- rec_eb$Pct.exp * 100
rec_eb$Pct.exp[rec_eb$Pct.exp < 0.1] <- NA
rec_eb$Group <- "Non-OSN lineage"
for (type in types_osn) rec_eb$Group[which(rec_eb$Type == type)] <- "OSN lineage"
rec_eb$Group <- factor(rec_eb$Group, levels=c("OSN lineage", "Non-OSN lineage"))
rec_eb$Gene <- factor(rec_eb$Gene, levels=marker_list)
rec_eb$Type <- factor(rec_eb$Type, levels=types)
pea <- ggplot(rec_ea, aes(x=Gene, y=Type, size=Pct.exp, color=Avg.exp))+
	geom_point()+scale_radius(range=c(0, 3))+
	labs(title="Illumina", x=NULL, y=NULL, color="Avg Exp.", size="Per Exp.")+
	scale_color_gradient(low="white", high="#440255")+
	facet_grid(Group ~ ., scales="free_y", space="free_y")+
	theme(panel.background=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=1), 
	strip.background=element_blank(), strip.text.y=element_blank(), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black"), 
	legend.box.spacing=unit(0, "pt"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, color="black"), plot.margin=margin(-10, -3, -10, -3))
peb <- ggplot(rec_eb, aes(x=Gene, y=Type, size=Pct.exp, color=Avg.exp))+
	geom_point()+scale_radius(range=c(0, 3))+
	labs(title="Nanopore", x=NULL, y=NULL, color="Avg Exp.", size="Per Exp.")+
	scale_color_gradient(low="white", high="#440255")+
	facet_grid(Group ~ ., scales="free_y", space="free_y")+
	theme(legend.position="none", panel.background=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=1), 
	strip.background=element_blank(), strip.clip="off", 
	strip.text=element_text(size=text_size, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black"), 
	axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, color="black"), plot.margin=margin(-10, -3, -10, -3))
clsList <- data.frame(Term=types, Info="Line", Group="Non-OSN lineage")
clsList$Group[match(types_osn, clsList$Term)] <- "OSN lineage"
clsList$Group <- factor(clsList$Group, levels=c("OSN lineage", "Non-OSN lineage"))
clsList$Term <- factor(clsList$Term, levels=types)
pec <- ggplot(clsList, aes(y=Term, x="", fill=Term))+geom_tile()+theme_minimal()+
	facet_grid(Group ~ ., scales="free_y", space="free_y")+
	labs(title=NULL, y=NULL, x=NULL, fill=NULL)+scale_fill_manual(values=cell_col)+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.ticks=element_blank(), strip.clip="off", 
	strip.text=element_blank(), strip.background=element_blank(), axis.text.x=element_blank(), 
	axis.text.y=element_text(size=text_size, color="black"), legend.position="none", plot.margin=margin(-10, -3, -10, -3))
pe <- wrap_elements(wrap_plots(list(pec, pea, peb), nrow=1, widths=c(1, 60, 60))+plot_layout(guides="collect"))+tag_thm

pp <- wrap_elements(ggdraw()+draw_image("fig01_a.png", scale=1.15))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(pp+plot_annotation(tag_levels=list(c("A")), theme=tag_thm)+tag_thm)+tag_thm, 
	wrap_elements(wrap_plots(list(pa, pb, pc), nrow=1, widths=c(1, 1, 1.1))+
	plot_annotation(tag_levels=list(c("B", "C", "E")), theme=tag_thm)+tag_thm)+tag_thm, 
	wrap_elements(pd+plot_annotation(tag_levels=list(c("D")), theme=tag_thm)+tag_thm)+tag_thm, 
	wrap_elements(pe+plot_annotation(tag_levels=list(c("F")), theme=theme(plot.margin=margin(-10, -3, -10, -3)))+tag_thm)+tag_thm), 
	ncol=1, heights=c(2.5, 3, 4, 4))+tag_thm, width=13, height=16, dpi=200, filename="oe_fl_fig_01.png", limitsize=F)




