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
library(SeuratWrappers)
library(org.Mm.eg.db)
library(GO.db)
library(ggalluvial)
library(ggtranscript)
library(parallel)
library(ggrepel)
library(ComplexUpset)
library(sparseMatrixStats)
library(scales)
library(extrafont)
library(DESeq2)
library(eulerr)
library(inflection)
library(scater)
library(labeling)
library(grid)
library(pheatmap)
library(ggpubr)
library(MASS)
library(cowplot)
library(magick)
library(ggh4x)
options(stringsAsFactors=FALSE)

col_list <- c(c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8"), brewer.pal(12,"Set3")[-c(2, 9)])
text_size <- 8
title_size <- 9
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, face="bold", color="black"), plot.margin=margin(0,-5,0,-5), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))
tag_thm2 <- theme(plot.tag=element_text(size=title_size, face="bold", color="black"), plot.margin=margin(-10,-5,-10,-5), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))

#gtf <- read.table("genes.gtf", sep="\t")
#gtf <- gtf[grep("^chr.*", gtf[, 1]),]
#colnames(gtf) <- c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
#gtf$transcript_id <- gsub(";.*", "", gsub(".*transcript_id ", "", gtf$attributes))
#gtf$gene_id <- gsub(";.*", "", gsub(".*gene_id ", "", gtf$attributes))
#gtf$gene_name <- gsub(";.*", "", gsub(".*gene_name ", "", gtf$attributes))
#gtf <- read.delim("gtf_info.tsv")

oe_rna <- readRDS("osn_rna.rds")
oe_cells_nano <- read.delim("oe_olfr_info.tsv", h=T)
oe_rna$Group <- "Other"
oe_rna$Group[match(oe_cells_nano$Cell[which(oe_cells_nano$Same > 0 & oe_cells_nano$Other == 0)], colnames(oe_rna))] <- "Multi"
oe_rna$Group[which(oe_rna$cell.subtype_fix != "mOSN")] <- "Other"
rec_ak <- data.frame(X=oe_rna@reductions$umap@cell.embeddings[,1], 
	Y=oe_rna@reductions$umap@cell.embeddings[,2], Group=oe_rna$Group)
rec_ak <- rec_ak[order(rec_ak$Group, decreasing=T),]
pak <- wrap_elements(ggplot(rec_ak, aes(x=X, y=Y, color=Group))+geom_point(size=1)+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour=NULL)+
	scale_color_manual(values=c("red", "gray60"), drop=F)+guides(colour="none")+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	annotate("text", x=-5, y=4, label="HBC → mOSN", color="black", size=3, hjust=0, vjust=1)+
	theme(axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	panel.border=element_rect(color="black", fill=NA, linewidth=0.35), legend.key.size=unit(10, "pt"), 
	plot.title=element_text(size=title_size, hjust=0.5, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	axis.text=element_text(size=text_size, colour="black"),  axis.title=element_text(size=title_size, colour="black"), 
	legend.position=c(0.19, 0.9)))+tag_thm

olfr_info_osn <- read.delim("osn_iso_info_all.tsv", h=T)
olfr_info_osn$Info <- paste(olfr_info_osn$Gene, olfr_info_osn$Symbol, olfr_info_osn$TypeC, olfr_info_osn$TypeE, olfr_info_osn$TSS, 
	olfr_info_osn$UTR, olfr_info_osn$CDS, olfr_info_osn$Exon)
olfr_info_oe <- read.delim("oe_iso_info_all.tsv", h=T)
olfr_info_oe$Info <- paste(olfr_info_oe$Gene, olfr_info_oe$Symbol, olfr_info_oe$TypeC, olfr_info_oe$TypeE, olfr_info_oe$TSS, 
	olfr_info_oe$UTR, olfr_info_oe$CDS, olfr_info_oe$Exon)
olfr_info <- data.frame(olfr_info_oe[match(intersect(olfr_info_osn$Info, olfr_info_oe$Info), olfr_info_oe$Info),], Sample="Merge")
olfr_info <- rbind(olfr_info, rbind(data.frame(olfr_info_osn[match(setdiff(olfr_info_osn$Info, olfr_info_oe$Info), olfr_info_osn$Info),], Sample="OSN"), 
	data.frame(olfr_info_oe[match(setdiff(olfr_info_oe$Info, olfr_info_osn$Info), olfr_info_oe$Info),], Sample="OE")))
isoforms <- rbind(data.frame(read.delim("osn_iso_info.tsv", h=T), Sample="OSN"), data.frame(read.delim("oe_iso_info.tsv", h=T), Sample="OE"))
cells_nano <- rbind(data.frame(read.delim("osn_olfr_info.tsv", h=T), Sample="OSN"), data.frame(read.delim("oe_olfr_info.tsv", h=T), Sample="OE"))
iso_trans_test <- rbind(data.frame(read.delim("osn_iso_trans_test.tsv", h=T), Sample="OSN"), data.frame(read.delim("oe_iso_trans_test.tsv", h=T), Sample="OE"))

rec_ccc <- data.frame(table(cells_nano$Other[which(cells_nano$Num > 1)]))
pcc <-wrap_elements(ggplot(rec_ccc, aes(x=Var1, y=Freq))+
	geom_bar(stat="identity", width=0.6, position=position_dodge(0.8), fill=col_list[1])+
	labs(title=NULL, x="Heterologous isoform", y="Number of OSNs")+
	scale_y_continuous(limits=c(0, max(rec_ccc$Freq)+20), expand=c(0, 0))+guides(fill="none")+
	geom_text(aes(x=Var1, y=Freq+2, label=Freq), size=3, vjust=0)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks.y=element_line(linewidth=0.35, color="black"), axis.ticks.x=element_blank(), 
	panel.background=element_rect(0, linetype=0), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

rec_aa <- data.frame(table(table(olfr_info$Gene)))
paa <- wrap_elements(ggplot(rec_aa, aes(x=Var1, y=Freq))+
	geom_bar(stat="identity", width=0.6, position=position_dodge(0.8), fill=col_list[1])+
	labs(title=NULL, x="Number of OR isoforms", y="Number of OR genes")+
	scale_y_continuous(limits=c(0, max(rec_aa$Freq)+70), expand=c(0, 0))+guides(fill="none")+
	geom_text(aes(x=Var1, y=Freq+8, label=Freq), size=3, vjust=0)+
	theme(panel.background=element_rect(0, linetype=0), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks.y=element_line(linewidth=0.35, color="black"), axis.ticks.x=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

# 858 OR gene, 139 have multi-isoforms
terms <- unique(olfr_info$Gene)
rec <- data.frame()
for (term in terms)
{
	ids <- which(olfr_info$Gene == term)[-1]
	if (length(ids) == 0) next
	rec <- rbind(rec, olfr_info[ids,])
}
rec_aba <- data.frame(Var1=c("AS", "ATSS"), Freq=as.numeric(table(rec$TSS)))
rec_abb <- data.frame(table(rec$TypeE[which(rec$TSS == "-")]))
levels(rec_abb[,1])[1] <- "UTR"
rec_abb[,1] <- factor(rec_abb[,1], levels=rev(c("UTR", "A3SS", "A5SS", "RI", "SE", "A3SS+A5SS", "A3SS+RI", "A5SS+RI")))
paba <- ggplot(rec_aba, aes(x=4, y=Freq, fill=Var1))+geom_col()+
	labs(title="Type of isoforms", x=NULL, y=NULL, fill=NULL)+
	geom_text(aes(label=Freq), size=3, position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+scale_x_continuous(limits=c(2.5, 4.5), expand=c(0, 0)) +
	scale_fill_manual(values=col_list[c(6,2)])+
	theme(plot.title=element_text(size=title_size, hjust=-0.5), panel.background=element_blank(), legend.key.size=unit(10, "pt"), 
	plot.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), axis.title=element_blank(), legend.title=element_text(colour="black", size=title_size), 
	legend.text=element_text(colour="black", size=text_size), plot.margin=margin(-10,-10,-10,-10), panel.spacing=unit(0, "pt"))
pa <- paa+inset_element(paba, 0.4, 0.3, 1, 1)+tag_thm

rec_sub <- rec[which(rec$TSS == "-"),]
res_mat <- data.frame(matrix(FALSE, nrow=nrow(rec_sub), ncol=4, dimnames=list(rec_sub$ID, c("A3SS", "A5SS", "RI", "SE"))))
res_mat$A3SS[grep("A3SS", rec_sub$TypeC)] <- TRUE
res_mat$A5SS[grep("A5SS", rec_sub$TypeC)] <- TRUE
res_mat$RI[grep("RI", rec_sub$TypeC)] <- TRUE
res_mat$SE[grep("SE", rec_sub$TypeC)] <- TRUE
res_mat <- res_mat[which(rec_sub$TypeC != "-"),]
pac <- upset(res_mat, colnames(res_mat), name="Overlap", width_ratio=0.25, wrap=T, stripes="white", 
	matrix=(intersection_matrix(geom=geom_point(size=1.5))), 
	base_annotations=list("Intersection size"=intersection_size(fill=col_list[2], color=col_list[2], 
	text=list(size=3))+scale_y_continuous(expand=c(0, 0))), 
	set_sizes=(upset_set_size()+scale_y_reverse(breaks=c(0, 40), expand=c(0, 0))+theme(panel.background=element_blank(), 
	panel.grid=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), 
	axis.title=element_text(size=title_size, color="black"), axis.text.x=element_text(size=text_size, color="black"), 
	axis.text.y=element_blank(), axis.ticks.x=element_line(linewidth=0.35, color="black"), axis.line.x=element_line(linewidth=0.35, color="black"), 
	plot.margin=margin(-5,-5,-5,-5), panel.spacing=unit(0, "pt"))), 
	#queries=list(upset_query(set="A3SS", fill=info_col["A3SS"], color=info_col["A3SS"]), 
	#upset_query(set="A5SS", fill=info_col["A5SS"], color=info_col["A5SS"]), 
	#upset_query(set="RI", fill=info_col["RI"], color=info_col["RI"]), 
	#upset_query(set="SE", fill=info_col["SE"], color=info_col["SE"])), 
	sort_sets="descending", sort_intersections="descending", sort_intersections_by=c("cardinality"), 
	themes=upset_modify_themes(list(intersections_matrix=list(theme_minimal(), theme(panel.background=element_blank(), 
	panel.grid=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), 
	axis.title.x=element_text(size=title_size, color="black"), axis.title.y=element_blank(), axis.text.x=element_blank(), 
	axis.text.y=element_text(size=text_size, color="black"), plot.margin=margin(-5,-5,-5,-5), panel.spacing=unit(0, "pt"))), 
	"Intersection size"=list(theme_minimal(),theme(panel.background=element_blank(), 
	panel.grid=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), 
	plot.margin=margin(-5,-5,-5,-5), panel.spacing=unit(0, "pt"), axis.text=element_text(size=text_size, color="black"), 
	axis.title=element_blank(), axis.text.x=element_blank(), axis.line.y=element_line(linewidth=0.35, color="black"), 
	axis.ticks.y=element_line(linewidth=0.35, color="black"))))))+tag_thm
#pad <- wrap_elements(wrap_plots(list(pac, wrap_elements(ggdraw()+
#	draw_image("figS04_b.png", scale=1)+theme(plot.margin=margin(-5,-5,-5,-20))))))+tag_thm
pad <- wrap_elements(ggdraw()+draw_image("figS04_b.png", scale=1))+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm
pae <- wrap_elements(ggdraw()+draw_image("figS04_c.png", scale=1))+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm
paf <- wrap_elements(ggdraw()+draw_image("figS04_d.png", scale=1))+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm

tag_thm2 <- theme(plot.tag=element_text(size=title_size, face="bold", color="black"), plot.margin=margin(-32,-5,-10,-5), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))
pad <- wrap_elements(ggdraw()+draw_image("figS04_b.png", scale=1)+tag_thm2)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm2
pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pac, pad, pak), nrow=1, widths=c(0.8, 1, 0.5, 0.8))+
	plot_annotation(tag_levels=list(c("A", "", "B", "", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pae, paf), nrow=1, widths=c(1, 1))+
	plot_annotation(tag_levels=list(c("C", "E")), theme=tag_thm))+tag_thm, 
	pblank), ncol=1, heights=c(1, 1, 2.6)), 
	width=210, height=297, dpi=300, units="mm", filename="oe_fl_fig_S04.pdf", limitsize=F)
