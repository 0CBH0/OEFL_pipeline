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
text_size <- 13
title_size <- 15
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, color="black"), plot.margin=margin(0,-5,0,-5), panel.spacing=unit(0, "pt"), 
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
	annotate("text", x=-4, y=4, label="HBC → mOSN", color="black", size=4, hjust=0, vjust=1)+
	theme(axis.line=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=0.8), 
	plot.title=element_text(size=title_size, hjust=0.5, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	axis.text=element_text(size=text_size, colour="black"),  axis.title=element_text(size=title_size, colour="black"), 
	legend.position=c(0.17, 0.9)))+tag_thm

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
	geom_text(aes(x=Var1, y=Freq+2, label=Freq), size=4, vjust=0)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

rec_aa <- data.frame(table(table(olfr_info$Gene)))
paa <- wrap_elements(ggplot(rec_aa, aes(x=Var1, y=Freq))+
	geom_bar(stat="identity", width=0.6, position=position_dodge(0.8), fill=col_list[1])+
	labs(title=NULL, x="Isoforms of gene", y="Number")+
	scale_y_continuous(limits=c(0, max(rec_aa$Freq)+70), expand=c(0, 0))+guides(fill="none")+
	geom_text(aes(x=Var1, y=Freq+8, label=Freq), size=4, vjust=0)+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
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
	geom_text(aes(label=Freq), size=4, position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+scale_x_continuous(limits=c(2.5, 4.5), expand=c(0, 0)) +
	scale_fill_manual(values=col_list[c(6,2)])+
	theme(plot.title=element_text(size=title_size, hjust=-0.5), panel.background=element_blank(), 
	plot.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), axis.title=element_blank(), legend.title=element_text(colour="black", size=title_size), 
	legend.text=element_text(colour="black", size=text_size), plot.margin=margin(-10,-10,-10,-10), panel.spacing=unit(0, "pt"))
pa <- wrap_elements(paa+inset_element(paba, 0.35, 0.3, 1, 1)+tag_thm)+tag_thm

rec_sub <- rec[which(rec$TSS == "-"),]
res_mat <- data.frame(matrix(FALSE, nrow=nrow(rec_sub), ncol=4, dimnames=list(rec_sub$ID, c("A3SS", "A5SS", "RI", "SE"))))
res_mat$A3SS[grep("A3SS", rec_sub$TypeC)] <- TRUE
res_mat$A5SS[grep("A5SS", rec_sub$TypeC)] <- TRUE
res_mat$RI[grep("RI", rec_sub$TypeC)] <- TRUE
res_mat$SE[grep("SE", rec_sub$TypeC)] <- TRUE
res_mat <- res_mat[which(rec_sub$TypeC != "-"),]
pac <- upset(res_mat, colnames(res_mat), name="Overlap", width_ratio=0.25, wrap=T, stripes="white", 
	base_annotations=list("Intersection size"=intersection_size(fill=col_list[2], color=col_list[2], 
	text=list(size=4))+scale_y_continuous(expand=c(0, 0))), 
	set_sizes=(upset_set_size()+scale_y_reverse(breaks=c(0, 40), expand=c(0, 0))+theme(panel.background=element_blank(), 
	panel.grid=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), 
	axis.title=element_text(size=title_size, color="black"), axis.text.x=element_text(size=text_size, color="black"), 
	axis.text.y=element_blank(), axis.ticks.x=element_line(color="black"), axis.line.x=element_line(color="black"), 
	plot.margin=margin(-5,-5,-5,-5), panel.spacing=unit(0, "pt"))), 
	#queries=list(upset_query(set="A3SS", fill=info_col["A3SS"], color=info_col["A3SS"]), 
	#upset_query(set="A5SS", fill=info_col["A5SS"], color=info_col["A5SS"]), 
	#upset_query(set="RI", fill=info_col["RI"], color=info_col["RI"]), 
	#upset_query(set="SE", fill=info_col["SE"], color=info_col["SE"])), 
	sort_sets="descending", sort_intersections="descending", sort_intersections_by=c("cardinality"), 
	themes=upset_modify_themes(list(intersections_matrix=list(theme_minimal(),theme(panel.background=element_blank(), 
	panel.grid=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), 
	axis.title.x=element_text(size=title_size, color="black"), axis.title.y=element_blank(), axis.text.x=element_blank(), 
	axis.text.y=element_text(size=text_size, color="black"), plot.margin=margin(-5,-5,-5,-5), panel.spacing=unit(0, "pt"))), 
	"Intersection size"=list(theme_minimal(),theme(panel.background=element_blank(), 
	panel.grid=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), 
	plot.margin=margin(-5,-5,-5,-5), panel.spacing=unit(0, "pt"), axis.text=element_text(size=text_size, color="black"), 
	axis.title=element_blank(), axis.text.x=element_blank(), axis.line.y=element_line(color="black"), 
	axis.ticks.y=element_line(color="black"))))))+tag_thm
pad <- wrap_elements(wrap_plots(list(pac, wrap_elements(ggdraw()+draw_image("figS04_b.png", scale=1)+tag_thm))))+tag_thm
pae <- wrap_elements(ggdraw()+draw_image("figS04_c.png", scale=1))+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8))+tag_thm
paf <- wrap_elements(ggdraw()+draw_image("figS04_d.png", scale=1))+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8))+tag_thm

pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pad, pak), nrow=1, widths=c(0.8, 1.5, 0.8))+
	plot_annotation(tag_levels=list(c("A", "B", "C")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pae, paf), nrow=1, widths=c(1, 1))+
	plot_annotation(tag_levels=list(c("D", "E")), theme=tag_thm))+tag_thm), 
	ncol=1, heights=c(1, 1)), width=13, height=8, dpi=200, filename="oe_fl_fig_S04.png", limitsize=F)


















cells <- colnames(osn_rna)
terms <- rownames(osn_rna[["RNA"]]@counts)[grep("^Olfr", rownames(osn_rna[["RNA"]]@counts))]
mat <- osn_rna[["RNA"]]@counts[terms, cells]
cells_rna <- cells[which(colSums(mat > 0) == 1)]
cells <- colnames(osn_sct_raw[["trans"]])
terms <- unique(osn_sct_raw[["features"]]$Symbol[grep("^Olfr", osn_sct_raw[["features"]]$Symbol)])
mat <- matrix(0, nrow=length(terms), ncol=length(cells), dimnames=list(terms, cells))
for (i in 1:length(terms)) mat[i,] <- colSums(osn_sct_raw[["trans"]][which(osn_sct_raw[["features"]]$Symbol == terms[i]), cells, drop=F])
cells_nano <- cells[which(colSums(mat > 0) == 1)]
cells <- intersect(cells_rna, cells_nano)

ids <- which(rowSums(osn_sct_raw[["trans"]][, cells]) > 0)
terms <- table(osn_sct_raw[["features"]]$Symbol[intersect(ids, grep("^Olfr", osn_sct_raw[["features"]]$Symbol))])
rec_ba <- data.frame(table(terms))
rec_ba
sum(rec_ba$Freq)

ids <- which(rowSums(osn_sct_raw[["trans"]][, cells]) > 0)
terms <- table(osn_sct_raw[["features"]]$Symbol[intersect(ids, grep("^Olfr", osn_sct_raw[["features"]]$Symbol))])
terms <- names(terms[which(terms > 1)])
res <- data.frame(Term=terms, Type=0, Num=0, Cells=0, Sed=0, Trans="", Count="", Anno="")
for (i in 1:length(terms))
{
	ids <- which(osn_sct_raw[["features"]]$Symbol == terms[i])
	rs <- rowSums(osn_sct_raw[["trans"]][ids, cells])
	ids <- ids[which(rs > 0)]
	res$Count[i] <- paste(rowSums(osn_sct_raw[["trans"]][ids, cells]), collapse=";")
	res$Type[i] <- max(as.numeric(names(table(colSums(osn_sct_raw[["trans"]][ids, cells] > 0)))))
	res$Num[i] <- length(ids)
	res$Cells[i] <- length(which(colSums(osn_sct_raw[["trans"]][ids, cells]) > 0))
	rmax <- rowMaxs(osn_sct_raw[["trans"]][ids, cells])
	rmax <- rmax[order(rmax, decreasing=T)]
	res$Sed[i] <- rmax[2]
	if (length(ids) == 2) res$Trans[i] <- paste0(names(rmax)[1], ";", names(rmax)[2]) else res$Trans[i] <- paste0(names(rmax)[1], ";", names(rmax)[2], ";", names(rmax)[3])
}
#res <- res[which(res$Type > 0 & res$Cells > 1 & res$Sed > 1),]
res <- res[which(res$Type > 0),]
write.csv(res, "osn_or_rec.csv", row.names=F, quote=F)














cells <- colnames(osn_sct_raw[["trans"]])[which(osn_rna$cell.subtype_fix == "Immature" | osn_rna$cell.subtype_fix == "Mature")]
cells <- colnames(osn_sct_raw[["trans"]])
terms <- unique(osn_sct_raw[["features"]]$Symbol[grep("^Olfr", osn_sct_raw[["features"]]$Symbol)])
mat <- matrix(0, nrow=length(terms), ncol=length(cells), dimnames=list(terms, cells))
for (i in 1:length(terms)) mat[i,] <- colSums(osn_sct_raw[["trans"]][which(osn_sct_raw[["features"]]$Symbol == terms[i]), cells, drop=F])
cells <- cells[which(colSums(mat > 1) == 1)]
write.csv(data.frame(BC=cells, Type="OSN"), "osn_olfr_bc.csv", row.names=F, quote=F)









#cells <- colnames(osn_sct_raw[["trans"]])[which(osn_rna$cell.subtype_fix == "Immature" | osn_rna$cell.subtype_fix == "Mature")]
cells <- colnames(osn_sct_raw[["trans"]])
terms <- unique(osn_sct_raw[["features"]]$Symbol[grep("^Olfr", osn_sct_raw[["features"]]$Symbol)])
mat <- matrix(0, nrow=length(terms), ncol=length(cells), dimnames=list(terms, cells))
for (i in 1:length(terms)) mat[i,] <- colSums(osn_sct_raw[["trans"]][which(osn_sct_raw[["features"]]$Symbol == terms[i]), cells, drop=F])
cells <- cells[which(colSums(mat > 1) == 1)]
write.csv(data.frame(BC=cells, Type="OSN"), "osn_olfr_bc.csv", row.names=F, quote=F)

ids <- which(rowSums(osn_sct_raw[["trans"]][, cells] > 1) > 0)
terms <- osn_sct_raw[["features"]][intersect(ids, grep("^Olfr", osn_sct_raw[["features"]]$Symbol)),]
terms <- data.frame(Term=terms$Symbol, Trans=terms$ID, Type="", Group="")
ids <- c()
for (t in unique(terms$Term))
{
	i <- which(terms$Term == t)
	if (length(i) > 1) ids <- c(ids, i)
}
terms <- terms[ids,]
#res_bb <- read.delim("or_iso_type.tsv")
#terms$Type[match(res_bb$Trans, terms$Trans)] <- res_bb$Type
#terms$Group[match(res_bb$Trans, terms$Trans)] <- res_bb$Group
#write.table(terms, "or_iso_type_all.tsv", row.names=F, quote=F, sep="\t")
terms <- read.delim("or_iso_type.tsv")




ids <- intersect(grep("NOVEL_", osn_sct_raw[["features"]]$Type), grep("^Olfr", osn_sct_raw[["features"]]$Symbol))


ids <- which(rowSums(osn_sct_raw[["trans"]][, cells] > 1) > 0)
terms <- table(osn_sct_raw[["features"]]$Symbol[intersect(ids, grep("^Olfr", osn_sct_raw[["features"]]$Symbol))])
rec_ba <- data.frame(table(terms))

pba <- ggplot(rec_ba, aes(x=terms, y=Freq))+
	geom_bar(stat="identity", width=0.8, fill=col_list[4])+
	scale_y_continuous(limits=c(0, 490), expand=c(0, 0))+
	geom_text(aes(x=terms, y=Freq+1, label=Freq), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Number of isoforms in OR", y="Count")+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())
ggsave(plot=pba, width=3, height=4, dpi=200, filename="or_type_test3.png", limitsize=F)

terms <- table(osn_sct_raw[["features"]]$Symbol[grep("^Olfr", osn_sct_raw[["features"]]$Symbol)])
terms <- names(terms[which(terms > 1)])
res <- data.frame(Term=terms, Type=0, Num=0, Cells=0, Sed=0, Trans="", Count="", Anno="")
for (i in 1:length(terms))
{
	ids <- names(which(rowSums(osn_sct_raw[["trans"]][which(osn_sct_raw[["features"]]$Symbol == terms[i]), cells] > 1) > 0))
	if (length(ids) < 2) next
	res$Count[i] <- paste(rowSums(osn_sct_raw[["trans"]][ids, cells]), collapse=";")
	res$Type[i] <- max(as.numeric(names(table(colSums(osn_sct_raw[["trans"]][ids, cells] > 1)))))
	res$Num[i] <- length(ids)
	res$Cells[i] <- length(which(colSums(osn_sct_raw[["trans"]][ids, cells] > 1) > 0))
	rmax <- rowMaxs(osn_sct_raw[["trans"]][ids, cells])
	rmax <- rmax[order(rmax, decreasing=T)]
	res$Sed[i] <- rmax[2]
	if (length(ids) == 2) res$Trans[i] <- paste0(names(rmax)[1], ";", names(rmax)[2]) else res$Trans[i] <- paste0(names(rmax)[1], ";", names(rmax)[2], ";", names(rmax)[3])
}
res <- res[which(res$Type > 0 & res$Cells > 1 & res$Sed > 1),]
#write.csv(res, "or_iso_test.csv", quote=F, row.names=F)

res <- read.csv("or_iso_test.csv", h=T)
rec <- data.frame(table(res$Type))
rec$Var1 <- factor(rec$Var1, levels=rec$Var1)
pca <- ggplot(rec, aes(x=Var1, y=Freq))+
	geom_bar(stat="identity", width=0.6, fill=col_list[4])+
	scale_y_continuous(limits=c(0, 23), expand=c(0, 0))+
	geom_text(aes(x=Var1, y=Freq+0.2, label=Freq), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Number of isoforms in cell", y="Count")+
	#scale_x_continuous(limits=c(0, max(term_total$dist)), breaks=seq(0, max(term_total$dist), 2))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())
ggsave(plot=pca, width=3, height=4, dpi=200, filename="or_type_test2.png", limitsize=F)


res <- data.frame(Type=c("same_CDS", "same_start", "diff_start"), Count=c(3, 39, 23))
res$Type <- factor(res$Type, levels=c("same_CDS", "same_start", "diff_start"))
ggsave(plot=ggplot(res, aes(x=Type, y=Count))+
	geom_bar(stat="identity", width=0.6, fill=col_list[4])+
	scale_y_continuous(limits=c(0, 43), expand=c(0, 0))+
	geom_text(aes(x=Type, y=Count+0.5, label=Count), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Type of OR isoforms", y="Count")+
	#scale_x_continuous(limits=c(0, max(term_total$dist)), breaks=seq(0, max(term_total$dist), 2))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()), 
	width=3, height=4, dpi=200, filename="or_type_test4.png", limitsize=F)


res_bb <- read.delim("or_iso_type_all.tsv")
res_bb <- res_bb[which(gsub(";.*", "", res_bb$Type) != "Full"),]
res_bb$Type <- factor(res_bb$Type, levels=rev(c("A3SS", "RI", "SE", "3UTR_RI", "5UTR_SE")))
res_bb <- data.frame(table(res_bb$Type))
pbb <- ggplot(res_bb, aes(x=Var1, y=Freq, fill=Var1))+
	geom_bar(stat="identity", width=0.6)+
	scale_fill_manual(values=col_list[c(2,5,6,4,8)])+
	scale_y_continuous(limits=c(0, 38), expand=c(0, 0))+
	geom_text(aes(x=Var1, y=Freq+0.5, label=Freq), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Type of OR isoforms", y="Count", fill="Type")+
	#scale_x_continuous(limits=c(0, max(term_total$dist)), breaks=seq(0, max(term_total$dist), 2))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())+tag_thm
pleb <- ggplot_gtable(ggplot_build(pbb))
plebb <- pleb$grobs[[which(sapply(pleb$grobs, function(x) x$name) == "guide-box")]]
pbb <- pbb+theme(legend.position="none")


























ics <- which(osn_sct_raw[["features"]]$Type == "NOVEL_IC")
ncs <- which(osn_sct_raw[["features"]]$Type == "NOVEL_NC")
rec_aa <- data.frame(Type=c("IC", "NC"), Number=c(length(ics), length(ncs)))

terms <- unique(osn_sct_raw[["features"]]$Gene[c(ics, ncs)])
rec <- data.frame()
for (x in terms)
{
	ids <- which(osn_sct_raw[["features"]]$Gene == x)
	cs <- colSums(osn_sct_raw[["trans"]][ids,, drop=F])
	rs <- rowSums(osn_sct_raw[["trans"]][ids,, drop=F])
	ids <- grep("NOVEL", names(rs))
	rec <- rbind(rec, data.frame(Term=names(rs)[ids], Count=rs[ids], Rate=rs[ids]*100/sum(rs), Cell=length(which(cs > 0))))
}

pba <- ggplot(rec, aes(x=log10(Count), y=Rate, size=Cell))+geom_point()+
	geom_pointdensity()+scale_color_viridis()+
	labs(title=NULL, x="Count (log10)", y="Propotion (%)", size="#Cells")+
	scale_x_continuous(expand=c(0, 0))+scale_y_continuous(expand=c(0, 0))+
	geom_vline(xintercept=1.3, color="red", linetype="dashed", linewidth=1)+
	geom_hline(yintercept=40, color="red", linetype="dashed", linewidth=1)+
	theme(axis.line=element_line(linetype=1, color='black'), panel.background=element_rect(0, linetype=0), 
	plot.title=element_text(size=title_size, hjust=0.5, colour="black", face="plain"), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, color="black"), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), legend.key.size=unit(15, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())+tag_thm
ggsave(plot=pba, width=6, height=4, dpi=200, "test_dist4.png", limitsize=F)

terms <- unique(osn_sct_raw[["features"]]$Symbol[match(rec$Term[which(rec$Count > 20 & rec$Rate > 40)], osn_sct_raw[["features"]]$ID)])
eid <- mapIds(org.Mm.eg.db, keys=terms, column="ENTREZID", keytype="SYMBOL", multiVals="first")
ego <- enrichGO(gene=eid, keyType="ENTREZID", OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH", readable=T)
ego <- ego@result[which((ego@result$p.adjust < 0.05 | ego@result$qvalue < 0.05) & ego@result$Count > 1),]
ego_id <- data.frame(ID=ego$ID, level=0)
id_filtered <- c()
for (i in 1:nrow(ego_id))
{
	if (!is.na(match(ego_id$ID[i], id_filtered))) next
	rec <- as.character(unlist(mget(ego_id$ID[i], GOBPPARENTS, ifnotfound=NA)))
	parents <- rec
	level <- 0
	while(length(parents) > 0)
	{
		parents_ori <- parents
		parents <- c()
		for (p in parents_ori)
		{
			if (p == "all" & ego_id$level[i] == 0) ego_id$level[i] <- level
			if (is.na(p) | p == "all") next
			rec <- c(rec, p)
			parents <- c(parents, as.character(unlist(mget(p, GOBPPARENTS, ifnotfound=NA))))
		}
		level <- level + 1
		parents <- unique(parents)
	}
	rec <- unique(rec)
	rec <- which(!is.na(match(ego_id$ID, rec)))
	if (length(rec) > 0) id_filtered <- c(id_filtered, ego_id$ID[rec])
}
ego_id$level[match(unique(id_filtered), ego_id$ID)] <- 0
ego_id <- ego_id[which(ego_id$level > 2),, drop=F]
ego <- ego[match(ego_id$ID, ego$ID),, drop=F]
ego$Level <- ego_id$level
ego <- ego[order(ego$pvalue),]
ego <- ego[order(ego$Count, decreasing=T),]
write.csv(ego, "test_go4.csv")

ego <- read.csv("test_go4.csv", r=1, h=T)
ego <- ego[1:10,]
ego$Rank <- factor(rev(1:nrow(ego)))
ego$Type <- "Novel isoforms with high proportion&counts"
rec_ngo <- ego
pngo <- ggplot(rec_ngo, aes(x=Count, y=Rank, color=pvalue, size=Count))+geom_point()+
	scale_color_gradient(low="#440255", high="#FFFFBF")+
	#scale_color_continuous(low=brewer.pal(11,"Spectral")[1], high=brewer.pal(11,"Spectral")[6], 
	#guide=guide_colorbar(reverse=TRUE))+
	scale_y_discrete(breaks=rec_ngo$Rank, labels=rec_ngo$Description, position="right")+
	facet_grid(Type~., scales="free_y", space="free_y", switch="y")+
	guides(color=guide_colorbar(order=1), size=guide_legend(order=0))+
	labs(title=NULL, x=NULL, y=NULL, color="p.val")+
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(colour="#800026", fill="#800026"), 
	strip.text=element_text(size=title_size*0.6, colour="white", face="bold"), 
	panel.spacing=unit(5, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	axis.line=element_blank(), panel.background=element_rect(fill='gray98'))+tag_thm
ggsave(plot=pngo, width=6, height=5, dpi=200, filename="test_go4.png", limitsize=F)


results <- parLapply(cl, terms, function(x) {
	ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
	rec_sub <- data.frame(Gene=x, RNA=0, Num=0, Count=0, Ctr=0, Pval=1, Cor=0, Diff=0, Trans="", Lv="")
	if (length(ids) < 2) return(rec_sub)
	mat <- as.matrix(osn_sct_tt[ids,])
	ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
	if (length(ts) < 2) return(rec_sub)
	rec_sub$RNA[1] <- sum(mat)/ncol(osn_rna)
	mat <- mat[ts,]
	rec_sub$Num[1] <- nrow(mat)
	rec_sub$Count[1] <- min(rowMaxs(mat))
	rec_sub$Ctr[1] <- min(rowSums(mat))*100/max(rowSums(mat))
	bk <- ceiling(bins/2)-1
	pp <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
	for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
	rk <- abs(mat[, 1] - mat[, bins])
	rec_sub$Diff[1] <- max(rk)
	rk <- rk[order(rk, decreasing=T)]
	rec_sub$Trans[1] <- paste(names(rk)[1:2], collapse=";")
	rec_sub$Lv[1] <- paste(round(mat[,5]+mat[,4]-mat[,2]-mat[,1])[names(rk)[1:2]], collapse=";")
	rec_sub$Pval[1] <- fisher.test(pp[names(rk)[1:2],])$p.value
	t <- rep(0, nrow(mat))
	for (k in 1:nrow(mat)) t[k] <- abs(summary(lm(mat[k,]~rec_val_tt))$adj.r.squared)
	rec_sub$Cor[1] <- max(t)
	rec_sub$Cor[1] <- max(abs(cor(t(mat[names(rk)[1:2],]), rec_val_tt))[, 1])
	return(rec_sub)
})
stopCluster(cl)
rec <- do.call("rbind", results)
rec$p.adj <- p.adjust(rec$Pval, "BH")




#grep("Olfr", osn_sct_raw[["features"]]$Symbol[c(ics, ncs)])
tt <- data.frame(Count=rowSums(osn_sct_raw[["trans"]][ics,]), Symbol=osn_sct_raw[["features"]]$Symbol[ics], Type=osn_sct_raw[["features"]]$Type[ics])
tt <- tt[order(tt[,1], decreasing=T),,]
tt <- tt[match(names(which(table(tt$Symbol) < 4)), tt$Symbol),]
genes <- tt$Symbol[1:300]
cell_num <- min(table(osn_rna$cell.subtype_fix))
col_types <- brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)]
names(col_types) <- types
genes = c("Acvr2b", "Birc5")
for (gene in genes)
{
	ids <- which(osn_sct_raw[["features"]]$Symbol == gene)
	if (length(ids) < 2) next
	#if (gene == "Nyap2") ids <- which(osn_sct_raw[["features"]]$Gene == gene & osn_sct_raw[["features"]]$ID != "ENSMUST00000068275")
	rec_tt <- matrix(0, nrow=length(ids), ncol=5, dimnames=list(osn_sct_raw[["features"]]$ID[ids], types))
	for (i in 1:5)
	{
		#cells <- sample(which(osn_rna$cell.subtype_fix == types[i]), cell_num)
		cells <- which(osn_rna$cell.subtype_fix == types[i])
		rec_tt[, i] <- rowSums(osn_sct_raw[["trans"]][ids, cells])*cell_num/length(cells)
	}
	ids_sub <- which(rowMaxs(rec_tt) > 1)
	if (length(ids_sub) < 2) next
	rec_tt <- rec_tt[ids_sub,]
	#for (i in 1:5) rec_tt[, i] <- rec_tt[, i]*100/max(1, sum(rec_tt[, i]))
	rec_ttt <- data.frame()
	for (i in which(colSums(rec_tt) > 0)) rec_ttt <- rbind(rec_ttt, data.frame(Type=types[i], 
		Isoform=rownames(rec_tt), Count=rec_tt[, i]))
	rec_ttt$Type <- factor(rec_ttt$Type, levels=types)
	rec_ttt$Isoform <- factor(rec_ttt$Isoform, levels=rev(rownames(rec_tt)))
	ppa <- ggplot(rec_ttt, aes(x=Type, y=Count, stratum=Isoform, alluvium=Isoform, fill=Isoform))+
		labs(title=NULL, x=NULL, y="Expr. (scaled)")+
		#geom_flow(alpha=1)+geom_stratum(alpha=1, color="white", linewidth=0)+
		geom_flow(alpha=0.8)+geom_stratum(alpha=1, color="white")+guides(fill="none")+
		scale_y_continuous(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
		scale_fill_manual(values=col_list[c(1:(nrow(rec_tt)))+5])+
		theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,colour="black"), 
		axis.title=element_text(size=title_size, colour="black"), axis.text.y=element_text(size=text_size, colour="black"), 
		axis.text.x=element_blank(), axis.ticks.x=element_blank())+tag_thm
	seq_id <- osn_sct_raw[["features"]]$Chr[ids[1]]
	strand <- osn_sct_raw[["features"]]$Strand[ids[1]]
	gene_id <- osn_sct_raw[["features"]]$Gene[ids[1]]
	info_exons <- data.frame()
	for (t in rownames(rec_tt))
	{
		end <- 0
		exons <- unlist(strsplit(paste0(osn_sct_raw[["features"]]$Exon[which(osn_sct_raw[["features"]]$ID == gsub("-", "_", t))], ";"), split=";"))
		for (e in unlist(strsplit(exons[1], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="UTR", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
		for (e in unlist(strsplit(exons[2], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="exon", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
		for (e in unlist(strsplit(exons[3], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="UTR", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
	}
	names(info_exons) <- c("seqnames", "start", "end", "strand", "type", "transcript_name", "gene_name")
	#for (t in 2:nrow(info_exons)) if (info_exons$start[t] == info_exons$end[t-1]) info_exons$start[t] <- info_exons$end[t-1] + 1
	info_exons <- info_exons[order(info_exons$type, decreasing=T),]
	info_exons <- info_exons[order(info_exons$end),]
	info_exons <- info_exons[order(info_exons$start),]
	info_exons <- info_exons[order(info_exons$transcript_name),]
	info_exons <- info_exons[which(info_exons$end - info_exons$start > 0),]
	info_ranges <- paste0(info_exons$seqnames[1], ": ", min(info_exons$start), " - ", max(info_exons$start))
	info_exons <- shorten_gaps(exons=info_exons, introns=to_intron(info_exons, "transcript_name"), group_var="transcript_name")
	info_introns <- info_exons[which(info_exons$type == "intron"),]
	info_utrs <- info_exons[which(info_exons$type == "UTR"),]
	info_exons <- info_exons[which(info_exons$type == "exon"),]
	info_introns$transcript_name <- factor(info_introns$transcript_name, levels=rownames(rec_tt))
	info_exons$transcript_name <- factor(info_exons$transcript_name, levels=rownames(rec_tt))
	info_utrs$transcript_name <- factor(info_utrs$transcript_name, levels=rownames(rec_tt))
	info_col <- rev(col_list[c(1:(nrow(rec_tt)))+5])
	names(info_col) <- rownames(rec_tt)
	ppb <- ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
		geom_intron(data=info_introns, aes(color=transcript_name, strand=strand), linewidth=0.5, arrow.min.intron.length=300)+
		geom_range(aes(fill=transcript_name, color=transcript_name), linewidth=0, height=0.3)+
		geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name, fill=transcript_name, 
		color=transcript_name), linewidth=0.5, height=0.1)+
		labs(title=NULL, x=NULL, y=NULL)+
		scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, labels=info_ranges)+
		geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), 
		x=0, hjust=0, vjust=0, nudge_y=0.3, size=0.3*text_size)+
		scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
		theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
		axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
		axis.text.x=element_text(size=text_size, colour="black"), legend.position="none")+tag_thm
	rec_ppp <- data.frame(Type=unique(rec_ttt$Type), Info="Line")
	rec_ppp$Type <- factor(rec_ppp$Type, levels=unique(rec_ttt$Type))
	ppp <- ggplot(rec_ppp, aes(y=1, x=Type, fill=Type))+geom_bar(stat="identity", width=0.32)+
		theme_minimal()+labs(title=NULL, y=NULL, x=NULL, fill=NULL)+
		scale_fill_manual(values=col_types[unique(rec_ttt$Type)])+
		scale_x_discrete(expand=c(0, 0))+scale_y_continuous(expand=c(0, 0))+
		theme(legend.position="none", panel.background=element_blank(), axis.title=element_blank(), 
		axis.text.y=element_blank(), axis.text.x=element_text(size=text_size, colour="black"), 
		axis.ticks=element_blank())+tag_thm
	ggsave(plot=wrap_elements(wrap_plots(A=ppa, B=ppb, C=ppp, design=c("AB\nCB"), heights=c(40, 1), widths=c(2, 3))+
		plot_annotation(gene, theme=theme(plot.title=element_text(size=title_size, colour="black", hjust=0.5)))), 
		width=6, height=3, dpi=200, filename=paste0(gene, "_test.png"), limitsize=F)
}



for (gene in genes)
{
	ids <- which(osn_sct_raw[["features"]]$Symbol == gene)
	if (length(ids) < 2) next
	#if (gene == "Nyap2") ids <- which(osn_sct_raw[["features"]]$Gene == gene & osn_sct_raw[["features"]]$ID != "ENSMUST00000068275")
	rec_tt <- matrix(0, nrow=length(ids), ncol=5, dimnames=list(osn_sct_raw[["features"]]$ID[ids], types))
	for (i in 1:5)
	{
		#cells <- sample(which(osn_rna$cell.subtype_fix == types[i]), cell_num)
		cells <- which(osn_rna$cell.subtype_fix == types[i])
		rec_tt[, i] <- rowSums(osn_sct_raw[["trans"]][ids, cells])*cell_num/length(cells)
	}
	ids_sub <- which(rowMaxs(rec_tt) > 1)
	if (length(ids_sub) < 2) next
	rec_tt <- rec_tt[ids_sub,]
	#for (i in 1:5) rec_tt[, i] <- rec_tt[, i]*100/max(1, sum(rec_tt[, i]))
	rec_ttt <- data.frame()
	for (i in which(colSums(rec_tt) > 0)) rec_ttt <- rbind(rec_ttt, data.frame(Type=types[i], 
		Isoform=rownames(rec_tt), Count=rec_tt[, i]))
	rec_ttt$Type <- factor(rec_ttt$Type, levels=types)
	rec_ttt$Isoform <- factor(rec_ttt$Isoform, levels=rev(rownames(rec_tt)))
	ppa <- ggplot(rec_ttt, aes(x=Type, y=Count, stratum=Isoform, alluvium=Isoform, fill=Isoform))+
		labs(title=NULL, x=NULL, y="Expr. (scaled)")+
		#geom_flow(alpha=1)+geom_stratum(alpha=1, color="white", linewidth=0)+
		geom_flow(alpha=0.8)+geom_stratum(alpha=1, color="white")+guides(fill="none")+
		scale_y_continuous(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
		scale_fill_manual(values=col_list[c(1:(nrow(rec_tt)))+5])+
		theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,colour="black"), 
		axis.title=element_text(size=title_size, colour="black"), axis.text.y=element_text(size=text_size, colour="black"), 
		axis.text.x=element_blank(), axis.ticks.x=element_blank())+tag_thm
	seq_id <- osn_sct_raw[["features"]]$Chr[ids[1]]
	strand <- osn_sct_raw[["features"]]$Strand[ids[1]]
	gene_id <- osn_sct_raw[["features"]]$Gene[ids[1]]
	info_exons <- data.frame()
	for (t in rownames(rec_tt))
	{
		end <- 0
		exons <- unlist(strsplit(paste0(osn_sct_raw[["features"]]$Exon[which(osn_sct_raw[["features"]]$ID == gsub("-", "_", t))], ";"), split=";"))
		for (e in unlist(strsplit(exons[1], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="UTR", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
		for (e in unlist(strsplit(exons[2], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="exon", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
		for (e in unlist(strsplit(exons[3], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="UTR", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
	}
	names(info_exons) <- c("seqnames", "start", "end", "strand", "type", "transcript_name", "gene_name")
	#for (t in 2:nrow(info_exons)) if (info_exons$start[t] == info_exons$end[t-1]) info_exons$start[t] <- info_exons$end[t-1] + 1
	info_exons <- info_exons[order(info_exons$type, decreasing=T),]
	info_exons <- info_exons[order(info_exons$end),]
	info_exons <- info_exons[order(info_exons$start),]
	info_exons <- info_exons[order(info_exons$transcript_name),]
	info_exons <- info_exons[which(info_exons$end - info_exons$start > 0),]
	info_ranges <- paste0(info_exons$seqnames[1], ": ", min(info_exons$start), " - ", max(info_exons$start))
	info_exons <- shorten_gaps(exons=info_exons, introns=to_intron(info_exons, "transcript_name"), group_var="transcript_name")
	info_introns <- info_exons[which(info_exons$type == "intron"),]
	info_utrs <- info_exons[which(info_exons$type == "UTR"),]
	info_exons <- info_exons[which(info_exons$type == "exon"),]
	info_introns$transcript_name <- factor(info_introns$transcript_name, levels=rownames(rec_tt))
	info_exons$transcript_name <- factor(info_exons$transcript_name, levels=rownames(rec_tt))
	info_utrs$transcript_name <- factor(info_utrs$transcript_name, levels=rownames(rec_tt))
	info_col <- rev(col_list[c(1:(nrow(rec_tt)))+5])
	names(info_col) <- rownames(rec_tt)
	ppb <- ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
		geom_range(aes(fill=transcript_name), linewidth=0, height=0.3, color="white")+
		geom_intron(data=info_introns, aes(color=transcript_name, strand=strand), 
		linewidth=1, arrow.min.intron.length=300)+
		#geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name), 
		#linewidth=0, height=0.3, fill="gray85", color="white")+
		labs(title=NULL, x=NULL, y=NULL)+
		scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, labels=info_ranges)+
		geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), 
		x=0, hjust=0, vjust=0, nudge_y=0.3, size=0.5*text_size)+
		scale_fill_manual(values=info_col)+
		scale_color_manual(values=info_col, drop=F)+
		theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
		axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
		axis.text.x=element_text(size=text_size, colour="black"), legend.position="none")+tag_thm
	ggsave(plot=wrap_elements(ppb+
		plot_annotation(gene, theme=theme(plot.title=element_text(size=title_size, colour="black", hjust=0.5)))), 
		width=6, height=4, dpi=200, filename=paste0(gene, "_test.png"), limitsize=F)
}


tt <- data.frame(rowSums(osn_sct_raw[["trans"]][ncs,]), Symbol=osn_sct_raw[["features"]]$Symbol[ncs], Type=osn_sct_raw[["features"]]$Type[ncs])
tt <- tt[order(tt[,1], decreasing=T),,]
tt <- tt[match(names(which(table(tt$Symbol) < 2)), tt$Symbol),]
genes <- tt$Symbol[1:50]
for (gene in genes)
{
	ids <- which(osn_sct_raw[["features"]]$Symbol == gene)
	if (length(ids) < 2) next
	#if (gene == "Nyap2") ids <- which(osn_sct_raw[["features"]]$Gene == gene & osn_sct_raw[["features"]]$ID != "ENSMUST00000068275")
	rec_tt <- matrix(0, nrow=length(ids), ncol=5, dimnames=list(osn_sct_raw[["features"]]$ID[ids], types))
	for (i in 1:5)
	{
		#cells <- sample(which(osn_rna$cell.subtype_fix == types[i]), cell_num)
		cells <- which(osn_rna$cell.subtype_fix == types[i])
		rec_tt[, i] <- rowSums(osn_sct_raw[["trans"]][ids, cells])*cell_num/length(cells)
	}
	ids_sub <- which(rowMaxs(rec_tt) > 1)
	if (length(ids_sub) < 2) next
	rec_tt <- rec_tt[ids_sub,]
	#for (i in 1:5) rec_tt[, i] <- rec_tt[, i]*100/max(1, sum(rec_tt[, i]))
	rec_ttt <- data.frame()
	for (i in which(colSums(rec_tt) > 0)) rec_ttt <- rbind(rec_ttt, data.frame(Type=types[i], 
		Isoform=rownames(rec_tt), Count=rec_tt[, i]))
	rec_ttt$Type <- factor(rec_ttt$Type, levels=types)
	rec_ttt$Isoform <- factor(rec_ttt$Isoform, levels=rev(rownames(rec_tt)))
	ppa <- ggplot(rec_ttt, aes(x=Type, y=Count, stratum=Isoform, alluvium=Isoform, fill=Isoform))+
		labs(title=NULL, x=NULL, y="Expr. (scaled)")+
		#geom_flow(alpha=1)+geom_stratum(alpha=1, color="white", linewidth=0)+
		geom_flow(alpha=0.8)+geom_stratum(alpha=1, color="white")+guides(fill="none")+
		scale_y_continuous(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
		scale_fill_manual(values=col_list[c(1:(nrow(rec_tt)))+5])+
		theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,colour="black"), 
		axis.title=element_text(size=title_size, colour="black"), axis.text.y=element_text(size=text_size, colour="black"), 
		axis.text.x=element_blank(), axis.ticks.x=element_blank())+tag_thm
	seq_id <- osn_sct_raw[["features"]]$Chr[ids[1]]
	strand <- osn_sct_raw[["features"]]$Strand[ids[1]]
	gene_id <- osn_sct_raw[["features"]]$Gene[ids[1]]
	info_exons <- data.frame()
	for (t in rownames(rec_tt))
	{
		end <- 0
		exons <- unlist(strsplit(paste0(osn_sct_raw[["features"]]$Exon[which(osn_sct_raw[["features"]]$ID == gsub("-", "_", t))], ";"), split=";"))
		for (e in unlist(strsplit(exons[1], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="UTR", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
		for (e in unlist(strsplit(exons[2], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="exon", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
		for (e in unlist(strsplit(exons[3], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seq_id=seq_id, start=e[1], end=e[2], strand=strand, type="UTR", transcript_id=t, gene_id=gene_id))
			end = e[2]
		}
	}
	names(info_exons) <- c("seqnames", "start", "end", "strand", "type", "transcript_name", "gene_name")
	#for (t in 2:nrow(info_exons)) if (info_exons$start[t] == info_exons$end[t-1]) info_exons$start[t] <- info_exons$end[t-1] + 1
	info_exons <- info_exons[order(info_exons$type, decreasing=T),]
	info_exons <- info_exons[order(info_exons$end),]
	info_exons <- info_exons[order(info_exons$start),]
	info_exons <- info_exons[order(info_exons$transcript_name),]
	info_exons <- info_exons[which(info_exons$end - info_exons$start > 0),]
	info_ranges <- paste0(info_exons$seqnames[1], ": ", min(info_exons$start), " - ", max(info_exons$start))
	info_exons <- shorten_gaps(exons=info_exons, introns=to_intron(info_exons, "transcript_name"), group_var="transcript_name")
	info_introns <- info_exons[which(info_exons$type == "intron"),]
	info_utrs <- info_exons[which(info_exons$type == "UTR"),]
	info_exons <- info_exons[which(info_exons$type == "exon"),]
	info_introns$transcript_name <- factor(info_introns$transcript_name, levels=rownames(rec_tt))
	info_exons$transcript_name <- factor(info_exons$transcript_name, levels=rownames(rec_tt))
	info_utrs$transcript_name <- factor(info_utrs$transcript_name, levels=rownames(rec_tt))
	info_col <- rev(col_list[c(1:(nrow(rec_tt)))+5])
	names(info_col) <- rownames(rec_tt)
	ppb <- ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
		geom_range(aes(fill=transcript_name), linewidth=0, height=0.3, color="white")+
		geom_intron(data=info_introns, aes(color=transcript_name, strand=strand), 
		linewidth=1, arrow.min.intron.length=300)+
		#geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name), 
		#linewidth=0, height=0.3, fill="gray85", color="white")+
		labs(title=NULL, x=NULL, y=NULL)+
		scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, labels=info_ranges)+
		geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), 
		x=0, hjust=0, vjust=0, nudge_y=0.3, size=0.5*text_size)+
		scale_fill_manual(values=info_col)+
		scale_color_manual(values=info_col, drop=F)+
		theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
		axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
		axis.text.x=element_text(size=text_size, colour="black"), legend.position="none")+tag_thm
	ggsave(plot=wrap_elements(ppb+
		plot_annotation(gene, theme=theme(plot.title=element_text(size=title_size, colour="black", hjust=0.5)))), 
		width=6, height=4, dpi=200, filename=paste0(gene, "_test.png"), limitsize=F)
}









rec_ab <- data.frame()
for (type in types)
{
	cs <- which(osn_rna$cell.subtype_fix == type)
	rec_ab <- rbind(rec_ab, data.frame(Type=type, Number=length(which(rowSums(osn_sct_raw[["trans"]][c(ics, ncs), cs] > 1) > 1))))
}
rec_ab$Type <- factor(rec_ab$Type, levels=types)
rec_ab$Group <- "In catalog"

rec_ac <- data.frame()
for (type in types)
{
	cs <- which(osn_rna$cell.subtype_fix == type)
	rec_ac <- rbind(rec_ac, data.frame(Type=type, Number=length(which(rowSums(osn_sct_raw[["trans"]][ncs, cs] > 1) > 1))))
}
rec_ac$Type <- factor(rec_ac$Type, levels=types)
rec_ac$Group <- "Novel junctions"
rec_ab <- rbind(rec_ab, rec_ac)

paa <- ggplot(rec_aa, aes(x=Type, y=Number, fill=Type))+
	geom_bar(stat="identity", width=0.6, position=position_dodge(0.8))+
	labs(title=NULL, x=NULL, y="Count")+
	scale_fill_manual(values=col_list[c(1,3)])+
	scale_x_discrete(labels=c("In catalog", "Novel junctions"))+
	scale_y_continuous(limits=c(0, 1180), expand=c(0, 0))+guides(fill="none")+
	geom_text(aes(x=Type, y=Number+2, label=Number), size=0.5*text_size, vjust=0)+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=title_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())
ggsave(plot=paa, width=4, height=4, dpi=200, filename="fig_04_A.png", limitsize=F)

pab <- ggplot(rec_ab, aes(x=Type, y=Number, fill=Group))+geom_bar(stat="identity", width=0.6)+
	labs(title=NULL, x=NULL, y="Number of novel isoforms", fill="Type")+
	scale_fill_manual(values=col_list[c(2, 6)], drop=F)+
	scale_y_continuous(expand=c(0, 0))+
	#geom_text(aes(x=Type, y=Number+2, label=Number), size=0.5*text_size, vjust=0)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"))+tag_thm
ggsave(plot=pab, width=5, height=4, dpi=200, filename="fig_04_B.png", limitsize=F)


#cells <- colnames(osn_sct_raw[["trans"]])[which(osn_rna$cell.subtype_fix == "Immature" | osn_rna$cell.subtype_fix == "Mature")]
cells <- colnames(osn_sct_raw[["trans"]])
terms <- unique(osn_sct_raw[["features"]]$Symbol[grep("^Olfr", osn_sct_raw[["features"]]$Symbol)])
mat <- matrix(0, nrow=length(terms), ncol=length(cells), dimnames=list(terms, cells))
for (i in 1:length(terms)) mat[i,] <- colSums(osn_sct_raw[["trans"]][which(osn_sct_raw[["features"]]$Symbol == terms[i]), cells, drop=F])
cells <- cells[which(colSums(mat > 1) == 1)]
write.csv(data.frame(BC=cells, Type="OSN"), "osn_olfr_bc.csv", row.names=F, quote=F)

ids <- which(rowSums(osn_sct_raw[["trans"]][, cells] > 1) > 0)
terms <- osn_sct_raw[["features"]][intersect(ids, grep("^Olfr", osn_sct_raw[["features"]]$Symbol)),]
terms <- data.frame(Term=terms$Symbol, Trans=terms$ID, Type="", Group="")
ids <- c()
for (t in unique(terms$Term))
{
	i <- which(terms$Term == t)
	if (length(i) > 1) ids <- c(ids, i)
}
terms <- terms[ids,]
#res_bb <- read.delim("or_iso_type.tsv")
#terms$Type[match(res_bb$Trans, terms$Trans)] <- res_bb$Type
#terms$Group[match(res_bb$Trans, terms$Trans)] <- res_bb$Group
#write.table(terms, "or_iso_type_all.tsv", row.names=F, quote=F, sep="\t")
terms <- read.delim("or_iso_type.tsv")




ids <- intersect(grep("NOVEL_", osn_sct_raw[["features"]]$Type), grep("^Olfr", osn_sct_raw[["features"]]$Symbol))


ids <- which(rowSums(osn_sct_raw[["trans"]][, cells] > 1) > 0)
terms <- table(osn_sct_raw[["features"]]$Symbol[intersect(ids, grep("^Olfr", osn_sct_raw[["features"]]$Symbol))])
rec_ba <- data.frame(table(terms))

pba <- ggplot(rec_ba, aes(x=terms, y=Freq))+
	geom_bar(stat="identity", width=0.8, fill=col_list[4])+
	scale_y_continuous(limits=c(0, 490), expand=c(0, 0))+
	geom_text(aes(x=terms, y=Freq+1, label=Freq), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Number of isoforms in OR", y="Count")+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())
ggsave(plot=pba, width=3, height=4, dpi=200, filename="or_type_test3.png", limitsize=F)

terms <- table(osn_sct_raw[["features"]]$Symbol[grep("^Olfr", osn_sct_raw[["features"]]$Symbol)])
terms <- names(terms[which(terms > 1)])
res <- data.frame(Term=terms, Type=0, Num=0, Cells=0, Sed=0, Trans="", Count="", Anno="")
for (i in 1:length(terms))
{
	ids <- names(which(rowSums(osn_sct_raw[["trans"]][which(osn_sct_raw[["features"]]$Symbol == terms[i]), cells] > 1) > 0))
	if (length(ids) < 2) next
	res$Count[i] <- paste(rowSums(osn_sct_raw[["trans"]][ids, cells]), collapse=";")
	res$Type[i] <- max(as.numeric(names(table(colSums(osn_sct_raw[["trans"]][ids, cells] > 1)))))
	res$Num[i] <- length(ids)
	res$Cells[i] <- length(which(colSums(osn_sct_raw[["trans"]][ids, cells] > 1) > 0))
	rmax <- rowMaxs(osn_sct_raw[["trans"]][ids, cells])
	rmax <- rmax[order(rmax, decreasing=T)]
	res$Sed[i] <- rmax[2]
	if (length(ids) == 2) res$Trans[i] <- paste0(names(rmax)[1], ";", names(rmax)[2]) else res$Trans[i] <- paste0(names(rmax)[1], ";", names(rmax)[2], ";", names(rmax)[3])
}
res <- res[which(res$Type > 0 & res$Cells > 1 & res$Sed > 1),]
#write.csv(res, "or_iso_test.csv", quote=F, row.names=F)

res <- read.csv("or_iso_test.csv", h=T)
rec <- data.frame(table(res$Type))
rec$Var1 <- factor(rec$Var1, levels=rec$Var1)
pca <- ggplot(rec, aes(x=Var1, y=Freq))+
	geom_bar(stat="identity", width=0.6, fill=col_list[4])+
	scale_y_continuous(limits=c(0, 23), expand=c(0, 0))+
	geom_text(aes(x=Var1, y=Freq+0.2, label=Freq), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Number of isoforms in cell", y="Count")+
	#scale_x_continuous(limits=c(0, max(term_total$dist)), breaks=seq(0, max(term_total$dist), 2))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())
ggsave(plot=pca, width=3, height=4, dpi=200, filename="or_type_test2.png", limitsize=F)


res <- data.frame(Type=c("same_CDS", "same_start", "diff_start"), Count=c(3, 39, 23))
res$Type <- factor(res$Type, levels=c("same_CDS", "same_start", "diff_start"))
ggsave(plot=ggplot(res, aes(x=Type, y=Count))+
	geom_bar(stat="identity", width=0.6, fill=col_list[4])+
	scale_y_continuous(limits=c(0, 43), expand=c(0, 0))+
	geom_text(aes(x=Type, y=Count+0.5, label=Count), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Type of OR isoforms", y="Count")+
	#scale_x_continuous(limits=c(0, max(term_total$dist)), breaks=seq(0, max(term_total$dist), 2))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()), 
	width=3, height=4, dpi=200, filename="or_type_test4.png", limitsize=F)


res_bb <- read.delim("or_iso_type_all.tsv")
res_bb <- res_bb[which(gsub(";.*", "", res_bb$Type) != "Full"),]
res_bb$Type <- factor(res_bb$Type, levels=rev(c("A3SS", "RI", "SE", "3UTR_RI", "5UTR_SE")))
res_bb <- data.frame(table(res_bb$Type))
pbb <- ggplot(res_bb, aes(x=Var1, y=Freq, fill=Var1))+
	geom_bar(stat="identity", width=0.6)+
	scale_fill_manual(values=col_list[c(2,5,6,4,8)])+
	scale_y_continuous(limits=c(0, 38), expand=c(0, 0))+
	geom_text(aes(x=Var1, y=Freq+0.5, label=Freq), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Type of OR isoforms", y="Count", fill="Type")+
	#scale_x_continuous(limits=c(0, max(term_total$dist)), breaks=seq(0, max(term_total$dist), 2))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())+tag_thm
pleb <- ggplot_gtable(ggplot_build(pbb))
plebb <- pleb$grobs[[which(sapply(pleb$grobs, function(x) x$name) == "guide-box")]]
pbb <- pbb+theme(legend.position="none")

data_raw <- read.delim("or_test_exam.tsv")
cols <- col_list[c(2,4,5,6,8)]
names(cols) <- c("A3SS", "A5SS", "RI", "SE", "UTR")
plsb <- lapply(c("A3SS", "RI", "SE"), function(type) {
	ids <- which(data_raw$Type == type)
	seq_id <- data_raw$Chr[ids[1]]
	strand <- data_raw$Strand[ids[1]]
	gene_id <- data_raw$Gene[ids[1]]
	info_exons <- data.frame()
	for (t in ids)
	{
		end <- 0
		exons <- unlist(strsplit(data_raw$Exon[t], split=";"))
		for (e in unlist(strsplit(exons[1], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seqnames=seq_id, start=e[1], end=e[2], strand=strand, type="UTR", transcript_name=data_raw$Trans[t], gene_name=gene_id))
			end = e[2]
		}
		for (e in unlist(strsplit(exons[2], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seqnames=seq_id, start=e[1], end=e[2], strand=strand, type="CDS", transcript_name=data_raw$Trans[t], gene_name=gene_id))
			end = e[2]
		}
		for (e in unlist(strsplit(exons[3], ",")))
		{
			e <- as.numeric(unlist(strsplit(e, "-")))
			if (e[1] == end) e[1] = e[1] + 1
			info_exons <- rbind(info_exons, data.frame(seqnames=seq_id, start=e[1], end=e[2], strand=strand, type="UTR", transcript_name=data_raw$Trans[t], gene_name=gene_id))
			end = e[2]
		}
	}
	info_exons <- info_exons[order(info_exons$type, decreasing=T),]
	info_exons <- info_exons[order(info_exons$end),]
	info_exons <- info_exons[order(info_exons$start),]
	info_exons <- info_exons[order(info_exons$transcript_name),]
	info_exons <- info_exons[which(info_exons$end - info_exons$start > 0),]
	info_ranges <- paste0(info_exons$seqnames[1], ": ", min(info_exons$start), " - ", max(info_exons$start))
	info_exons <- shorten_gaps(exons=info_exons, introns=to_intron(info_exons, "transcript_name"), group_var="transcript_name")
	info_introns <- info_exons[which(info_exons$type == "intron"),]
	info_utrs <- info_exons[which(info_exons$type == "UTR"),]
	info_cds <- info_exons[which(info_exons$type == "CDS"),]
	info_introns$transcript_name <- factor(info_introns$transcript_name, levels=data_raw$Trans[ids])
	info_utrs$transcript_name <- factor(info_utrs$transcript_name, levels=data_raw$Trans[ids])
	info_cds$transcript_name <- factor(info_cds$transcript_name, levels=data_raw$Trans[ids])
	return(ggplot(info_cds, aes(xstart=start, xend=end, y=transcript_name))+
		geom_intron(data=info_introns, aes(strand=strand), color=cols[type], linewidth=0.5, arrow.min.intron.length=200)+
		geom_range(fill=cols[type], color=cols[type], linewidth=0, height=0.3)+
		geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name), fill=cols[type], color=cols[type], linewidth=0.5, height=0.1)+
		labs(title=NULL, x=NULL, y=type)+
		#scale_x_continuous(breaks=(max(c(info_cds$end, info_utrs$end))+min(c(info_cds$start, info_utrs$start)))/2, labels=info_ranges)+
		geom_text(data=info_introns, aes(y=transcript_name, label=transcript_name), color=cols[type], x=0, hjust=0, vjust=0, nudge_y=0.3, size=0.3*text_size)+
		#scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
		theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
		axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(size=title_size, colour="black"), 
		legend.position="none")+tag_thm)
})
pttt <- wrap_plots(A=pbb, B=plsb[[1]], C=plsb[[2]], D=plsb[[3]], E=plebb, design=c("ABE\nACE\nADE"), widths=c(5,6,1.2))+tag_thm
ggsave(pttt, width=11, height=4, dpi=200, filename="test_or4e.png", limitsize=F)



#data_raw <- read.delim("or_iso_test_aa_cal.tsv")
#data_raw <- data_raw[grep("_NOVEL", data_raw$Term),]
#for (i in 1:nrow(data_raw)) data_raw$Term[i] <- paste(unlist(strsplit(data_raw$Term[i], split="_"))[5:6], collapse="_")
#res_bb <- read.delim("or_iso_type.tsv")
#data_raw$Group <- res_bb$Type[match(data_raw$Term, res_bb$Trans)]
#data_raw <- data_raw[which(data_raw$Group != "UTR"),]
#write.table(data_raw, "or_iso_test_aa_part.tsv", quote=F, sep="\t")
data_raw <- read.delim("or_iso_test_aa_part.tsv")
data_raw <- data_raw[grep("_NOVEL", data_raw$Term),]
data_df <- data.frame()
for (i in 3:9) data_df <- rbind(data_df, data.frame(Term=data_raw[, 1], Type=data_raw[, 2], Pos=colnames(data_raw)[i], Rate=data_raw[, i]))
ptd <- ggplot(data_df, aes(x=Pos, y=Term, fill=Rate))+geom_tile()+geom_text(aes(label=Rate), colour="black")+
	facet_grid(Type~., scales="free_y", space="free_y")+
	scale_fill_gradientn(colours=colorRampPalette(c("white", "royalblue4"))(100))+
	#scale_x_discrete(breaks=c(1, seq(0, pos_limit, 10)), position="top")+scale_y_discrete(expand=c(0, 0))+
	labs(title=NULL, x=NULL, y=NULL, fill="Match\nRate (%)")+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), 
	panel.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text.x=element_blank(), axis.text.y=element_text(size=10, colour="black"), 
	strip.text=element_text(size=title_size*0.8, colour="white", face="bold"), 
	panel.spacing=unit(5, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	legend.title=element_text(size=12, colour="black"), legend.text=element_text(size=10, colour="black"))

rec_de <- data.frame(Type=c("Domain1", "Domain2", "Domain3", "Domain4", "Domain5", "Domain6", "Domain7"), Info="Line")
pte <- ggplot(rec_de, aes(y="", x=Type, fill=Type))+geom_tile()+theme_minimal()+
	geom_text(aes(label=Type), size=title_size*0.3, colour="white")+
	labs(title=NULL, y=NULL, x=NULL, fill=NULL)+
	scale_fill_manual(values=c("#ffa24b", "#39ff7f", "#3ba8b4", "#ff73c5", "#ff2443", "#3e3eff", "#ffd23b"))+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), panel.spacing=unit(5, "pt"), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	strip.text.y=element_text(size=text_size, colour="black"), legend.position="none", 
	strip.background=element_rect(colour="white", fill="white"), axis.text=element_blank(), 
	axis.ticks=element_blank())+tag_thm

ptdd <- ggplot_gtable(ggplot_build(ptd))
stripr <- which(grepl('strip-r', ptdd$layout$name))
fills <- cols[c(1, 3, 4)]
k <- 1
for (i in stripr)
{
	j <- which(grepl('rect', ptdd$grobs[[i]]$grobs[[1]]$childrenOrder))
	ptdd$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
	k <- k+1
}

library("ggplotify")
ptdd <- as.ggplot(ptdd)

ggsave(plot=wrap_elements(wrap_plots(list(pte,ptd), ncol=1, heights=c(1,20)))+tag_thm, 
	width=12, height=5, dpi=200, filename="test.png", limitsize=F)


data_raw <- read.delim("../oe_full_length/or_iso_test_aa_test.tsv")
#data_raw <- data_raw[grep("_NOVEL", data_raw$Term),]
data_df <- data.frame()
for (i in 3:9) data_df <- rbind(data_df, data.frame(Term=data_raw[, 1], Type=data_raw[, 2], Pos=colnames(data_raw)[i], Rate=data_raw[, i]))
data_df$Type <- factor(data_df$Type, levels=c("SE", "A3SS", "RI"), labels=c("Sudo", "Shift", "Gap"))
data_df$Term <- factor(data_df$Term, levels=data_raw$Term)
data_df$Rate[which(data_df$Rate > 20)] <- 100
data_df$Rate[which(data_df$Rate <= 20)] <- 0
data_df$Rate <- factor(data_df$Rate, levels=c(0, 100), labels=c("Lost", "Exist"))
ptd <- ggplot(data_df, aes(x=Pos, y=Term, fill=Rate))+geom_tile()+#geom_text(aes(label=Rate), colour="black")+
	facet_grid(Type~., scales="free_y", space="free_y")+
	scale_fill_manual(values=c("white", "royalblue4"))+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	labs(title=NULL, x=NULL, y=NULL, fill="Domain")+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), 
	panel.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), 
	strip.text=element_text(size=title_size*0.8, colour="white", face="bold"), 
	panel.spacing=unit(5, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	legend.title=element_text(size=12, colour="black"), legend.text=element_text(size=10, colour="black"))

rec_de <- data.frame(Type=c("D1", "D2", "D3", "D4", "D5", "D6", "D7"), Info="Line")
pte <- ggplot(rec_de, aes(y="", x=Type, fill=Type))+geom_tile()+theme_minimal()+
	geom_text(aes(label=Type, color=Type), size=title_size*0.3)+
	labs(title=NULL, y=NULL, x=NULL, fill=NULL)+
	scale_fill_manual(values=c("#ffa24b", "#39ff7f", "#3ba8b4", "#ff73c5", "#ff2443", "#3e3eff", "#ffd23b"))+
	scale_color_manual(values=c("#ffa24b", "#39ff7f", "#3ba8b4", "#ff73c5", "#ff2443", "#3e3eff", "#ffd23b"))+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), panel.spacing=unit(5, "pt"), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	strip.text.y=element_text(size=text_size, colour="black"), legend.position="none", 
	strip.background=element_rect(colour="white", fill="white"), axis.text=element_blank(), 
	axis.ticks=element_blank())+tag_thm

ggsave(plot=wrap_elements(wrap_plots(list(pte,ptd), ncol=1, heights=c(1,20)))+tag_thm, 
	width=6, height=5, dpi=200, filename="test.png", limitsize=F)




res_da <- read.delim("or_iso_type.tsv")
res_da <- res_da[which(res_da$Type != "Full"),]
res_da$Group <- factor(res_da$Group, levels=c("Pseudo", "Shift", "Same"), labels=c("Pseudo", "Shift", "Gap"))
res_da <- data.frame(table(res_da$Group))
pda <- ggplot(res_da, aes(x=Var1, y=Freq, fill=Var1))+
	geom_bar(stat="identity", width=0.6)+
	scale_fill_manual(values=col_list[c(6,2,5)])+
	scale_y_continuous(limits=c(0, 20), expand=c(0, 0))+
	geom_text(aes(x=Var1, y=Freq+0.5, label=Freq), size=0.5*text_size, vjust=0)+
	labs(title=NULL, x="Type of OR isoforms", y="Count", fill="Type")+
	#scale_x_continuous(limits=c(0, max(term_total$dist)), breaks=seq(0, max(term_total$dist), 2))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())
pled <- ggplot_gtable(ggplot_build(pda))
pledd <- pled$grobs[[which(sapply(pled$grobs, function(x) x$name) == "guide-box")]]
pda <- pda+theme(legend.position="none")

pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(paa+tag_thm, pab+tag_thm, pngo+tag_thm), nrow=1, widths=c(3, 4, 1))+plot_annotation(tag_levels=list(c("A", "B", "C"))))+tag_thm, 
	wrap_elements(wrap_plots(list(pba+tag_thm, pttt+tag_thm), nrow=1, widths=c(3, 10))+plot_annotation(tag_levels=list(c("D", "E"))))+tag_thm, 
	wrap_elements(wrap_plots(list(pca+tag_thm, pab+tag_thm, pngo+tag_thm), nrow=1, heights=c(3, 3, 1))+plot_annotation(tag_levels=list(c("F", "G", ""))))+tag_thm, 
	wrap_elements(wrap_plots(list(pda+tag_thm, pab+tag_thm, pngo+tag_thm), nrow=1, heights=c(3, 3, 1))+plot_annotation(tag_levels=list(c("H", "I", ""))))+tag_thm), 
	#wrap_elements(wrap_plots(A=pca, B=pca, C=pca, D=pca, E=pca, 
	#design=c("ABC\nADE"), widths=c(1, 3, 3))+plot_annotation(tag_levels=list(c("E", "F", "", "G", ""))))+tag_thm), 
	ncol=1, heights=c(1, 1, 1, 1)), width=16, height=16, dpi=200, filename="oe_fl_fig_04.png", limitsize=F)


ggsave(plot=pda, width=4, height=4, dpi=200, filename="oe_fl_fig_04_part_h.png", limitsize=F)



type <- paste0(sample, "+BPlasma")
res_total <- data.frame()
bcs <- c()
terms <- c()
bcs <- c(bcs, colnames(sce)[which(sce$cell.temp == type)])
for (bc in bcs)
{
	ids <- which(bcrinfo$barcode == bc)
	if (length(ids) > 0) terms <- c(terms, ids)
}
bcrinfo <- bcrinfo[terms,]
res <- unique(bcrinfo$raw_clonotype_id)
res <- data.frame(Var1=res, Freq=bcrinfo$cc[match(res, bcrinfo$raw_clonotype_id)])
res$group <- "Unique"
res$group[which(res$Freq >= 2)] <- "2-4"
res$group[which(res$Freq >= 5)] <- "5-9"
res$group[which(res$Freq >= 10)] <- "≥10"
res$group <- factor(res$group, levels=c("Unique", "2-4", "5-9", "≥10"))
res <- data.frame(table(res$group))
res$Var1 <- factor(res$Var1, levels=c("Unique", "2-4", "5-9", "≥10"))


res_da <- read.delim("or_iso_type.tsv")
res_da <- res_da[which(res_da$Type != "Full"),]
res_da$Group <- factor(res_da$Group, levels=c("Sudo", "Shift", "Same"), labels=c("Sudo", "Shift", "Gap"))
res_da <- data.frame(table(res_da$Group))
pda <- ggplot(res_da, aes(x=4, y=Freq, fill=Var1))+geom_col()+
	labs(title="Type of OR isoforms", x=NULL, y=NULL, fill="Type")+
	geom_text(aes(label=Freq), position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+xlim(c(2.2, 4.5))+
	scale_fill_manual(values=col_list[c(6,2,5)])+
	theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
	axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), 
	axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
	legend.text=element_text(colour="black", size=12))+tag_thm
ggsave(plot=pda, width=4, height=4, dpi=200, filename="oe_fl_fig_04_part_h.png", limitsize=F)






pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(paa+tag_thm, pab+tag_thm, pngo+tag_thm), nrow=1, widths=c(3, 4, 1))+plot_annotation(tag_levels=list(c("A", "B", "C"))))+tag_thm, 
	wrap_elements(wrap_plots(list(pba+tag_thm, pbb+tag_thm, pblank), nrow=1, widths=c(3, 4, 6))+plot_annotation(tag_levels=list(c("D", "E", "F"))))+tag_thm, 
	wrap_elements(wrap_plots(list(pca+tag_thm, pab+tag_thm, pngo+tag_thm), nrow=1, heights=c(3, 3, 1))+plot_annotation(tag_levels=list(c("G", "H", ""))))+tag_thm), 
	#wrap_elements(wrap_plots(A=pca, B=pca, C=pca, D=pca, E=pca, 
	#design=c("ABC\nADE"), widths=c(1, 3, 3))+plot_annotation(tag_levels=list(c("E", "F", "", "G", ""))))+tag_thm), 
	ncol=1, heights=c(1, 1, 1)), width=16, height=12, dpi=200, filename="oe_fl_fig_04.png", limitsize=F)




info_introns$transcript_name <- factor(info_introns$transcript_name, levels=data_raw$Trans)
info_introns$type <- factor(data_raw$Type[match(info_introns$transcript_name, data_raw$Trans)], levels=types)
info_utrs$transcript_name <- factor(info_utrs$transcript_name, levels=data_raw$Trans)
info_utrs$type <- factor(data_raw$Type[match(info_utrs$transcript_name, data_raw$Trans)], levels=types)
info_cds$transcript_name <- factor(info_cds$transcript_name, levels=data_raw$Trans)
info_cds$type <- factor(data_raw$Type[match(info_cds$transcript_name, data_raw$Trans)], levels=types)

pttt <- ggplot(info_cds, aes(xstart=start, xend=end, y=transcript_name))+
	geom_intron(data=info_introns, aes(color=type, strand=strand), linewidth=0.5, arrow.min.intron.length=100)+
	geom_range(aes(fill=transcript_name, color=type), linewidth=0, height=0.3)+
	geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name, fill=type, color=type), linewidth=0.5, height=0.1)+
	labs(title=NULL, x=NULL, y=NULL)+
	facet_grid(type ~ ., scales="free", space="free_x")+
	#scale_x_continuous(breaks=(max(c(info_cds$end, info_utrs$end))+min(c(info_cds$start, info_utrs$start)))/2, labels=info_ranges)+
	#geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), x=0, hjust=0, vjust=0, nudge_y=0.3, size=0.3*text_size)+
	#scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
	theme(panel.border=element_blank(), panel.background=element_blank(), 
	axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black"), legend.position="none")+tag_thm
ggsave(pttt, width=5, height=6, dpi=200, filename="test.png", limitsize=F)



		theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
		axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
		axis.text.x=element_text(size=text_size, colour="black"), legend.position="none")+tag_thm










A3SS
Olfr1427	ENSMUST00000216069
chr19	12098006-12098705;12098705-12099637;12099637-12099665,12103132-12103334
Olfr1427	ENSMUSG00000067525_NOVEL00
chr19	12098680-12098702;12098702-12098969;12098969-12098979,12103132-12103238;
A5SS(+)
Olfr53	ENSMUST00000084456
chr7	140646452-140646588,140651895-140651981;140651981-140652916;140652916-140653085
Olfr53	ENSMUSG00000094819_NOVEL00
chr7	140646492-140646588;140652658-652743;652743-140653093
RI
Olfr1026	ENSMUST00000219615
chr2	85921558-85921638,85923240-85923270;85923270-85924190;85924190-85924623
Olfr1026	ENSMUSG00000042863_NOVEL00
chr2	85921607-85921638,85923240-85923270;85923270-85923390,85923940-85924190;85924190-85924618

SE
Olfr2	ENSMUST00000208147
chr7	106996259-107000878;107000878-107001858;107001858-107001981,107002463-107002736,107005519-107005561,107005921-107006058
Olfr2	ENSMUSG00000070417_NOVEL00
chr7	106995362-106995842,106995974-106996039;;107005519-107005561,107005921-107005949

UTR



	ids <- which(osn_sct_raw[["features"]]$Symbol == gene & osn_sct_raw[["features"]]$Type == "protein_coding")

osn_sct_raw[["features"]][which(osn_sct_raw[["features"]]$ID == "ENSMUST00000216069"),]
osn_sct_raw[["features"]][which(osn_sct_raw[["features"]]$ID == "ENSMUSG00000067525_NOVEL00"),]












