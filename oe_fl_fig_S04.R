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
library(factR)
library(BSgenome.Mmusculus.UCSC.mm10)
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

olfr_info_osn <- data.frame(read.delim("osn_iso_info_all.tsv", h=T), Sample="OSN")
olfr_info_oe <- data.frame(read.delim("oe_iso_info_all.tsv", h=T), Sample="OE")
terms <- intersect(olfr_info_osn$Exon, olfr_info_oe$Exon)
terms_osn <- setdiff(olfr_info_osn$Exon, olfr_info_oe$Exon)
terms_oe <- setdiff(olfr_info_oe$Exon, olfr_info_osn$Exon)
olfr_info <- olfr_info_osn[match(terms, olfr_info_osn$Exon),]
olfr_info$Sample <- "Both"
olfr_info <- rbind(olfr_info, rbind(olfr_info_osn[match(terms_osn, olfr_info_osn$Exon),], olfr_info_oe[match(terms_oe, olfr_info_oe$Exon),]))
olfr_info <- olfr_info[order(olfr_info$UTR),]
olfr_info <- olfr_info[order(olfr_info$TSS),]
olfr_info <- olfr_info[order(olfr_info$TypeE),]
olfr_info <- olfr_info[order(olfr_info$TypeC),]
olfr_info <- olfr_info[order(olfr_info$Gene),]

isoforms_osn <- data.frame(read.delim("osn_iso_info.tsv", h=T), Sample="OSN")
isoforms_oe <- data.frame(read.delim("oe_iso_info.tsv", h=T), Sample="OE")
terms <- intersect(isoforms_osn$Exon, isoforms_oe$Exon)
terms_osn <- setdiff(isoforms_osn$Exon, isoforms_oe$Exon)
terms_oe <- setdiff(isoforms_oe$Exon, isoforms_osn$Exon)
isoforms <- isoforms_osn[match(terms, isoforms_osn$Exon),]
isoforms$Count <- isoforms$Count + isoforms_oe$Count[match(terms, isoforms_oe$Exon)]
isoforms$Sample <- "Both"
isoforms <- rbind(isoforms, rbind(isoforms_osn[match(terms_osn, isoforms_osn$Exon),], isoforms_oe[match(terms_oe, isoforms_oe$Exon),]))
cells_nano <- rbind(data.frame(read.delim("osn_olfr_info.tsv", h=T), Sample="OSN"), data.frame(read.delim("oe_olfr_info.tsv", h=T), Sample="OE"))
iso_trans_test <- rbind(data.frame(read.delim("osn_iso_trans_test.tsv", h=T), Sample="OSN"), data.frame(read.delim("oe_iso_trans_test.tsv", h=T), Sample="OE"))
iso_trans_sub_osn <- data.frame(read.delim("osn_iso_trans_test.tsv", h=T), Sample="OSN")
iso_trans_sub_osn <- iso_trans_sub_osn[which(iso_trans_sub_osn$TypeC != "-"),]
trans_test_nmd <- as.data.frame(predictNMD(importGTF("/home/cbh/work/oe_fl_re/osn_iso_trans_test.gtf")))
iso_trans_sub_osn$stop_to_lastEJ <- trans_test_nmd$stop_to_lastEJ[match(iso_trans_sub_osn$ID, trans_test_nmd$transcript)]
iso_trans_sub_osn$is_NMD <- trans_test_nmd$is_NMD[match(iso_trans_sub_osn$ID, trans_test_nmd$transcript)]
motif_info <- read.delim("/mnt/md0/oe_full_length/output/OSNfulllength/iso_trans_test_aa_cal.tsv")
iso_trans_sub_osn <- cbind(iso_trans_sub_osn, motif_info[match(iso_trans_sub_osn$ID, gsub(".*[+-]_", "", motif_info$Term)), -1])
# nmdectectiveb/nmdetect osn_iso_trans_test.gtf osn_iso_trans_test
nmd_test <- read.delim("osn_iso_trans_test.tidtable.tsv", h=F)
iso_trans_sub_osn$NMDect <- gsub(" \\(.*$", "", nmd_test[match(iso_trans_sub_osn$ID, nmd_test[, 1]), 3])
iso_trans_sub_oe <- data.frame(read.delim("oe_iso_trans_test.tsv", h=T), Sample="OE")
iso_trans_sub_oe <- iso_trans_sub_oe[which(iso_trans_sub_oe$TypeC != "-"),]
trans_test_nmd <- as.data.frame(predictNMD(importGTF("/home/cbh/work/oe_fl_re/oe_iso_trans_test.gtf")))
iso_trans_sub_oe$stop_to_lastEJ <- trans_test_nmd$stop_to_lastEJ[match(iso_trans_sub_oe$ID, trans_test_nmd$transcript)]
iso_trans_sub_oe$is_NMD <- trans_test_nmd$is_NMD[match(iso_trans_sub_oe$ID, trans_test_nmd$transcript)]
motif_info <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/iso_trans_test_aa_cal.tsv")
iso_trans_sub_oe <- cbind(iso_trans_sub_oe, motif_info[match(iso_trans_sub_oe$ID, gsub(".*[+-]_", "", motif_info$Term)), -1])
# nmdectectiveb/nmdetect oe_iso_trans_test.gtf oe_iso_trans_test
nmd_test <- read.delim("oe_iso_trans_test.tidtable.tsv", h=F)
iso_trans_sub_oe$NMDect <- gsub(" \\(.*$", "", nmd_test[match(iso_trans_sub_oe$ID, nmd_test[, 1]), 3])
terms <- intersect(iso_trans_sub_osn$Exon, iso_trans_sub_oe$Exon)
terms_osn <- setdiff(iso_trans_sub_osn$Exon, iso_trans_sub_oe$Exon)
terms_oe <- setdiff(iso_trans_sub_oe$Exon, iso_trans_sub_osn$Exon)
iso_trans_sub <- iso_trans_sub_osn[match(terms, iso_trans_sub_osn$Exon),]
iso_trans_sub$Count <- iso_trans_sub$Count + iso_trans_sub_oe$Count[match(terms, iso_trans_sub_oe$Exon)]
iso_trans_sub$Sample <- "Both"
iso_trans_sub <- rbind(iso_trans_sub, rbind(iso_trans_sub_osn[match(terms_osn, iso_trans_sub_osn$Exon),], iso_trans_sub_oe[match(terms_oe, iso_trans_sub_oe$Exon),]))
iso_trans_sub$TypeC <- factor(iso_trans_sub$TypeC, levels=c("A5SS", "RI", "A5SS+RI"))
iso_trans_sub$ID <- paste(iso_trans_sub$ID, iso_trans_sub$Sample)

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

pad <- wrap_elements(ggdraw()+draw_image("figS04_b.png", scale=1))+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm
pae <- wrap_elements(ggdraw()+draw_image("figS04_c.png", scale=1))+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm

tag_thm3 <- theme(plot.tag=element_text(size=title_size, face="bold", color="black"), plot.margin=margin(-32,-5,-10,-5), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))
pad <- wrap_elements(ggdraw()+draw_image("figS04_b.png", scale=1)+tag_thm3)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm3

tag_thm4 <- theme(plot.tag=element_text(size=title_size, face="bold", color="black"), plot.margin=margin(-10,0,0,0), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))
paf <- wrap_elements(ggdraw()+draw_image("figS04_d.png", scale=1)+tag_thm4)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm
pxf <- wrap_elements(ggdraw()+draw_image("aaa.png", scale=1)+tag_thm4)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.4))+tag_thm

rxes <- data.frame()
terms <- c("Last exon", "Start-proximal", "Long exon", "50 nt rule", "Trigger NMD")
for (term in terms) rxes <- rbind(rxes, data.frame(Group=term, Count=length(which(iso_trans_sub$NMDect == term))))
rxes$Rate <- rxes$Count*100/sum(rxes$Count)
rxes$Group <- factor(rxes$Group, levels=rxes$Group)
pxes <- ggplot(rxes, aes(x=4, y=Rate, fill=Group))+geom_col()+
	labs(title=NULL, x=NULL, y=NULL, fill=NULL)+
	geom_text(aes(label=paste0(Rate, "%")), size=3, position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+scale_x_continuous(limits=c(2.5, 4.5), expand=c(0, 0)) +
	scale_fill_manual(values=col_list[c(4,3,5,6,7)])+guides(fill=guide_legend(ncol=1))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	plot.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), axis.title=element_blank(), legend.title=element_text(colour="black", size=title_size), 
	legend.text=element_text(colour="black", size=text_size), legend.position="bottom")+tag_thm
pxeb <- ggplot(iso_trans_sub, aes(x=ID, y=-log2(-stop_to_lastEJ)))+geom_bar(stat="identity", width=0.6, fill=col_list[4])+
	labs(title=NULL, x="OR Isoforms", y="stop_to_lastEJ (bp, log2)")+
	scale_y_continuous(limits=c(-13, 5.5), breaks=c(-12, -8, -4, 0, 5), expand=c(0, 0))+
	geom_hline(yintercept=5, linetype="dashed", color=col_list[3], linewidth=0.35)+
	geom_hline(yintercept=0, color="black", linewidth=0.35)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), axis.ticks.x=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black"), axis.text.x=element_blank(), 
	axis.title=element_text(size=title_size, colour="black"), 
	axis.line.y=element_line(linewidth=0.35, color="black"), axis.line.x=element_blank(), 
	axis.ticks.y=element_line(linewidth=0.35, color="black"), panel.spacing=unit(0, "pt"), 
	plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))+tag_thm

pxe <- wrap_elements(wrap_plots(list(wrap_elements(pxes)+tag_thm, wrap_elements(pxeb)), nrow=1, widths=c(1, 2))+
	tag_thm2)+tag_thm

gtf <- read.delim("gtf_info.tsv")
types <- c("HBC", "GBC", "INP", "iOSN", "mOSN")
osn_sct <- readRDS("osn_sct.rds")
osn_rna <- readRDS("osn_rna.rds")
osn_rna$nIsoform <- as.numeric(colSums(osn_sct[["trans"]] > 0))
osn_rna$nGene<- as.numeric(colSums(osn_sct[["genes"]] > 0))
data_raw <- H5File$new("/mnt/md0/oe_full_length/output/OEfulllength/osn_trans_ass.h5", mode="r")
sct <- list()
sparse.mat <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
	x=as.numeric(data_raw[["matrix/data"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
rownames(sparse.mat) <- data_raw[["matrix/genes"]][]
colnames(sparse.mat) <- data_raw[["matrix/barcodes"]][]
sparse.utr3 <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
	x=as.numeric(data_raw[["matrix/utr3"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
rownames(sparse.utr3) <- data_raw[["matrix/genes"]][]
colnames(sparse.utr3) <- data_raw[["matrix/barcodes"]][]
sparse.utr5 <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
	x=as.numeric(data_raw[["matrix/utr5"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
rownames(sparse.utr5) <- data_raw[["matrix/genes"]][]
colnames(sparse.utr5) <- data_raw[["matrix/barcodes"]][]
sct[["matrix"]] <- as.sparse(sparse.mat)
sct[["utr3"]] <- as.sparse(sparse.utr3)
sct[["utr5"]] <- as.sparse(sparse.utr5)
sct[["features"]] <- data.frame(Name=data_raw[["matrix/features/name"]][], ID=data_raw[["matrix/features/id"]][], 
	Chr=data_raw[["matrix/features/chr"]][], Strand=data_raw[["matrix/features/strand"]][], Gene=data_raw[["matrix/features/gene"]][], 
	Type=data_raw[["matrix/features/type"]][], UTR3=data_raw[["matrix/features/utr3"]][], UTR5=data_raw[["matrix/features/utr5"]][], 
	Body=data_raw[["matrix/features/body"]][], Exon=data_raw[["matrix/features/exon"]][])
sct[["features"]]$Count <- rowSums(sct[["matrix"]])
sct[["features"]]$Symbol <- gtf$gene_name[match(sct[["features"]]$Gene, gtf$gene_id)]
data_raw$close_all()
osn_sct_raw <- list()
osn_sct_raw[["trans"]] <- sct[["matrix"]]
osn_sct_raw[["utr3"]] <- sct[["utr3"]]
osn_sct_raw[["utr5"]] <- sct[["utr5"]]
osn_sct_raw[["features"]] <- sct[["features"]]
osn_sct_raw[["features"]]$ID <- sct[["features"]]$Name

olfr_info <- read.csv("/mnt/md0/oe_full_length/output/OEfulllength/olfr_info.csv", h=T)
olfr_info <- olfr_info[match(osn_sct_raw[["features"]]$ID[grep("^Olfr", osn_sct_raw[["features"]]$Symbol)], olfr_info$ID),]
olfr_info <- olfr_info[order(olfr_info$UTR),]
olfr_info <- olfr_info[order(olfr_info$TSS),]
olfr_info <- olfr_info[order(olfr_info$TypeE),]
olfr_info <- olfr_info[order(olfr_info$TypeC),]
olfr_info <- olfr_info[order(olfr_info$Gene),]
olfr_info$CSS <- "Same"
olfr_info$CSS[grep("A3SS", olfr_info$TypeC)] <- "Diff"
olfr_info$CSS[grep("SE", olfr_info$TypeC)] <- "Diff"
olfr_sct <- osn_sct_raw[["trans"]][match(olfr_info$ID, osn_sct_raw[["features"]]$ID),]
olfr_sct_info <- olfr_info[match(rownames(olfr_sct), olfr_info$ID),]
olfr_sct_info$Count <- rowSums(olfr_sct)

rec <- data.frame(X=osn_rna@reductions$umap@cell.embeddings[,1], 
	Y=osn_rna@reductions$umap@cell.embeddings[,2], Type=osn_rna$cell.subtype_fix, 
	Counts=colSums(osn_rna[["RNA"]]$counts[grep("^Olfr", rownames(osn_rna[["RNA"]]$counts)),] > 1))
rec$Counts[which(rec$Counts == 0)] <- "0"
rec$Counts[which(rec$Counts == 1)] <- "1"
rec$Counts[which(rec$Counts == 2)] <- "2"
rec$Counts[which(rec$Counts == 3)] <- "3"
rec$Counts[which(rec$Counts > 3)] <- "≥4"
rec$Counts <- factor(rec$Counts, levels=c("0", "1", "2", "3", "≥4"))
res <- data.frame()
for (type in types[3:5])
{
	ids <- which(rec$Type == type & rec$Counts != "0")
	res <- rbind(res, data.frame(Type=type, Total=length(ids), table(rec$Counts[ids])))
}
res <- res[which(res$Var1 != "0"),]
res$Rate <- res$Freq*100/res$Total
res$Type <- factor(res$Type, levels=types[3:5])
res$Var1 <- factor(as.character(res$Var1), levels=c("1", "2", "3", "≥4"))
res$Group <- "A"
res_tta <- res
ptta <- ggplot(res_tta, aes(x=Type, y=Rate, fill=Var1))+geom_bar(stat="identity", position=position_dodge(0.7), width=0.6)+
	labs(title="# OR genes\n(by illumina)", x=NULL, y="Percentage (%)", fill="Num.")+
	scale_fill_manual(values=brewer.pal(9,"YlGnBu")[c(9,7,5,3)])+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_y_continuous(breaks=c(25, 50, 75), expand=c(0, 0))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), axis.ticks.x=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title.y=element_text(size=title_size, colour="black"), 
	axis.title.x=element_text(size=title_size, colour="black", hjust=1.3), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks.y=element_line(linewidth=0.35, color="black"))+tag_thm

ids <- match(colnames(olfr_sct), colnames(osn_rna))
rec <- data.frame(Cell=colnames(olfr_sct), X=osn_rna@reductions$umap@cell.embeddings[ids, 1], 
	Y=osn_rna@reductions$umap@cell.embeddings[ids, 2], Type=osn_rna$cell.subtype_fix[ids], Gene=0, Isoform=0)
for (i in 1:ncol(olfr_sct))
{
	terms <- olfr_sct_info[which(olfr_sct[, i] > 0),, drop=F]
	rec$Gene[i] <- length(table(terms$Gene))
	rec$Isoform[i] <- length(table(terms$ID))
}
rec$Counts <- rec$Gene
rec$Counts[which(rec$Counts == 0)] <- "0"
rec$Counts[which(rec$Counts == 1)] <- "1"
rec$Counts[which(rec$Counts == 2)] <- "2"
rec$Counts[which(rec$Counts == 3)] <- "3"
rec$Counts[which(rec$Counts > 3)] <- "≥4"
rec$Counts <- factor(rec$Counts, levels=c("0", "1", "2", "3", "≥4"))
res <- data.frame()
for (type in types[3:5])
{
	ids <- which(rec$Type == type & rec$Counts != "0")
	res <- rbind(res, data.frame(Type=type, Total=length(ids), table(rec$Counts[ids])))
}
res <- res[which(res$Var1 != "0"),]
res$Rate <- res$Freq*100/res$Total
res$Type <- factor(res$Type, levels=types[3:5])
res$Var1 <- factor(as.character(res$Var1), levels=c("1", "2", "3", "≥4"))
res$Group <- "B"
res_ttb <- res
pttb <- ggplot(res_ttb, aes(x=Type, y=Rate, fill=Var1))+geom_bar(stat="identity", position=position_dodge(0.7), width=0.6)+
	labs(title="# OR genes\n(by nanopore)", x=NULL, y=NULL, fill="Num.")+
	scale_fill_manual(values=brewer.pal(9,"YlGnBu")[c(9,7,5,3)])+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_y_continuous(breaks=c(25, 50, 75), expand=c(0, 0))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), axis.ticks.x=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title.y=element_text(size=title_size, colour="black"), 
	axis.title.x=element_text(size=title_size, colour="black", hjust=1.3), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks.y=element_line(linewidth=0.35, color="black"))+tag_thm

ids <- match(colnames(olfr_sct), colnames(osn_rna))
rec <- data.frame(Cell=colnames(olfr_sct), X=osn_rna@reductions$umap@cell.embeddings[ids, 1], 
	Y=osn_rna@reductions$umap@cell.embeddings[ids, 2], Type=osn_rna$cell.subtype_fix[ids], Gene=0, Isoform=0)
for (i in 1:ncol(olfr_sct))
{
	terms <- olfr_sct_info[which(olfr_sct[, i] > 0),, drop=F]
	rec$Gene[i] <- length(table(terms$Gene))
	rec$Isoform[i] <- length(table(terms$ID))
}
rec$Counts <- rec$Isoform/rec$Gene
rec$Counts[which(is.nan(rec$Counts))] <- 0
rec$Counts[which(rec$Counts == 0)] <- "0"
rec$Counts[which(rec$Counts == 1)] <- "1"
rec$Counts[which(rec$Counts == 2)] <- "2"
rec$Counts[which(rec$Counts == 3)] <- "3"
rec$Counts[which(rec$Counts > 3)] <- "≥4"
rec$Counts <- factor(rec$Counts, levels=c("0", "1", "2", "3", "≥4"))
res <- data.frame()
for (type in types[3:5])
{
	ids <- which(rec$Type == type & rec$Counts != "0")
	res <- rbind(res, data.frame(Type=type, Total=length(ids), table(rec$Counts[ids])))
}
res <- res[which(res$Var1 != "0"),]
res$Rate <- res$Freq*100/res$Total
res$Type <- factor(res$Type, levels=types[3:5])
res$Var1 <- factor(as.character(res$Var1), levels=c("1", "2", "3", "≥4"))
res$Group <- "C"
res_ttc <- res
pttc <- ggplot(res_ttc, aes(x=Type, y=Rate, fill=Var1))+geom_bar(stat="identity", position=position_dodge(0.7), width=0.6)+
	labs(title="# OR isoform\n(per gene)", x=NULL, y=NULL, fill="Num.")+
	scale_fill_manual(values=brewer.pal(9,"YlGnBu")[c(9,7,5,3)])+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_y_continuous(breaks=c(25, 50, 75), expand=c(0, 0))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), axis.line=element_line(linewidth=0.35, color="black"), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(linewidth=0.35, color="black"))+tag_thm

res_ttm <- rbind(res_tta, res_ttb, res_ttc)
res_ttm$Group <- factor(res_ttm$Group, levels=c("A", "B", "C"), 
	labels=c("# OR genes\n(by illumina)", "# OR genes\n(by nanopore)", "# OR isoform\n(per gene)"))
pttm <- wrap_elements(ggplot(res_ttm, aes(x=Type, y=Rate, fill=Var1))+
	geom_bar(stat="identity", position=position_dodge(0.7), width=0.6)+
	facet_wrap(~Group, nrow=1)+labs(title=NULL, x=NULL, y=NULL, fill="Num.")+
	scale_fill_manual(values=brewer.pal(9,"YlGnBu")[c(9,7,5,3)])+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_y_continuous(breaks=seq(20, 80, 20), expand=c(0, 0))+
	theme(panel.border=element_rect(linewidth=0.35, fill='transparent',colour='black'), panel.background=element_blank(), 
	legend.key.size=unit(12, "pt"), legend.background=element_blank(), 
	strip.background=element_blank(), strip.text=element_text(size=title_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), axis.title.y=element_text(size=title_size, colour="black"), 
	axis.title.x=element_text(size=title_size, colour="black", hjust=1.3), axis.line=element_blank(), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(linewidth=0.35, color="black"))+tag_thm)+tag_thm

pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pac, pad, pak), nrow=1, widths=c(0.8, 1, 0.5, 0.8))+
	plot_annotation(tag_levels=list(c("A", "", "B", "", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pae, pxe), nrow=1, widths=c(1, 1))+
	plot_annotation(tag_levels=list(c("C", "E")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(paf, pttm), nrow=1, widths=c(1, 1))+
	plot_annotation(tag_levels=list(c("F", "G")), theme=tag_thm))+tag_thm, 
	pblank), ncol=1, heights=c(1, 1, 1, 1.6)), 
	width=210, height=297, dpi=300, units="mm", filename="oe_fl_fig_S04.png", limitsize=F)

ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pac, pad, pak), nrow=1, widths=c(0.8, 1, 0.5, 0.8))+
	plot_annotation(tag_levels=list(c("A", "", "B", "", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pae, pxe), nrow=1, widths=c(1, 1))+
	plot_annotation(tag_levels=list(c("C", "E")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(paf, pttm), nrow=1, widths=c(1, 1))+
	plot_annotation(tag_levels=list(c("F", "G")), theme=tag_thm))+tag_thm, 
	pblank), ncol=1, heights=c(1, 1, 1, 1.6)), 
	width=210, height=297, dpi=300, units="mm", filename="oe_fl_fig_S04.pdf", limitsize=F)

