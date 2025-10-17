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
tag_thm <- theme(plot.tag=element_text(size=title_size, color="black"), plot.margin=margin(-5,-5,-5,-5), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))

#gtf <- read.table("genes.gtf", sep="\t")
#gtf <- gtf[grep("^chr.*", gtf[, 1]),]
#colnames(gtf) <- c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
#gtf$transcript_id <- gsub(";.*", "", gsub(".*transcript_id ", "", gtf$attributes))
#gtf$gene_id <- gsub(";.*", "", gsub(".*gene_id ", "", gtf$attributes))
#gtf$gene_name <- gsub(";.*", "", gsub(".*gene_name ", "", gtf$attributes))
#gtf <- read.delim("gtf_info.tsv")
#
#osn_rna <- readRDS("osn_rna2_fix.rds")
##rec <- data.frame(bc=colnames(osn_rna), type="mOSN)
##write.table(rec, "/mnt/md0/oe_full_length/output/OSNfulllength/osn_bc_info.tsv", col.names=F, row.names=F, quote=F, sep="\t")
#data_raw <- H5File$new("/mnt/md0/oe_full_length/output/OSNfulllength/osn_trans_ass.h5", mode="r")
#sct <- list()
#sparse.mat <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
#	x=as.numeric(data_raw[["matrix/data"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
#rownames(sparse.mat) <- data_raw[["matrix/genes"]][]
#colnames(sparse.mat) <- data_raw[["matrix/barcodes"]][]
#sparse.utr3 <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
#	x=as.numeric(data_raw[["matrix/utr3"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
#rownames(sparse.utr3) <- data_raw[["matrix/genes"]][]
#colnames(sparse.utr3) <- data_raw[["matrix/barcodes"]][]
#sparse.utr5 <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
#	x=as.numeric(data_raw[["matrix/utr5"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
#rownames(sparse.utr5) <- data_raw[["matrix/genes"]][]
#colnames(sparse.utr5) <- data_raw[["matrix/barcodes"]][]
#sct[["matrix"]] <- as.sparse(sparse.mat)
#sct[["utr3"]] <- as.sparse(sparse.utr3)
#sct[["utr5"]] <- as.sparse(sparse.utr5)
#sct[["features"]] <- data.frame(Name=data_raw[["matrix/features/name"]][], ID=data_raw[["matrix/features/id"]][], 
#	Chr=data_raw[["matrix/features/chr"]][], Strand=data_raw[["matrix/features/strand"]][], Gene=data_raw[["matrix/features/gene"]][], 
#	Type=data_raw[["matrix/features/type"]][], UTR3=data_raw[["matrix/features/utr3"]][], UTR5=data_raw[["matrix/features/utr5"]][], 
#	Body=data_raw[["matrix/features/body"]][], Exon=data_raw[["matrix/features/exon"]][])
#sct[["features"]]$Count <- rowSums(sct[["matrix"]])
#sct[["features"]]$Symbol <- gtf$gene_name[match(sct[["features"]]$Gene, gtf$gene_id)]
#data_raw$close_all()
#osn_sct_raw <- list()
#osn_sct_raw[["trans"]] <- sct[["matrix"]]
#osn_sct_raw[["utr3"]] <- sct[["utr3"]]
#osn_sct_raw[["utr5"]] <- sct[["utr5"]]
#osn_sct_raw[["features"]] <- sct[["features"]]
#osn_sct_raw[["features"]]$ID <- sct[["features"]]$Name
#olfr_info <- read.csv("/mnt/md0/oe_full_length/output/OSNfulllength/olfr_info.csv", h=T)
#olfr_info <- olfr_info[match(osn_sct_raw[["features"]]$ID[grep("^Olfr", osn_sct_raw[["features"]]$Symbol)], olfr_info$ID),]
#olfr_info <- olfr_info[order(olfr_info$UTR),]
#olfr_info <- olfr_info[order(olfr_info$TSS),]
#olfr_info <- olfr_info[order(olfr_info$TypeE),]
#olfr_info <- olfr_info[order(olfr_info$TypeC),]
#olfr_info <- olfr_info[order(olfr_info$Gene),]
#olfr_info$CSS <- "Same"
#olfr_info$CSS[grep("A3SS", olfr_info$TypeC)] <- "Diff"
#olfr_info$CSS[grep("SE", olfr_info$TypeC)] <- "Diff"
#write.table(olfr_info, "osn_iso_info_all.tsv", sep="\t", quote=F, row.names=F)
#olfr_sct <- osn_sct_raw[["trans"]][match(olfr_info$ID, osn_sct_raw[["features"]]$ID),]
#olfr_sct <- olfr_sct[, which(colSums(olfr_sct) > 0)]
#olfr_sct_info <- olfr_info[match(rownames(olfr_sct), olfr_info$ID),]
#olfr_sct_info$Count <- rowSums(olfr_sct)
#
#count <- 0
#rec <- c()
#isoforms <- c()
#cells_nano <- data.frame()
#for (i in 1:ncol(olfr_sct))
#{
#	terms <- olfr_sct_info[which(olfr_sct[, i] > 0),, drop=F]
#	ids <- which(terms$TypeC == "-")
#	if (length(ids) == 0)
#	{
#		if (length(unique(terms$Gene)) == 1) cells_nano <- rbind(cells_nano, data.frame(Cell=colnames(olfr_sct)[i], 
#			Gene=terms$Gene[1], Symbol=terms$Symbol[1], Num=nrow(terms), Same=nrow(terms), Other=0, Main=0, Sub=0, 
#			Isoform=paste0(terms$ID, collapse=",")))
#		next
#	}
#	if (length(unique(terms$Gene[ids])) > 1)
#	{
#		count <- count + 1
#		next
#	}
#	nf <- terms$ID[ids[which.max(olfr_sct[terms$ID[ids], i])]]
#	ns <- terms$ID[which(terms$Gene == terms$Gene[ids[1]] & terms$ID != nf)]
#	nf <- olfr_sct[nf, i]
#	nsc <- 0
#	if (length(ns) > 0) nsc <- sum(olfr_sct[ns, i])
#	cells_nano <- rbind(cells_nano, data.frame(Cell=colnames(olfr_sct)[i], Gene=terms$Gene[ids[1]], Symbol=terms$Symbol[ids[1]], 
#		Num=nrow(terms), Same=length(ns), Other=nrow(terms)-length(ns)-1, Main=nf, Sub=nsc, Isoform=paste0(terms$ID, collapse=",")))
#	ids <- which(terms$TSS == "-")
#	if (length(ids) > 0) isoforms <- c(isoforms, terms$ID[-1])
#	ids <- which(terms$Gene == terms$Gene[ids[1]] & terms$CSS == "Same" & terms$Count > 0)
#	if (length(ids) > 1) rec <- c(rec, terms$ID[ids])
#}
#isoforms <- olfr_sct_info[match(unique(isoforms), olfr_sct_info$ID),]
#rec <- olfr_sct_info[match(unique(rec), olfr_sct_info$ID),]
#ids <- match(rec$ID, osn_sct_raw[["features"]]$ID)
#rec$Chr <- osn_sct_raw[["features"]]$Chr[ids]
#rec$Strand <- osn_sct_raw[["features"]]$Strand[ids]
#rec$UL <- osn_sct_raw[["features"]]$UTR3[ids]+osn_sct_raw[["features"]]$UTR5[ids]
#for (i in 1:nrow(rec))
#{
#	if (rec$UL[i] != 0) next
#	ids <- which(olfr_info$Gene == rec$Gene[i] & olfr_info$CDS != "")
#	if (length(ids) == 0) print("???")
#	strand <- substr(rec$Gene[i], nchar(rec$Gene[i]), nchar(rec$Gene[i]))
#	ss <- as.numeric(unlist(strsplit(unlist(strsplit(olfr_info$CDS[ids[1]], split="-")), split=";")))
#	exons <- unlist(strsplit(rec$Exon[i], split=";"))
#	cds <- c()
#	if (strand == "+")
#	{
#		ss <- ss[1]
#		for (e in exons)
#		{
#			e <- as.numeric(unlist(strsplit(e, split="-")))
#			if (e[2] < ss) next
#			cds <- c(cds, paste0(max(e[1], ss), "-", e[2]))
#		}
#	} else {
#		ss <- ss[length(ss)]
#		for (e in exons)
#		{
#			e <- as.numeric(unlist(strsplit(e, split="-")))
#			if (e[1] > ss) break
#			cds <- c(cds, paste0(e[1], "-", min(e[2], ss)))
#		}
#	}
#	rec$CDS[i] <- paste0(cds, collapse=";")
#}
#iso_trans_test <- rec
#write.table(isoforms, "osn_iso_info.tsv", sep="\t", quote=F, row.names=F)
#write.table(cells_nano, "osn_olfr_info.tsv", sep="\t", quote=F, row.names=F)
#write.table(iso_trans_test, "osn_iso_trans_test.tsv", sep="\t", quote=F, row.names=F)
################################################################################
#sce <- CreateSeuratObject(Read10X("/mnt/md0/oe_full_length/output/OSNfulllength/OSNfulllength.gene_raw_feature_bc_matrix"), project="nano")
#length(which(colSums(sce[["RNA"]]$counts[grep("^Olfr", rownames(sce[["RNA"]]$counts)), colnames(osn_sct_raw[["trans"]])]) == 0))
#172
#sce <- CreateSeuratObject(Read10X("/mnt/md0/oe_full_length/output/OSNfulllength/OSNfulllength.transcript_processed_feature_bc_matrix"), project="nano")
#length(which(colSums(sce[["RNA"]]$counts[grep("^Olfr", gtf$gene_name[match(rownames(sce), gtf$transcript_id)]), 
#	intersect(colnames(osn_sct_raw[["trans"]]), colnames(sce))]) == 0))
#1847
#length(which(colSums(osn_rna[["RNA"]]$counts[grep("^Olfr", rownames(osn_rna[["RNA"]]$counts)),]) == 0))
#length(which(colSums(osn_sct_raw[["trans"]][grep("^Olfr", osn_sct_raw[["features"]]$Symbol),]) == 0))
#58/6206/9848 (RNA/Nano/Total)
#ids_rna <- grep("^Olfr", rownames(osn_rna[["RNA"]]$counts))
#ra_mat <- osn_rna[["RNA"]]$counts[ids_rna, cells_nano$Cell]
#ids <- match(cells_nano$Symbol, rownames(ra_mat))
#rec <- c()
#for (i in 1:nrow(cells_nano))
#{
#	if (is.na(ids[i])) next
#	if (ra_mat[ids[i], i] > 0) rec <- c(rec, i)
#}
#3224/3394 (same/total)
#6207/247/3394/9848 (none/multi/single/total)
################################################################################
# None	8205	60.5%
# Multi	284	2.1%
# Single	5068	37.4%
# Total	13557	-
#cells <- setdiff(names(which(colSums(sce[["RNA"]]$counts[grep("^Olfr", gtf$gene_name[match(rownames(sce), gtf$transcript_id)]), 
#	intersect(colnames(osn_sct_raw[["trans"]]), colnames(sce))])> 0)), cells_nano$Cell)
#write.table(cells, "osn_null_test.tsv", row.names=F, col.names=F, quote=F)
#null_info <- read.delim("osn_null_res.tsv", h=T)
#null_info_unique <- null_info[which(null_info$type == "unique"),]
#or_test_ra <- data.frame(Group=rep(c("Detected", "None"), each=3), Type=rep(c("None", "Multi", "Single"), 2), Number=c(5587, 284, 5068, 2618, 0, 0))
#or_test_ra$Type <- factor(or_test_ra$Type, levels=c("None", "Multi", "Single"))
#or_test_rb <- data.frame(Type=c("Ambiguous", "Inconsistent", "Intro", "Mono_no_TSS"), Number=c(119425, 37901, 14036, 13726))
#or_test_rb$Rate <- or_test_rb$Number*100/sum(or_test_rb$Number)
#or_test_rb$Type <- factor(or_test_rb$Type, levels=or_test_rb$Type)
#p_ora <- ggplot(or_test_ra, aes(x=Type, y=Number, fill=Group))+geom_bar(stat="identity", width=0.6)+
#	labs(title=NULL, x="OR types of isoforms", y="Number of cells", fill="OR Gene")+
#	scale_fill_manual(values=col_list[c(2, 6)], drop=F)+scale_y_continuous(limits=c(0, 8600), expand=c(0, 0))+
#	geom_text(data=data.frame(Type=c("None", "Multi", "Single"), Number=c(8205, 284, 5068), Group="Detected"), 
#	aes(x=Type, y=Number+30, label=Number), size=5, vjust=0)+
#	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
#	legend.title=element_text(size=title_size, colour="black"), 
#	legend.text=element_text(size=text_size, colour="black"), 
#	legend.key=element_blank(), legend.background=element_blank(), 
#	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
#	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
#	axis.line=element_line(colour="black"), legend.position=c(0.8, 0.9))
#p_orb <- ggplot(or_test_rb, aes(y=Type, x=Rate))+geom_bar(stat="identity", width=0.6, fill=col_list[3])+
#	labs(title=NULL, y=NULL, x="Percentage of OR reads (%)")+scale_x_continuous(expand=c(0, 0))+
#	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
#	legend.title=element_text(size=title_size, colour="black"), 
#	legend.text=element_text(size=text_size, colour="black"), 
#	legend.key=element_blank(), legend.background=element_blank(), 
#	axis.ticks.y=element_blank(), axis.ticks.x=element_line(colour="black"), 
#	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
#	axis.line=element_line(colour="black"))
#ggsave(plot=p_ora+p_orb, width=8, height=5, dpi=200, filename="osn_null_res.png", limitsize=F)
################################################################################
isoforms <- rbind(data.frame(read.delim("osn_iso_info.tsv", h=T), Sample="OSN"), data.frame(read.delim("oe_iso_info.tsv", h=T), Sample="OE"))
cells_nano <- rbind(data.frame(read.delim("osn_olfr_info.tsv", h=T), Sample="OSN"), data.frame(read.delim("oe_olfr_info.tsv", h=T), Sample="OE"))
iso_trans_test <- rbind(data.frame(read.delim("osn_iso_trans_test.tsv", h=T), Sample="OSN"), data.frame(read.delim("oe_iso_trans_test.tsv", h=T), Sample="OE"))

rec_aa <- data.frame(table(cells_nano$Num))
paa <- ggplot(rec_aa, aes(x=Var1, y=Freq))+
	geom_bar(stat="identity", width=0.6, position=position_dodge(0.8), fill=col_list[1])+
	labs(title=NULL, x="OR isoform of OSNs", y="Number of cells")+
	scale_y_continuous(limits=c(0, max(rec_aa$Freq)+300), expand=c(0, 0))+guides(fill="none")+
	geom_text(aes(x=Var1, y=Freq+20, label=Freq), size=4, vjust=0)+
	theme(axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank())+tag_thm
rec_ab <- data.frame(Var1=c("AS", "ATSS"), Freq=as.numeric(table(isoforms$TSS)))
pab <- ggplot(rec_ab, aes(x=4, y=Freq, fill=Var1))+geom_col()+
	labs(title=NULL, x=NULL, y=NULL, fill=NULL)+
	geom_text(aes(label=Freq), size=4, position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+scale_x_continuous(limits=c(2.5, 4.5), expand=c(0, 0)) +
	scale_fill_manual(values=col_list[c(6,2)])+
	theme(plot.title=element_text(size=title_size, hjust=0), panel.background=element_blank(), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	legend.position="bottom", legend.direction="vertical", 
	plot.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), axis.title=element_blank(), legend.title=element_text(colour="black", size=title_size), 
	legend.text=element_text(colour="black", size=text_size), plot.margin=margin(-10,-10,-10,-10), panel.spacing=unit(0, "pt"))
rec_ac <- data.frame(Var1=c("UTR", "CDS"), Freq=c(length(which(isoforms$TypeC == "-")), length(which(isoforms$TypeC != "-"))))
pac <- ggplot(rec_ac, aes(x=4, y=Freq, fill=Var1))+geom_col()+
	labs(title="Type of isoforms", x=NULL, y=NULL, fill=NULL)+
	geom_text(aes(label=Freq), size=4, position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+scale_x_continuous(limits=c(2.5, 4.5), expand=c(0, 0)) +
	scale_fill_manual(values=col_list[c(4,7)])+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	legend.position="bottom", legend.direction="vertical", 
	plot.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), axis.title=element_blank(), legend.title=element_text(colour="black", size=title_size), 
	legend.text=element_text(colour="black", size=text_size), plot.margin=margin(-10,-10,-10,-10), panel.spacing=unit(0, "pt"))
rec_ad <- data.frame(table(cells_nano$Other[which(cells_nano$Num > 1)] > 0))
rec_ad$Var1 <- c("Same", "Diff.")
pad <- ggplot(rec_ad, aes(x=4, y=Freq, fill=Var1))+geom_col()+
	labs(title=NULL, x=NULL, y=NULL, fill=NULL)+
	geom_text(aes(label=Freq), size=4, position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+scale_x_continuous(limits=c(2.5, 4.5), expand=c(0, 0)) +
	scale_fill_manual(values=col_list[c(3,12)])+
	theme(plot.title=element_text(size=title_size, hjust=0), panel.background=element_blank(), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	legend.position="bottom", legend.direction="vertical", 
	plot.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), axis.title=element_blank(), legend.title=element_text(colour="black", size=title_size), 
	legend.text=element_text(colour="black", size=text_size), plot.margin=margin(-10,-10,-10,-10), panel.spacing=unit(0, "pt"))
pae <- wrap_plots(list(pab, pac, pad), nrow=1)+theme(plot.margin=margin(-10,-10,-10,-10), panel.spacing=unit(0, "pt"))
pa <- wrap_elements(paa+inset_element(pae, 0.2, 0.1, 1, 1)+tag_thm)+tag_thm

ids <- which(cells_nano$Sub > 0)
rec_ba <- data.frame(Cell=cells_nano$Cell[ids], Main=cells_nano$Main[ids], Sub=cells_nano$Sub[ids])
rec_ba$Count <- rec_ba$Main + rec_ba$Sub
rec_ba$Main <- rec_ba$Main*100/rec_ba$Count
rec_ba$Sub <- rec_ba$Sub*100/rec_ba$Count
rec_ba <- rec_ba[order(rec_ba$Main),]
rec_ba$Cell <- factor(rec_ba$Cell, levels=rec_ba$Cell)
rec_bb <- data.frame(Cell=rep(rec_ba$Cell, 2), Group=rep(c("Essential ", "Trivial"), each=nrow(rec_ba)), Rate=c(rec_ba$Main, rec_ba$Sub))
paca <- ggplot(rec_ba, aes(x=Cell, y=log2(Count)))+geom_bar(stat="identity", fill="gray50", color="gray50")+
	labs(title=NULL, x=NULL, y="# (log2)")+scale_y_continuous(expand=c(0, 0))+guides(fill="none")+
	theme(axis.line.x=element_blank(), axis.line.y=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.x=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.margin=margin(-3,0,0,0))
pacb <- ggplot(rec_bb, aes(x=Cell, y=Rate, fill=Group, color=Group))+geom_bar(stat="identity")+
	labs(title=NULL, x="OR isoform of OSNs", y="Percentage (%)", fill=NULL, color=NULL)+
	scale_y_continuous(breaks=seq(0, 90, 20), expand=c(0, 0))+
	scale_fill_manual(values=col_list[c(4,6)])+scale_color_manual(values=col_list[c(4,6)])+
	theme(axis.line.x=element_blank(), axis.line.y=element_line(linetype=1, colour='black'), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.x=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.margin=margin(3,0,0,0))
pb <- wrap_elements(ggarrange(paca, pacb, align="v", common.legend=T, legend="top", ncol=1, heights=c(1,3)))+tag_thm

info_col <- col_list[c(2,4,5,6,8)]
names(info_col) <- c("A3SS", "A5SS", "RI", "SE", "UTR")
res_mat <- data.frame(matrix(FALSE, nrow=nrow(isoforms), ncol=4, dimnames=list(isoforms$ID, c("A3SS", "A5SS", "RI", "SE"))))
res_mat$A3SS[grep("A3SS", isoforms$TypeC)] <- TRUE
res_mat$A5SS[grep("A5SS", isoforms$TypeC)] <- TRUE
res_mat$RI[grep("RI", isoforms$TypeC)] <- TRUE
res_mat$SE[grep("SE", isoforms$TypeC)] <- TRUE
res_mat <- res_mat[which(isoforms$TypeC != "-"),]
pc <- upset(res_mat, colnames(res_mat), name="Overlap", width_ratio=0.25, wrap=T, stripes="white", 
	base_annotations=list("Intersection size"=intersection_size(fill=col_list[4], color=col_list[4], 
	text=list(size=4))+scale_y_continuous(expand=c(0, 0))), 
	set_sizes=(upset_set_size()+scale_y_reverse(breaks=c(0, 20, 40), expand=c(0, 0))+theme(panel.background=element_blank(), 
	panel.grid=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), 
	axis.title=element_text(size=title_size, color="black"), axis.text.x=element_text(size=text_size, color="black"), 
	axis.text.y=element_blank(), axis.ticks.x=element_line(color="black"), axis.line.x=element_line(color="black"), 
	plot.margin=margin(-5,-5,-5,-5), panel.spacing=unit(0, "pt"))), 
	queries=list(upset_query(set="A3SS", fill=info_col["A3SS"], color=info_col["A3SS"]), 
	upset_query(set="A5SS", fill=info_col["A5SS"], color=info_col["A5SS"]), 
	upset_query(set="RI", fill=info_col["RI"], color=info_col["RI"]), 
	upset_query(set="SE", fill=info_col["SE"], color=info_col["SE"])), 
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

data_raw <- read.delim("or_test_exam.tsv")
plsb <- lapply(c("A3SS", "A5SS", "RI", "SE"), function(type) {
	ids <- which(data_raw$Type == type)
	seq_id <- data_raw$Chr[ids[1]]
	strand <- data_raw$Strand[ids[1]]
	gene_id <- data_raw$Gene[ids[1]]
	info_exons <- data.frame()
	trans <- data_raw$ID[ids]
	for (t in trans)
	{
		end <- 0
		exons <- unlist(strsplit(paste0(data_raw$Exon[which(data_raw$ID == gsub("-", "_", t))], ";"), split=";"))
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
	regions <- c(min(info_exons$start), max(info_exons$end))
	#for (t in 2:nrow(info_exons)) if (info_exons$start[t] == info_exons$end[t-1]) info_exons$start[t] <- info_exons$end[t-1] + 1
	info_exons <- info_exons[order(info_exons$type, decreasing=T),]
	info_exons <- info_exons[order(info_exons$end),]
	info_exons <- info_exons[order(info_exons$start),]
	info_exons <- info_exons[order(info_exons$transcript_name),]
	info_exons <- info_exons[which(info_exons$end - info_exons$start > 0),]
	info_ranges <- paste0(info_exons$seqnames[1], strand, ": ", min(info_exons$start), " - ", max(info_exons$start))
	info_exons_raw <- info_exons
	info_exons <- shorten_gaps(exons=info_exons, introns=to_intron(info_exons, "transcript_name"), group_var="transcript_name", 
		target_gap_width=as.integer(median(info_exons$end - info_exons$start)))
	coord_trans <- data.frame(r=c(info_exons_raw$start, info_exons_raw$end), 
		i=c(info_exons$start[which(info_exons$type != "intron")], info_exons$end[which(info_exons$type != "intron")]))
	coord_trans <- coord_trans[!duplicated(coord_trans),]
	coord_trans <- coord_trans[order(coord_trans$r),]
	info_introns <- info_exons[which(info_exons$type == "intron"),]
	info_utrs <- info_exons[which(info_exons$type == "UTR"),]
	info_exons <- info_exons[which(info_exons$type == "exon"),]
	info_introns$transcript_name <- factor(info_introns$transcript_name, levels=trans)
	info_exons$transcript_name <- factor(info_exons$transcript_name, levels=trans)
	info_utrs$transcript_name <- factor(info_utrs$transcript_name, levels=trans)
	info_introns$rs <- coord_trans$r[match(info_introns$start, coord_trans$i)]
	info_introns$re <- coord_trans$r[match(info_introns$end, coord_trans$i)]
	info_utrs$rs <- coord_trans$r[match(info_utrs$start, coord_trans$i)]
	info_utrs$re <- coord_trans$r[match(info_utrs$end, coord_trans$i)]
	info_exons$rs <- coord_trans$r[match(info_exons$start, coord_trans$i)]
	info_exons$re <- coord_trans$r[match(info_exons$end, coord_trans$i)]
	return(ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
		geom_intron(data=info_introns, aes(strand=strand), color=info_col[type], linewidth=0.5, arrow.min.intron.length=100)+
		geom_range(fill=info_col[type], color=info_col[type], linewidth=0, height=0.3)+
		geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name), 
		fill=info_col[type], color=info_col[type], linewidth=0.5, height=0.1)+
		labs(title=NULL, x=NULL, y=type)+scale_y_discrete(expand=expansion(mult=c(0.5,1)))+
		#scale_x_continuous(breaks=(max(c(info_cds$end, info_utrs$end))+min(c(info_cds$start, info_utrs$start)))/2, labels=info_ranges)+
		geom_text(data=info_introns, aes(y=transcript_name, label=transcript_name), 
		color=info_col[type], x=0, hjust=0, vjust=0, nudge_y=0.3, size=4)+
		theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
		axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(size=title_size, colour="black"), 
		legend.position="none", plot.margin=margin(2,-3,2,-3), panel.spacing=unit(0, "pt")))
})
pd <- wrap_elements(wrap_plots(plsb, ncol=1)+tag_thm)+tag_thm

iso_trans_sub <- iso_trans_test[which(iso_trans_test$TypeC != "-"),]
iso_trans_sub$TypeC <- factor(iso_trans_sub$TypeC, levels=c("A5SS", "RI", "A5SS+RI"))
iso_trans_sub$ID <- paste(iso_trans_sub$ID, iso_trans_sub$Sample)
data_raw <- rbind(data.frame(read.delim("/mnt/md0/oe_full_length/output/OSNfulllength/iso_trans_test_aa_cal.tsv"), Sample="OSN"), 
	data.frame(read.delim("/mnt/md0/oe_full_length/output/OEfulllength/iso_trans_test_aa_cal.tsv"), Sample="OE"))
data_raw$Term <- paste(data_raw$Term, data_raw$Sample)
data_raw <- data_raw[match(iso_trans_sub$ID, gsub(".*[+-]_", "", data_raw$Term)),]
ranks <- data.frame(Term=data_raw$Term, ID=iso_trans_sub$ID, Type=iso_trans_sub$TypeC, Rank=0)
for (i in 2:8) ranks$Rank[which(data_raw[, i] > 40)] <- i
ranks <- ranks[order(ranks$Rank),]
data_df <- data.frame()
for (i in 2:8) data_df <- rbind(data_df, data.frame(Term=iso_trans_sub$ID, Type=iso_trans_sub$TypeC, 
	Pos=colnames(data_raw)[i], Rate=data_raw[, i]))
data_df$Type <- factor(data_df$Type, levels=c("A5SS", "RI", "A5SS+RI"))
data_df$Term <- factor(data_df$Term, levels=ranks$ID)
data_df$Rate[which(data_df$Rate > 40)] <- 100
data_df$Rate[which(data_df$Rate <= 40)] <- 0
data_df$Rate <- factor(data_df$Rate, levels=c(0, 100), labels=c("Lost", "Exist"))
ptd <- ggplot(data_df, aes(x=Pos, y=Term, fill=Rate))+geom_tile()+#geom_text(aes(label=Rate), colour="black")+
	facet_grid(Type~., scales="free_y", space="free_y")+
	scale_fill_manual(values=c("gray90", col_list[2]))+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	labs(title=NULL, x=NULL, y=NULL, fill="Domain")+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), 
	panel.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), strip.text=element_blank(), strip.background=element_blank(), 
	panel.spacing=unit(3, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), legend.text=element_text(size=text_size, colour="black"), 
	plot.margin=margin(0,0,0,0))
rec_te <- data.frame(Type=c("Domain1", "Domain2", "Domain3", "Domain4", "Domain5", "Domain6", "Domain7"), Info="Line")
pte <- ggplot(rec_te, aes(y="", x=Type, fill=Type))+geom_tile()+theme_minimal()+
	#geom_text(aes(label=Type), size=title_size*0.3, colour="white")+
	labs(title=NULL, y=NULL, x=NULL, fill=NULL)+
	scale_fill_manual(values=c("#ffa24b", "#39ff7f", "#3ba8b4", "#ff73c5", "#ff2443", "#3e3eff", "#ffd23b"))+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), panel.spacing=unit(3, "pt"), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	strip.text.y=element_text(size=text_size, colour="black"), legend.position="none", 
	strip.background=element_rect(colour="white", fill="white"), axis.text=element_blank(), 
	axis.ticks=element_blank(), plot.margin=margin(0,0,0,0))
ptf <- ggplot(iso_trans_sub, aes(x="", y=ID, fill=TypeC))+geom_tile()+theme_minimal()+
	labs(title=NULL, y=NULL, x=NULL, fill="Type")+
	facet_grid(TypeC~., scales="free_y", space="free_y")+
	scale_fill_manual(values=c("#6994b3", "#d1934b", "#9D947F"))+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), panel.spacing=unit(3, "pt"), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks=element_blank(), 
	strip.text=element_blank(), strip.background=element_blank(), axis.text=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), legend.text=element_text(size=text_size, colour="black"), 
	plot.margin=margin(0,0,0,0))
ptc <- wrap_elements(ggdraw()+draw_image("fig04_f.png", scale=1))+tag_thm
pf <- wrap_elements(wrap_plots(A=ptc, B=pte, C=ptd, D=ptf, design="AA\nB#\nCD", 
	heights=c(8,1,20), widths=c(20,1))+plot_layout(guides="collect"))+tag_thm

pe <- wrap_elements(ggdraw()+draw_image("fig04_e.png", scale=1.05))+theme(plot.tag=element_text(size=title_size, color="black"), 
	plot.margin=margin(-20,-5,0,-20), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA))

pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pb, pe), nrow=1, widths=c(1.4, 0.9, 1))+
	plot_annotation(tag_levels=list(c("A", "B", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pc, pd, pf), nrow=1, widths=c(1.2, 1.2, 1))+
	plot_annotation(tag_levels=list(c("C", "", "E")), theme=tag_thm))+tag_thm), 
	ncol=1, heights=c(1, 1)), width=13, height=8, dpi=200, filename="oe_fl_fig_04.png", limitsize=F)


