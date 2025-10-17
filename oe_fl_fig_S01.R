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
library(extrafont)
options(stringsAsFactors=FALSE)

col_list <- c(c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8"), brewer.pal(12,"Set3")[-c(2, 9)])
text_size <- 13
title_size <- 15
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, colour="black"), plot.margin=margin(-3,-3,-3,-3), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", colour=NA),  plot.background=element_rect(fill="transparent", colour=NA), 
	legend.box.spacing=unit(0, "pt"))


gtf <- read.delim("gtf_info.tsv")
types <- c("Monocyte", "Macrophage", "BC", "RC", "Microvillar", "Neutrophils", "Basophil", 
	"Osteogenic", "Pericytes", "Ensheathing", "Sustentacular", "Bowman", 
	"HBC", "GBC", "INP", "Immature", "Mature")

cell_col <- col_list[1:length(types)]
names(cell_col) <- types
types_osn <- c("HBC", "GBC", "INP", "Immature", "Mature")
sce_rna <- readRDS("~/work/OEPac/OER_fix.rds")
sce_nano <- readRDS("oe_nano_fix.rds")
cells <- intersect(colnames(sce_rna), colnames(sce_nano))
sce_rna <- subset(sce_rna, cells=cells)
sce_rna$cell.subtype_fix <- as.character(sce_rna$cell.subtype_fix)
sce_rna$cell.subtype_fix[which(sce_rna$cell.subtype_fix == "LANSC")] <- "BC"
sce_rna$cell.subtype_fix <- factor(sce_rna$cell.subtype_fix, levels=types)
sce_nano <- subset(sce_nano, cells=cells)
sce_nano$cell.subtype_fix <- as.character(sce_nano$cell.subtype_fix)
sce_nano$cell.subtype_fix[which(sce_nano$cell.subtype_fix == "Brush")] <- "GBC"
sce_nano$cell.subtype_fix[which(sce_nano$cell.subtype_fix == "LANSC")] <- "BC"
sce_nano$cell.subtype_rna[which(sce_nano$cell.subtype_rna == "LANSC")] <- "BC"
sce_nano$cell.subtype_fix <- factor(sce_nano$cell.subtype_fix, levels=types)

cs_info <- data.frame(sce_rna@meta.data)
pa <- wrap_elements(ggplot(cs_info, aes(x="Term", y=nFeature_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Genes", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), panel.background=element_blank()))+tag_thm
pb <- wrap_elements(ggplot(cs_info, aes(x="Term", y=nCount_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Counts", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), panel.background=element_blank()))+tag_thm
pc <- wrap_elements(ggplot(cs_info, aes(x="Term", y=cell.mt))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Mito. (%)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), panel.background=element_blank()))+tag_thm

res_info <- read.delim("ass_type_info.tsv", h=T)
res_info$sample <- "ONT_OE"
res_info$count[1] <- res_info$count[1] + res_info$count[2]
res_info <- res_info[-2,]
res_info$rate <- res_info$count*100/sum(res_info$count)
res_info$type <- factor(res_info$type, levels=rev(res_info$type), labels=rev(c("FSM", "ISM3'", "ISM5'", "ISM\ninternal", "Mono\nexonic")))
pd <- wrap_elements(ggplot(res_info, aes(x=type, y=rate))+
	geom_bar(stat="identity", position=position_dodge(0.8), fill=col_list[4])+
	labs(title=NULL, x=NULL, y="Reads percentage", fill="Type")+
	#scale_fill_manual(values=rev(c("#db5f56", "#8ddd3e", "#fdb54e", "#7dcac1", "#5684da")), guide=guide_legend(reverse=T), drop=F)+
	scale_y_continuous(expand=c(0, 0))+
	theme(legend.position="none", axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

#reads_info <- read.csv("/mnt/md0/oe_full_length/output/OEfulllength/tagged_len.csv", h=F)
#reads_info_sub <- reads_info[sample(1:nrow(reads_info), 1000000),]
#terms <- c()
#cell_info <- data.frame(Cell=colnames(sce_nano), Len=0, Isoform=sce_nano$nFeature_RNA)
#for (i in 1:nrow(cell_info))
#{
#	if (i %% 1000 == 0) print(i)
#	ids <- which(reads_info_sub[, 3] == rownames(cell_info)[i])
#	terms <- c(terms, ids)
#	cell_info$Len[i] <- mean(reads_info[ids, 2])
#}
#terms <- reads_info_sub[terms, 2]
#write.csv(terms, "sce_nano_len_sub.csv", quote=F)
#write.csv(cell_info, "sce_nano_info.csv", quote=F)
terms <- read.csv("sce_nano_len_sub.csv", r=1)[, 1]
cell_info <- read.csv("sce_nano_info.csv", r=1, h=T)
pe <- wrap_elements(ggplot(, aes(x=terms, fill=type, color=type))+geom_density(fill=col_list[2], color=col_list[2], alpha=0.7, linewidth=1)+
	scale_y_continuous(expand=c(0, 0))+scale_x_continuous(breaks=seq(0, 3000, 500), limits=c(0, 2000), expand=c(0, 0))+
	labs(title=NULL, x="Length of tagged reads", y="Density\n")+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title=element_text(size=title_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), axis.line=element_line(colour="black")))+tag_thm
pf <- wrap_elements(ggplot(cell_info, aes(x="Term", y=Len))+geom_violin(fill=col_list[2], color=col_list[2])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Length (mean)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), panel.background=element_blank()))+tag_thm
pg <- wrap_elements(ggplot(cell_info, aes(x="Term", y=Isoform))+geom_violin(fill=col_list[2], color=col_list[2])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Isoforms", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), panel.background=element_blank()))+tag_thm

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

rs <- rowSums(osn_sct_raw[["trans"]] > 0)
res_ph <- data.frame(table(osn_sct_raw[["features"]]$Gene[intersect(which(rs > 2), which(osn_sct_raw[["features"]]$Type == "protein_coding"))]))
res_ph <- data.frame(table(res_ph$Freq))
res_ph$Rate <- res_ph$Freq * 100 / sum(res_ph$Freq)
ph <- wrap_elements(ggplot(res_ph, aes(x=Var1, y=Rate))+
	geom_bar(stat="identity", position=position_dodge(0.8), fill=col_list[4])+
	labs(title=NULL, x="# isoforms", y="Percentage of genes (%)")+
	scale_y_continuous(limits=c(0, 85), expand=c(0, 0))+
	geom_text(aes(x=Var1, y=Rate+2, label=Freq), size=4, vjust=0)+guides(fill="none")+
	theme(legend.position="none", axis.line=element_line(linetype=1, colour='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, E=pe, F=pf, G=pg, H=ph, design="ABCD\nEFGH", widths=c(2,2,2,3.5))+
	plot_annotation(tag_levels=list(c("A", "B", "C", "D", "E", "F", "G", "H")))+tag_thm, 
	width=13, height=8, dpi=200, filename="oe_fl_fig_S01.png", limitsize=F)

