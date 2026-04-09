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
library(ggh4x)
options(stringsAsFactors=FALSE)

col_list <- c(c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8"), brewer.pal(12,"Set3")[-c(2, 9)])
text_size <- 8
title_size <- 9
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, face="bold", colour="black"), plot.margin=margin(-3,0,3,0), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", colour=NA),  plot.background=element_rect(fill="transparent", colour=NA), 
	legend.box.spacing=unit(0, "pt"))

gtf <- read.table("genes.gtf", sep="\t")
gtf <- gtf[grep("^chr.*", gtf[, 1]),]
colnames(gtf) <- c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gtf$transcript_id <- gsub(";.*", "", gsub(".*transcript_id ", "", gtf$attributes))
gtf$gene_id <- gsub(";.*", "", gsub(".*gene_id ", "", gtf$attributes))
gtf$gene_name <- gsub(";.*", "", gsub(".*gene_name ", "", gtf$attributes))
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

data_raw <- H5File$new("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_info.h5", mode="r")
scg <- list()
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
sparse.body <- sparseMatrix(i=data_raw[["matrix/indices"]][]+1, p=data_raw[["matrix/indptr"]][], 
	x=as.numeric(data_raw[["matrix/body"]][]), dims=data_raw[["matrix/shape"]][], repr = "T")
rownames(sparse.body) <- data_raw[["matrix/genes"]][]
colnames(sparse.body) <- data_raw[["matrix/barcodes"]][]
scg[["matrix"]] <- as.sparse(sparse.mat)
scg[["utr3"]] <- as.sparse(sparse.utr3)
scg[["utr5"]] <- as.sparse(sparse.utr5)
scg[["body"]] <- as.sparse(sparse.body)
scg[["features"]] <- data.frame(Name=data_raw[["matrix/features/name"]][], ID=data_raw[["matrix/features/id"]][], 
	Chr=data_raw[["matrix/features/chr"]][], Strand=data_raw[["matrix/features/strand"]][], Gene=data_raw[["matrix/features/gene"]][], 
	Trans=data_raw[["matrix/features/trans"]][])
scg[["features"]]$Count <- rowSums(scg[["matrix"]])
scg[["features"]]$Symbol <- gtf$gene_name[match(scg[["features"]]$Gene, gtf$gene_id)]
data_raw$close_all()

sce_trans <- CreateSeuratObject(osn_sct[["genes"]], project="OEnano")
sce_trans$cell.subtype_fix <- as.character(osn_rna$cell.subtype_fix[match(colnames(sce_trans), colnames(osn_rna))])
sce_trans <- SCTransform(sce_trans, method="glmGamPoi", variable.features.n=nrow(sce_trans))
mat_utr5 <- round(scg[["utr5"]]*scg[["body"]]/(scg[["body"]]+0.1))
mat_utr5 <- scg[["utr5"]][, match(colnames(osn_rna), colnames(scg[["utr5"]]))]
mat_utr5 <- mat_utr5[which(rowSums(mat_utr5 > 10) > 10),]
mat_utr5 <- round(scg[["utr5"]]*scg[["body"]]/(scg[["body"]]+0.1))
mat_utr5 <- scg[["utr5"]][, match(colnames(osn_rna), colnames(scg[["utr5"]]))]
mat_utr5 <- mat_utr5[which(rowSums(mat_utr5 > 10) > 10),]
rec <- data.frame(Term=rownames(mat_utr5), Gene=scg[["features"]]$Gene[match(rownames(mat_utr5), scg[["features"]]$Name)], 
	Trans=scg[["features"]]$Trans[match(rownames(mat_utr5), scg[["features"]]$Name)], 
	Symbol=scg[["features"]]$Symbol[match(rownames(mat_utr5), scg[["features"]]$Name)], Pval.e=1, FC.e=0, Pval=1, FC=0, Len=0, Diff=0)
rec <- rec[which(!is.na(match(rec$Gene, rownames(sce_trans)))),]
term_ids <- match(rec$Term, rownames(mat_utr5))
osn_rna$random <- as.numeric(sample(osn_rna$Pseudotime, ncol(osn_rna)))
ids_a <- which(osn_rna$cell.subtype_fix != "iOSN" & osn_rna$cell.subtype_fix != "mOSN")
ids_b <- which(osn_rna$cell.subtype_fix == "iOSN" | osn_rna$cell.subtype_fix == "mOSN")
cl <- makeCluster(16, type="FORK")
results <- parLapply(cl, 1:nrow(rec), function(i) 
{
	rec_sub <- rec[i,, drop=F]
	#tta <- sce_trans[["SCT"]]$scale.data[rec$Gene, ids_a]
	#ttb <- sce_trans[["SCT"]]$scale.data[rec$Gene, ids_b]
	#ttr = t.test(tta, ttb)
	#rec_sub$FC.e[1] <- log2(ttr$estimate[2]+0.00001)-log2(ttr$estimate[1]+0.00001)
	#rec_sub$Pval.e[1] <- ttr$p.value
	tta <- mat_utr5[term_ids[i], ids_a]
	ttb <- mat_utr5[term_ids[i], ids_b]
	ia <- which(tta > 10)
	ib <- which(ttb > 10)
	if (length(ia) < 10 | length(ib) < 10) return(rec_sub)
	tta <- tta[ia]
	ttb <- ttb[ib]
	if (length(unique(tta)) == 1 & length(unique(ttb)) == 1 & min(unique(tta)) == min(unique(ttb))) return(rec_sub)
	ttr <- wilcox.test(tta, ttb)
	#ttr <- t.test(tta, ttb)
	tt <- rbind(data.frame(Group="Early", Len=tta), data.frame(Group="Later", Len=ttb))
	ptt <- ggplot(tt, aes(x=Len, color=Group, fill=Group))+geom_density(alpha=0.4)+
		labs(title=paste0(scg[["features"]]$Trans[which(scg[["features"]]$Name == rec_sub$Term)], " (", rec_sub$Symbol, ")"), 
		x=paste0("Length (", length(tta), "/", length(ttb), ")"), y="Density", color="Group")+
		scale_fill_manual(values=c("#FC8D62", "#8DA0CB"))+
		scale_color_manual(values=c("#FC8D62", "#8DA0CB"))+
		scale_x_continuous(limits=c(min(tt$Len)-2, max(tt$Len)+2))+
		theme(plot.title=element_text(size=14, hjust=0.5, colour="black"), panel.background=element_blank(), 
		axis.line=element_line(colour="black"))
	#ggsave(plot=ptt, width=5, height=4, dpi=200, paste0("utr5_", rec_sub$Symbol, "_dist.png"), limitsize=F)
	rec_sub$FC[1] <- log2(mean(ttb)+0.00001)-log2(mean(tta)+0.00001)
	rec_sub$Diff[1] <- mean(ttb) - mean(tta)
	rec_sub$Pval[1] <- ttr$p.value
	rec_sub$Len[1] <- round(max(ttr$estimate))
	return(rec_sub)
})
stopCluster(cl)
rec <- do.call("rbind", results)

sce_trans$Group <- "B"
sce_trans$Group[ids_a] <- "A"
group_df <- FindMarkers(sce_trans, ident.1="A", ident.2="B", group.by="Group")
rec$Pval.e <- group_df$p_val[match(rec$Gene, rownames(group_df))]
rec$Padj.e <- group_df$p_val_adj[match(rec$Gene, rownames(group_df))]
rec$FC.e <- group_df$avg_log2FC[match(rec$Gene, rownames(group_df))]
rec$Padj <- p.adjust(rec$Pval, method="BH")
rec <- rec[order(rec$FC, decreasing=T),]
rec <- rec[order(rec$Padj),]
deg_info <- rec
deg_info$Group <- "N.S."
deg_info$Group[which(deg_info$Padj < 0.05 & deg_info$Diff > 30)] <- "Increased"
deg_info$Group[which(deg_info$Padj < 0.05 & deg_info$Diff < -30)] <- ""
#deg_info$X <- log2(abs(deg_info$Diff))*deg_info$Diff/abs(deg_info$Diff)
deg_info$X <- deg_info$Diff
deg_info$Y <- -log10(deg_info$Padj)
deg_info$Y[which(deg_info$Y > 10)] <- 10
deg_info$Rank <- abs(deg_info$X)
deg_info <- deg_info[order(deg_info$Rank, decreasing=T),]
deg_info$Anno <- ""
deg_info$Anno[which(deg_info$Symbol == "Scg2")] <- "Scg2"
deg_info$Group <- factor(deg_info$Group, levels=c("Increased", "", "N.S."))
deg_info <- deg_info[which(abs(deg_info$Diff) > 10),]
pba <- ggplot(deg_info, aes(x=X, y=Y, color=Group))+geom_point(size=2)+
	labs(title=NULL, x="Difference of 5'UTR length (bp)", y="P.adj (-log10)", color=NULL)+
	scale_color_manual(values=c("#FC8D62", "#8DA0CB", "gray80"), drop=F)+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	scale_x_continuous(breaks=c(-100, -30, 0, 30, 100), labels=c("-100", "-30     ", "0", "    30", "100"))+
	scale_y_continuous(limits=c(0, 11), breaks=c(1.3, seq(0, 10, 2.5)), labels=c("1.3", "0.0", "2.5", "5.0", "7.5", "≥\n10.0"), expand=c(0, 0))+
	geom_text_repel(data=deg_info[which(deg_info$Anno != ""),,drop=F], aes(label=Anno), color="black", size=3, 
	segment.size=0.5, direction="y", , nudge_y=2, nudge_x=0.05, hjust=0)+
	geom_hline(yintercept=-log10(0.05), colour="gray30", linetype="dashed", linewidth=0.5)+
	geom_vline(xintercept=c(-30, 30), colour="gray30", linetype="dashed", linewidth=0.5)+
	theme(plot.title=element_text(size=title_size, hjust=0.5, colour="black"), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	#legend.box.background=element_rect(linewidth=0.35, color="black"),
	legend.title=element_text(size=title_size, colour="black"), legend.text=element_text(size=text_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	legend.position=c(0.28, 0.83))+tag_thm

reads_info <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_len.tsv", h=T)
reads_info$type <- factor(reads_info$type, levels=types)
reads_info$len2 <- reads_info$utr3 + reads_info$utr5 + reads_info$body
reads_info <- reads_info[which(reads_info$utr3 > 10 & reads_info$utr5 > 10 & reads_info$body > 10),]
pbb <- ggplot(reads_info, aes(x=len2, color=type)) + stat_ecdf(linewidth=0.5)+
	labs(title=NULL, x="Length of isoforms", y="ECDF", color=NULL)+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	scale_x_continuous(limits=c(0, 2500), breaks=c(0, 500, 1500, 2500), expand=c(0, 0))+
	#annotate("text", x=100, y=1, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	#color="black", size=4, hjust=0, vjust=1)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	legend.position="none")+tag_thm
pbc <- ggplot(reads_info, aes(x=utr3, color=type)) + stat_ecdf(linewidth=0.5)+
	labs(title=NULL, x="Length of 3'UTRs", y=NULL, color=NULL)+
	scale_x_continuous(limits=c(0, 1500), expand=c(0, 0))+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	#annotate("text", x=50, y=1, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	#color="black", size=4, hjust=0, vjust=1)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), legend.key.size=unit(12, "pt"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	legend.position="none")+tag_thm
pbd <- ggplot(reads_info, aes(x=body, color=type)) + stat_ecdf(linewidth=0.5)+
	labs(title=NULL, x="Length of CDS", y=NULL, color=NULL)+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	scale_x_continuous(limits=c(0, 1500), breaks=c(0, 800, 1500), expand=c(0, 0))+
	#annotate("text", x=100, y=1, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	#color="black", size=4, hjust=0, vjust=1)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), legend.key.size=unit(12, "pt"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	legend.position="none")+tag_thm
pbe <- ggplot(reads_info, aes(x=utr5, color=type)) + stat_ecdf(linewidth=0.5)+
	labs(title=NULL, x="Length of 5'UTRs", y=NULL, color=NULL)+
	scale_x_continuous(limits=c(0, 1000), breaks=c(0, 400, 800), expand=c(0, 0))+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	#annotate("text", x=50, y=1, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	#color="black", size=4, hjust=0, vjust=1)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), legend.key.size=unit(12, "pt"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	legend.position=c(0.7, 0.23))+tag_thm

#gene <- "Cenpq"
#gene <- "Micos13"
#gene <- "Pigp"
#lims <- c(95, 140, 250)
#gene <- "Npc2"
#lims <- c(95, 140, 250)
gene <- "Scg2"
lims <- c(115, 115, 250)
term <- deg_info$Term[which(deg_info$Symbol == gene)]
trans <- scg[["features"]]$Trans[match(term, scg[["features"]]$Name)]
gene <- scg[["features"]]$Symbol[match(term, scg[["features"]]$Name)]
term_data <- data.frame(Cell=colnames(mat_utr5), Type=osn_rna$cell.subtype_fix, Length=mat_utr5[term,])
#term_data <- term_data[which(term_data$Length > max(mean(term_data$Length)*2/3, mean(term_data$Length)-100)),]
term_data <- term_data[which(term_data$Length > 0),]
term_data$Group <- "Early"
term_data$Group[which(term_data$Type == "mOSN" | term_data$Type == "iOSN")] <- "Later"
term_data$Type <- factor(term_data$Type, levels=types)
term_data$Group <- factor(term_data$Group, levels=c("Early", "Later"))
term_data$LG <- "M"
term_data$LG[which(term_data$Length < lims[1])] <- "S"
term_data$LG[which(term_data$Length > lims[2])] <- "L"
rec <- data.frame()
for (type in c("Early", "Later")) rec <- rbind(rec, data.frame(Type=type, 
	S=length(which(term_data$Length < lims[1] & term_data$Group == type)), 
	M=length(which(term_data$Length >= lims[1] & term_data$Length <= lims[2] & term_data$Group == type)), 
	L=length(which(term_data$Length > lims[2] & term_data$Group == type))))
rec$RS <- rec$S*100/(rec$S+rec$M+rec$L)
rec$RL <- rec$L*100/(rec$S+rec$M+rec$L)
rec_caa <- rbind(data.frame(Type=rec$Type, Group="RS", Count=rec$RS), 
	data.frame(Type=rec$Type, Group="RL", Count=rec$RL))
rec_caa$Type <- factor(rec_caa$Type, levels=c("Early", "Later"))
rec_caa$Group <- factor(rec_caa$Group, levels=c("RL", "RS"))
pcaa <- ggplot(rec_caa, aes(x=Type, y=Count, fill=Group))+geom_bar(stat="identity", width=0.6)+
	labs(title=NULL, x=gene, y="Percentage (%)", fill=NULL)+
	scale_fill_manual(values=col_list[c(4,3)], drop=F)+
	scale_y_continuous(breaks=c(25, 50, 75), expand=c(0, 0))+guides(fill=guide_legend(nrow=1))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), legend.position="none", axis.ticks.x=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title.y=element_text(size=title_size, colour="black"), 
	axis.title.x=element_text(size=title_size, colour="black", hjust=1.3), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks.y=element_line(linewidth=0.35, color="black"))+tag_thm
dens <- data.frame(density(term_data$Length)[c("x", "y")])
dens <- dens[which(dens$x > 0 & dens$x < lims[3]),]
dens$g <- "M"
dens$g[which(dens$x < lims[1])] <- "S" 
dens$g[which(dens$x > lims[2])] <- "L" 
#ss <- which(dens$g == "M")
#ss <- dens[c(ss[1], ss[length(ss)]),]
#ss$g <- c("S", "L")
#dens <- rbind(dens, ss)
dens$g <- factor(dens$g, levels=c("S", "L"), labels=c("Short", "Long")) 
pcab <- ggplot(dens, aes(x=x, y=y, color=g, fill=g))+geom_area(alpha=1)+
	labs(title=NULL, x="Length of 3'UTR (bp)", y=NULL, color=NULL, fill=NULL)+
	scale_color_manual(values=col_list[c(3,4)], drop=F)+
	scale_fill_manual(values=col_list[c(3,4)], drop=F)+
	guides(color=guide_legend(reverse=T), fill=guide_legend(reverse=T))+
	scale_x_continuous(breaks=seq(0, 250, 50))+
	scale_y_continuous(breaks=c((max(dens$y)+min(dens$y))/2), labels=c("Distribution"), expand=c(0, 0))+coord_flip()+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), legend.key.size=unit(12, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), legend.text=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black", angle=90, hjust=0.5, vjust=0.5), 
	axis.title=element_text(size=title_size, colour="black"), 
	axis.ticks.y=element_line(linewidth=0.35, color="black"), axis.ticks.x=element_blank(), 
	axis.line.y=element_line(linewidth=0.35, color="black"), axis.line.x=element_blank(), 
	legend.position=c(0.6, 0.8))+tag_thm
pca <- wrap_elements(wrap_plots(list(pcaa, pcab), nrow=1)+tag_thm)+tag_thm
#pca<- ggplot(term_data, aes(x=type, y=Length, fill=Group, color=Group))+geom_violin()+
#	geom_boxplot(width=0.2, outlier.alpha=0, fill="white")+scale_y_continuous()+
#	scale_x_discrete(drop=F)+scale_fill_manual(values=brewer.pal(9,"YlGnBu")[c(6,9)], drop=F)+
#	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(6,9)], drop=F)+
#	#annotate("text", x=0.5, y=1500, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
#	#color="black", size=4, hjust=0, vjust=1)+
#	labs(title=NULL, x="3'UTR length of Scg2", y="Length of 5'UTR (bp)")+
#	theme(axis.title=element_text(size=title_size, colour="black"), axis.text=element_text(size=text_size, colour="black"), 
#	axis.line=element_line(colour="black"), legend.position="none")+tag_thm

ego <- read.csv("num_go.csv", r=1, h=T)
ego <- ego[c(1,2,4:11),]
ego$Rank <- factor(rev(1:nrow(ego)))
ego$Type <- "Genes with increased isoforms"
rec_aa <- ego
paa <- ggplot(rec_aa, aes(x=Count, y=Rank, color=pvalue, size=Count))+geom_point()+
	scale_color_gradient(low="#440255", high="#FFFFBF")+
	#scale_color_continuous(low=brewer.pal(11,"Spectral")[1], high=brewer.pal(11,"Spectral")[6], 
	#guide=guide_colorbar(reverse=TRUE))+
	scale_size_continuous(breaks=c(25,28,30,32,35))+
	scale_y_discrete(breaks=rec_aa$Rank, labels=rec_aa$Description, position="right")+
	facet_grid(Type~., scales="free_y", space="free_y", switch="y")+
	guides(color=guide_colorbar(order=1), size=guide_legend(order=0))+
	labs(title=NULL, x=NULL, y=NULL, color="p.val")+
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), legend.key.size=unit(12, "pt"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(colour=col_list[2], fill=col_list[2]), 
	strip.text=element_text(size=title_size, color="white", face="bold"), 
	panel.spacing=unit(5, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	axis.line=element_blank(), panel.background=element_rect(fill='gray98'))+tag_thm

ego <- read.csv("test_go4.csv", r=1, h=T)
ego <- ego[1:10,]
ego$Rank <- factor(rev(1:nrow(ego)))
ego$Type <- "High expression novel isoforms"
rec_ab <- ego
pab <- ggplot(rec_ab, aes(x=Count, y=Rank, color=pvalue, size=Count))+geom_point()+
	scale_color_gradient(low="#440255", high="#FFFFBF")+
	#scale_color_continuous(low=brewer.pal(11,"Spectral")[1], high=brewer.pal(11,"Spectral")[6], 
	#guide=guide_colorbar(reverse=TRUE))+
	scale_size_continuous(breaks=c(5,7,9,11))+
	scale_y_discrete(breaks=rec_ab$Rank, labels=rec_ab$Description, position="right")+
	facet_grid(Type~., scales="free_y", space="free_y", switch="y")+
	guides(color=guide_colorbar(order=1), size=guide_legend(order=0))+
	labs(title=NULL, x=NULL, y=NULL, color="p.val")+
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), legend.key.size=unit(12, "pt"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(color=col_list[2], fill=col_list[2]), 
	strip.text=element_text(size=title_size, colour="white", face="bold"), 
	panel.spacing=unit(5, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	axis.line=element_blank(), panel.background=element_rect(fill='gray98'))+tag_thm

ego_info <- read.csv("cmp_raw_go2.csv", r=1, h=T)
ego_total <- data.frame()
for (j in 1:4)
{
	ids <- which(ego_info$Type == j)
	if (length(ids) < 1) next
	ego_total <- rbind(ego_total, ego_info[ids[1:min(8, length(ids))],])
}
ego_total$Rank <- factor(rev(1:nrow(ego_total)))
ego_total$Type <- factor(ego_total$Type, levels=4:1, labels=rev(c("HBC>GBC", "GBC>INP", "INP>iOSN", "iOSN>mOSN")))
rec_bca <- ego_total
cols <- rev(colorRampPalette(brewer.pal(9,"YlOrRd")[c(2,5,9)])(4))
pcb <- ggplot(rec_bca, aes(x=Count, y=Rank, color=pvalue, size=Count))+geom_point()+
	scale_color_gradient(low="#440255", high="#FFFFBF")+
	scale_y_discrete(breaks=rec_bca$Rank, labels=rec_bca$Description, position="right")+
	facet_grid2(Type~., scales="free_y", space="free_y", switch="y", 
	strip=strip_themed(background_y=elem_list_rect(color=cols, fill=cols)))+
	guides(color=guide_colorbar(order=1), size=guide_legend(order=0))+
	labs(title=NULL, x=NULL, y=NULL, color="p.val")+
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), legend.key.size=unit(12, "pt"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.text=element_text(size=title_size, colour="white", face="bold"), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	axis.line=element_blank(), panel.background=element_rect(fill='gray98'))+tag_thm

pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(paa,pab), nrow=1, widths=c(1,1))+
	plot_annotation(tag_levels=list(c("A", "B")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pba,pbb,pbc,pbd,pbe), nrow=1, widths=c(1,0.5,0.5,0.5,0.5))+
	plot_annotation(tag_levels=list(c("C", "E", "", "", "")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pca, pcb), nrow=1, widths=c(6,1))+
	plot_annotation(tag_levels=list(c("D", "F")), theme=tag_thm))+tag_thm, 
	pblank), ncol=1, heights=c(1,1,1,1.6)), 
	width=210, height=297, dpi=300, units="mm", filename="oe_fl_fig_S02.png", limitsize=F)

ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(paa,pab), nrow=1, widths=c(1,1))+
	plot_annotation(tag_levels=list(c("A", "B")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pba,pbb,pbc,pbd,pbe), nrow=1, widths=c(1,0.5,0.5,0.5,0.5))+
	plot_annotation(tag_levels=list(c("C", "E", "", "", "")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pca, pcb), nrow=1, widths=c(6,1))+
	plot_annotation(tag_levels=list(c("D", "F")), theme=tag_thm))+tag_thm, 
	pblank), ncol=1, heights=c(1,1,1,1.6)), 
	width=210, height=297, dpi=300, units="mm", filename="oe_fl_fig_S02.pdf", limitsize=F)

