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
text_size <- 13
title_size <- 15
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, colour="black"), plot.margin=margin(-3,0,3,0), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", colour=NA),  plot.background=element_rect(fill="transparent", colour=NA), 
	legend.box.spacing=unit(0, "pt"))

#gtf <- read.table("genes.gtf", sep="\t")
#gtf <- gtf[grep("^chr.*", gtf[, 1]),]
#colnames(gtf) <- c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
#gtf$transcript_id <- gsub(";.*", "", gsub(".*transcript_id ", "", gtf$attributes))
#gtf$gene_id <- gsub(";.*", "", gsub(".*gene_id ", "", gtf$attributes))
#gtf$gene_name <- gsub(";.*", "", gsub(".*gene_name ", "", gtf$attributes))
gtf <- read.delim("gtf_info.tsv")

types <- c("HBC", "GBC", "INP", "iOSN", "mOSN")
osn_sct <- readRDS("osn_sct.rds")
osn_rna <- readRDS("osn_rna.rds")
osn_rna$nIsoform <- as.numeric(colSums(osn_sct[["trans"]] > 0))
osn_rna$nGene<- as.numeric(colSums(osn_sct[["genes"]] > 0))
#rec <- data.frame(bc=colnames(osn_rna), type=osn_rna$cell.subtype_fix)
#write.table(rec, "/mnt/md0/oe_full_length/output/OEfulllength/neuron_bc_info.tsv", col.names=F, row.names=F, quote=F, sep="\t")

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

reads_info <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_len.tsv", h=T)
reads_info$type <- factor(reads_info$type, levels=types)
reads_info_aa <- reads_info[which(reads_info$len > 10),]
r = summary(aov(len~type, data=reads_info_aa))
paa <- ggplot(reads_info_aa, aes(x=body, color=type)) + stat_ecdf(linewidth=0.8)+
	labs(title=NULL, x="Length of CDS", y="ECDF", color=NULL)+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	scale_x_continuous(limits=c(0, 3000), expand=c(0, 0))+
	#annotate("text", x=100, y=1, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	#color="black", size=4, hjust=0, vjust=1)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position = c(0.82, 0.23))+tag_thm

reads_info <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_len.tsv", h=T)
reads_info$type <- factor(reads_info$type, levels=types)
reads_info_ab <- reads_info[which(reads_info$utr5 > 10 & reads_info$body > 10),]
r = summary(aov(len~type, data=reads_info_ab))
pab <- ggplot(reads_info_ab, aes(x=utr5, color=type)) + stat_ecdf(linewidth=0.8)+
	labs(title=NULL, x="Length of 5'UTRs", y="ECDF", color=NULL)+
	scale_x_continuous(limits=c(0, 1500), expand=c(0, 0))+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	#annotate("text", x=50, y=1, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	#color="black", size=4, hjust=0, vjust=1)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position = c(0.82, 0.23))+tag_thm

sce_trans <- CreateSeuratObject(scg[["matrix"]], project="OEnano")
sce_trans$cell.subtype_fix <- as.character(osn_rna$cell.subtype_fix[match(colnames(sce_trans), colnames(osn_rna))])
sce_trans <- SCTransform(sce_trans, method="glmGamPoi")
reads_info <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_len.tsv", h=T)
reads_info$type <- factor(reads_info$type, levels=types)
rec <- data.frame()
for (type in types[-4])
{
	terms <- FindMarkers(sce_trans, ident.1=type, group.by="cell.subtype_fix")
	terms <- rownames(terms)[which(terms$p_val_adj < 0.01 & terms$avg_log2FC > 2)]
	rec <- rbind(rec, reads_info[match(paste0(type, ",", terms), paste0(reads_info$type, ",", gsub("_", "-", reads_info$trans))),])
}
rec_ac <- rec[which(rec$utr5 > 10 & rec$body > 10 & rec$utr5 < 2000),]
rec_ac$type <- factor(rec_ac$type, levels=types[-4], labels=c(types[1:3], "OSN"))
r = summary(aov(utr5~type, data=rec_ac))
pac <- ggplot(rec_ac, aes(x=type, y=utr5, fill=type, color=type))+geom_violin()+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white")+scale_y_continuous()+
	scale_x_discrete(drop=F)+scale_fill_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,9)], drop=F)+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,9)], drop=F)+
	#annotate("text", x=0.5, y=900, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	#color="black", size=4, hjust=0, vjust=1)+
	labs(title=NULL, x="5'UTR length of\nmarker isoforms", y="Length of 5'UTR")+
	theme(axis.title=element_text(size=title_size, colour="black"), axis.text=element_text(size=text_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position="none")+tag_thm

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
deg_info$Group[which(deg_info$Padj < 0.05 & deg_info$Diff < -30)] <- "Decreased"
#deg_info$X <- log2(abs(deg_info$Diff))*deg_info$Diff/abs(deg_info$Diff)
deg_info$X <- deg_info$Diff
deg_info$Y <- -log10(deg_info$Padj)
deg_info$Rank <- abs(deg_info$X)
deg_info <- deg_info[order(deg_info$Rank, decreasing=T),]
deg_info$Anno <- ""
#deg_info$Anno[which(deg_info$Symbol == "Calm1")] <- "Calm1"
#deg_info$Anno <- deg_info$Symbol
#deg_info$Anno[which(deg_info$Group == "N.S.")] <- ""
#deg_info$Anno[which(deg_info$Anno != "")[-c(1:10)]] <- ""
deg_info$Group <- factor(deg_info$Group, levels=c("Increased", "Decreased", "N.S."))
deg_info <- deg_info[which(abs(deg_info$Diff) > 10),]
pad <- ggplot(deg_info, aes(x=X, y=Y, color=Group))+geom_point(size=2)+
	labs(title=NULL, x="Difference of 5'UTR length", y="P.adj (-log10)", color=NULL)+
	scale_color_manual(values=c("#FC8D62", "#8DA0CB", "gray80"), drop=F)+scale_y_continuous(expand=c(0, 0))+
	#scale_x_continuous(breaks=seq(-90, 100, 30))+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	#geom_text(aes(y=Y+1, label=Anno), size=3, col="black", vjust=0)+
	geom_text_repel(data=deg_info[which(deg_info$Anno != ""),,drop=F], aes(label=Anno), color="black", size=4, 
	segment.size=0.5, direction="y", , nudge_y=2, nudge_x=0.05, hjust=0)+
	geom_hline(yintercept=-log10(0.05), colour="gray30", linetype="dashed", linewidth=0.5)+
	geom_vline(xintercept=c(-30, 30), colour="gray30", linetype="dashed", linewidth=0.5)+
	theme(plot.title=element_text(size=title_size, hjust=0.5, colour="black"), panel.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), legend.text=element_text(size=text_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"))+tag_thm

ego <- read.csv("num_go.csv", r=1, h=T)
ego <- ego[c(1,2,4:11),]
ego$Rank <- factor(rev(1:nrow(ego)))
ego$Type <- "Genes with increased isoforms"
rec_ba <- ego
pba <- ggplot(rec_ba, aes(x=Count, y=Rank, color=pvalue, size=Count))+geom_point()+
	scale_color_gradient(low="#440255", high="#FFFFBF")+
	#scale_color_continuous(low=brewer.pal(11,"Spectral")[1], high=brewer.pal(11,"Spectral")[6], 
	#guide=guide_colorbar(reverse=TRUE))+
	scale_y_discrete(breaks=rec_ba$Rank, labels=rec_ba$Description, position="right")+
	facet_grid(Type~., scales="free_y", space="free_y", switch="y")+
	guides(color=guide_colorbar(order=1), size=guide_legend(order=0))+
	labs(title=NULL, x=NULL, y=NULL, color="p.val")+
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(colour=col_list[2], fill=col_list[2]), 
	strip.text=element_text(size=text_size, color="white", face="bold"), 
	panel.spacing=unit(5, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	axis.line=element_blank(), panel.background=element_rect(fill='gray98'))+tag_thm

ego <- read.csv("test_go4.csv", r=1, h=T)
ego <- ego[1:10,]
ego$Rank <- factor(rev(1:nrow(ego)))
ego$Type <- "High expression novel isoforms"
rec_bb <- ego
pbb <- ggplot(rec_bb, aes(x=Count, y=Rank, color=pvalue, size=Count))+geom_point()+
	scale_color_gradient(low="#440255", high="#FFFFBF")+
	#scale_color_continuous(low=brewer.pal(11,"Spectral")[1], high=brewer.pal(11,"Spectral")[6], 
	#guide=guide_colorbar(reverse=TRUE))+
	scale_y_discrete(breaks=rec_bb$Rank, labels=rec_bb$Description, position="right")+
	facet_grid(Type~., scales="free_y", space="free_y", switch="y")+
	guides(color=guide_colorbar(order=1), size=guide_legend(order=0))+
	labs(title=NULL, x=NULL, y=NULL, color="p.val")+
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(color=col_list[2], fill=col_list[2]), 
	strip.text=element_text(size=text_size, colour="white", face="bold"), 
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
pca <- ggplot(rec_bca, aes(x=Count, y=Rank, color=pvalue, size=Count))+geom_point()+
	scale_color_gradient(low="#440255", high="#FFFFBF")+
	scale_y_discrete(breaks=rec_bca$Rank, labels=rec_bca$Description, position="right")+
	facet_grid2(Type~., scales="free_y", space="free_y", switch="y", 
	strip=strip_themed(background_y=elem_list_rect(color=cols, fill=cols)))+
	guides(color=guide_colorbar(order=1), size=guide_legend(order=0))+
	labs(title=NULL, x=NULL, y=NULL, color="p.val")+
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.text=element_text(size=text_size, colour="white", face="bold"), 
	panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	axis.line=element_blank(), panel.background=element_rect(fill='gray98'))+tag_thm

pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(paa,pab,pac,pad), nrow=1, widths=c(0.6,0.6,1,1))+
	plot_annotation(tag_levels=list(c("A", "B", "C", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pba,pbb), nrow=1, widths=c(1,1))+
	plot_annotation(tag_levels=list(c("E", "F")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pca, pblank), nrow=1, widths=c(1,6))+
	plot_annotation(tag_levels=list(c("G")), theme=tag_thm))+tag_thm), 
	ncol=1, heights=c(1,1,1)), width=13, height=12, dpi=200, filename="oe_fl_fig_S02.png", limitsize=F)


