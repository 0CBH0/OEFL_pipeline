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
options(stringsAsFactors=FALSE)

col_list <- c(c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8"), brewer.pal(12,"Set3")[-c(2, 9)])
text_size <- 13
title_size <- 15
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, color="black"), plot.margin=margin(0,-3,0,-3), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))

#gtf <- read.table("genes.gtf", sep="\t")
#gtf <- gtf[grep("^chr.*", gtf[, 1]),]
#colnames(gtf) <- c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
#gtf$transcript_id <- gsub(";.*", "", gsub(".*transcript_id ", "", gtf$attributes))
#gtf$gene_id <- gsub(";.*", "", gsub(".*gene_id ", "", gtf$attributes))
#gtf$gene_name <- gsub(";.*", "", gsub(".*gene_name ", "", gtf$attributes))
gtf <- read.delim("gtf_info.tsv")

spatial_genes <- c("Nqo1", "Acsm4", "Ncam2", "Nfix", "Nfib", "Plxna1", "Sema3a", "Nrp1")
types <- c("HBC", "GBC", "INP", "iOSN", "mOSN")
osn_rna <- readRDS("osn_rna2_fix.rds")
#osn_rna$sp <- colSums(osn_rna[["SCT"]]$scale.data[sp_genes,] > 0)
#osn_rna <- subset(osn_rna, cells=colnames(osn_rna)[which(osn_rna$cell.subtype_fix == "Mature")])
#osn_rna <- SCTransform(osn_rna, method="glmGamPoi")
#osn_rna <- RunUMAP(osn_rna, dims=1:6)

genes_d <- intersect(read.table("25_dorsal_features.txt")[, 1], rownames(osn_rna[["SCT"]]$scale.data))
genes_v <- intersect(read.table("48_ventral_features.txt")[, 1], rownames(osn_rna[["SCT"]]$scale.data))
genes_a <- intersect(read.table("50_anterior_features.txt")[, 1], rownames(osn_rna[["SCT"]]$scale.data))
genes_p <- intersect(read.table("57_posterior_features.txt")[, 1], rownames(osn_rna[["SCT"]]$scale.data))
genes_d <- genes_d[which(cor(t(osn_rna[["SCT"]]$scale.data[genes_d,]), osn_rna@reductions[["umap"]]@cell.embeddings[, 1]) < -0.1)]
genes_v <- genes_v[which(cor(t(osn_rna[["SCT"]]$scale.data[genes_v,]), osn_rna@reductions[["umap"]]@cell.embeddings[, 1]) > 0.1)]
genes_a <- genes_a[which(cor(t(osn_rna[["SCT"]]$scale.data[genes_a,]), osn_rna@reductions[["umap"]]@cell.embeddings[, 2]) > 0.1)]
genes_p <- genes_p[which(cor(t(osn_rna[["SCT"]]$scale.data[genes_p,]), osn_rna@reductions[["umap"]]@cell.embeddings[, 2]) < -0.1)]
sp_genes <- union(union(union(genes_d, genes_v), genes_a), genes_p)

gene_mat <- as.matrix(osn_rna[["SCT"]]$scale.data[genes_d,])
for (i in 1:nrow(gene_mat))
{
	gene_mat[i,] <- gene_mat[i,] - min(gene_mat[i,])
	gene_mat[i,] <- gene_mat[i,] * 2 / max(gene_mat[i,]) - 1
}
osn_rna$vds <- colSums(gene_mat)
osn_rna$vds <- osn_rna$vds - min(osn_rna$vds)
osn_rna$vds <- osn_rna$vds * 2 / max(osn_rna$vds) - 1

gene_mat <- as.matrix(osn_rna[["SCT"]]$scale.data[genes_v,])
for (i in 1:nrow(gene_mat))
{
	gene_mat[i,] <- gene_mat[i,] - min(gene_mat[i,])
	gene_mat[i,] <- gene_mat[i,] * 2 / max(gene_mat[i,]) - 1
}
osn_rna$vvs <- colSums(gene_mat)
osn_rna$vvs <- osn_rna$vvs - min(osn_rna$vvs)
osn_rna$vvs <- osn_rna$vvs * 2 / max(osn_rna$vvs) - 1

gene_mat <- as.matrix(osn_rna[["SCT"]]$scale.data[genes_a,])
for (i in 1:nrow(gene_mat))
{
	gene_mat[i,] <- gene_mat[i,] - min(gene_mat[i,])
	gene_mat[i,] <- gene_mat[i,] * 2 / max(gene_mat[i,]) - 1
}
osn_rna$vas <- colSums(gene_mat)
osn_rna$vas <- osn_rna$vas - min(osn_rna$vas)
osn_rna$vas <- osn_rna$vas * 2 / max(osn_rna$vas) - 1

gene_mat <- as.matrix(osn_rna[["SCT"]]$scale.data[genes_p,])
for (i in 1:nrow(gene_mat))
{
	gene_mat[i,] <- gene_mat[i,] - min(gene_mat[i,])
	gene_mat[i,] <- gene_mat[i,] * 2 / max(gene_mat[i,]) - 1
}
osn_rna$vps <- colSums(gene_mat)
osn_rna$vps <- osn_rna$vps - min(osn_rna$vps)
osn_rna$vps <- osn_rna$vps * 2 / max(osn_rna$vps) - 1

osn_rna$vdv <- (osn_rna$vvs - osn_rna$vds) / 2
osn_rna$vap <- (osn_rna$vps - osn_rna$vas) / 2
rec <- data.frame(BC=colnames(osn_rna), Group="S")
rec$Group[which(osn_rna$vdv <= quantile(osn_rna$vdv, 0.3))] <- "D"
rec$Group[which(osn_rna$vdv >= quantile(osn_rna$vdv, 0.7))] <- "V"
rec <- rec[which(rec$Group != "S"),]
write.csv(rec, "barcode_osn_dv2.csv", row.names=F, quote=F)
rec <- data.frame(BC=colnames(osn_rna), Group="S")
rec$Group[which(osn_rna$vap <= quantile(osn_rna$vap, 0.3))] <- "A"
rec$Group[which(osn_rna$vap >= quantile(osn_rna$vap, 0.7))] <- "P"
rec <- rec[which(rec$Group != "S"),]
write.csv(rec, "barcode_osn_ap2.csv", row.names=F, quote=F)

data_raw <- H5File$new("/mnt/md0/oe_full_length/output/OSNfulllength/osn_trans_ass.h5", mode="r")
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
ids <- match(colnames(osn_rna), colnames(osn_sct_raw[["trans"]]))
osn_sct_raw[["trans"]] <- osn_sct_raw[["trans"]][, ids]
osn_sct_raw[["utr3"]] <- osn_sct_raw[["utr3"]][, ids]
osn_sct_raw[["utr5"]] <- osn_sct_raw[["utr5"]][, ids]

pa <- wrap_elements(ggdraw()+draw_image("fig03_a.png", scale=1.2))+tag_thm
termb <- rbind(rbind(data.frame(Term="A.score", Val=osn_rna$vas), data.frame(Term="P.score", Val=-osn_rna$vps)), 
	rbind(data.frame(Term="D.score", Val=osn_rna$vds), data.frame(Term="V.score", Val=-osn_rna$vvs)))
termb$Term <- factor(termb$Term, levels=c("D.score", "V.score", "A.score", "P.score"))
pb <- ggplot(termb, aes(x=Val, fill=Term, color=Term))+geom_density(linewidth=0.8, alpha=0.02)+
	labs(title="Spatial score", x="Score", y="Density", color=NULL, fill=NULL)+
	scale_fill_manual(values=col_list[c(1,2,3,6)])+scale_color_manual(values=col_list[c(1,2,3,6)])+
	scale_y_continuous(limits=c(0, 3.4), expand=c(0, 0))+
	theme(axis.line=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=0.8), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black", face="plain"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), 
	panel.background=element_rect(0, linetype=0), axis.text=element_text(size=text_size, color="black"), 
	axis.title=element_text(size=title_size, color="black"), legend.position=c(0.76, 0.76))+tag_thm

osn_rna_sub <- subset(osn_rna, cells=colnames(osn_rna)[which(osn_rna$vdv != 0)])
rec_ac <- data.frame(osn_rna_sub@meta.data)
rec_ac$X <- osn_rna_sub@reductions[["umap"]]@cell.embeddings[, 1]
rec_ac$Y <- osn_rna_sub@reductions[["umap"]]@cell.embeddings[, 2]
pc <- ggplot(rec_ac, aes(x=X, y=Y, color=-vdv))+geom_point(size=1)+
	labs(title="Score of D > V", x="UMAP1", y="UMAP2", color="Score")+
	scale_color_gradientn(colors=rev(brewer.pal(11, "Spectral")))+
	theme(axis.line=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=0.8), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black", face="plain"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), 
	panel.background=element_rect(0, linetype=0), axis.text=element_blank(), axis.ticks=element_blank(), 
	axis.title=element_text(size=title_size, color="black"))+tag_thm
osn_rna_sub <- subset(osn_rna, cells=colnames(osn_rna)[which(osn_rna$vap != 0)])
rec_ad <- data.frame(osn_rna_sub@meta.data)
rec_ad$X <- osn_rna_sub@reductions[["umap"]]@cell.embeddings[, 1]
rec_ad$Y <- osn_rna_sub@reductions[["umap"]]@cell.embeddings[, 2]
pd <- ggplot(rec_ad, aes(x=X, y=Y, color=-vap))+geom_point(size=1)+
	labs(title="Score of A > P", x="UMAP1", y="UMAP2", color="Score")+
	scale_color_gradientn(colors=rev(brewer.pal(11, "Spectral")))+
	theme(axis.line=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=0.8), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black", face="plain"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), 
	panel.background=element_rect(0, linetype=0), axis.text=element_blank(), axis.ticks=element_blank(), 
	axis.title=element_text(size=title_size, color="black"))+tag_thm

terms <- names(which(table(osn_sct_raw[["features"]]$Gene) > 1))
bins <- 40
rec_total_counts <- data.frame()
rec_total_rate <- data.frame()
bin_step <- ncol(osn_sct_raw[["trans"]])/(bins+0.25)
bin_size <- 1.25*bin_step
rec_ppp <- data.frame(Cell=colnames(osn_rna), Value=osn_rna$vdv)
rec_ppp <- rec_ppp[order(rec_ppp$Value),]
rec_ppp$GS <- 1:nrow(rec_ppp)
rec_ppp$LH <- ""
rec_ppp$LH[1:round(nrow(rec_ppp)/4)] <- "L"
rec_ppp$LH[round(nrow(rec_ppp)*3/4):nrow(rec_ppp)] <- "H"
rec_ppp <- rec_ppp[match(colnames(osn_rna), rec_ppp$Cell),]
id_l <- match(rec_ppp$Cell[which(rec_ppp$LH == "L")], colnames(osn_sct_raw[["trans"]]))
id_h <- match(rec_ppp$Cell[which(rec_ppp$LH == "H")], colnames(osn_sct_raw[["trans"]]))
osn_sct_dv_lh <- data.frame(matrix(0, nrow=nrow(osn_sct_raw[["trans"]]), ncol=2, dimnames=list(rownames(osn_sct_raw[["trans"]]), c("L", "H"))))
osn_sct_dv_lh$L <- rowSums(osn_sct_raw[["trans"]][,which(rec_ppp$LH == "L")])
osn_sct_dv_lh$H <- rowSums(osn_sct_raw[["trans"]][,which(rec_ppp$LH == "H")])
rec_val_dv <- rep(0, bins)
osn_sct_dv <- data.frame(matrix(0, nrow=nrow(osn_sct_raw[["trans"]]), ncol=bins, dimnames=list(rownames(osn_sct_raw[["trans"]]), 1:bins)))
for (i in 1:bins)
{
	s <- round((i-1)*bin_step+1)
	e <- round((i-1)*bin_step+bin_size)
	if (i == bins) e <- ncol(osn_sct_raw[["trans"]])
	#print(paste(s, e, e - s + 1))
	ids <- which(rec_ppp$GS >= s & rec_ppp$GS <= e)
	rec_val_dv[i] <- mean(rec_ppp$Value[ids])
	osn_sct_dv[, i] <- rowSums(osn_sct_raw[["trans"]][,ids])
}
cl <- makeCluster(16, type="FORK")
results <- parLapply(cl, terms, function(x) {
	ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
	rec_sub <- data.frame(Gene=x, RNA=0, Num=0, Count=0, CM=0, Ctr=0, Pval=1, R2=0, Cor=0, CP=0, Diff=0, DS=0, TP=1, Trans="", Lv="", PV=0, SV="", SR=0)
	if (length(ids) < 2) return(rec_sub)
	mat <- as.matrix(osn_sct_dv[ids,])
	ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
	if (length(ts) < 2) return(rec_sub)
	for (i in ids) rec_sub$TP[1] <- min(rec_sub$TP[1], t.test(osn_sct_raw[["trans"]][i, id_l], osn_sct_raw[["trans"]][i, id_h])$p.value)
	rec_sub$RNA[1] <- sum(mat)/ncol(osn_rna)
	mat <- mat[ts,]
	rec_sub$Num[1] <- nrow(mat)
	rec_sub$Count[1] <- min(rowMaxs(mat))
	rec_sub$CM[1] <- max(mat)
	rec_sub$Ctr[1] <- min(rowSums(mat))*100/max(rowSums(mat))
	for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
	rk <- abs(mat[, 1] - mat[, bins])
	rec_sub$Diff[1] <- max(rk)
	rk <- rk[order(rk, decreasing=T)]
	rec_sub$Trans[1] <- paste(names(rk)[1:2], collapse=";")
	pp <- osn_sct_dv_lh[ids,]
	ppc <- pp
	for (i in 1:ncol(pp)) pp[, i] <- pp[, i]*100/sum(pp[, i])
	rec_sub$Lv[1] <- paste(round(pp[names(rk)[1:2],2]-pp[names(rk)[1:2],1]), collapse=";")
	rec_sub$PV[1] <- max(abs(pp[names(rk)[1:2],2]-pp[names(rk)[1:2],1]))
	ia <- round(pp[names(rk)[1:2],1])
	ib <- round(pp[names(rk)[1:2],2])
	ia <- as.numeric(round(ia[1]*100/sum(ia), 2))
	ib <- as.numeric(round(ib[1]*100/sum(ib), 2))
	rec_sub$SV[1] <- paste0(ia, ";", ib)
	rec_sub$SR[1] <- round(ia/ib, 2)
	rec_sub$Pval[1] <- fisher.test(pp[names(rk)[1:2],])$p.value
	t <- rep(0, nrow(mat))
	for (k in 1:nrow(mat)) t[k] <- abs(summary(lm(mat[k,]~rec_val_dv))$adj.r.squared)
	rec_sub$R2[1] <- max(t)
	rec_sub$DS[1] <- 0
	ra <- cor.test(mat[names(rk)[1],], rec_val_dv, method="spearman")
	rb <- cor.test(mat[names(rk)[2],], rec_val_dv, method="spearman")
	rec_sub$Cor[1] <- max(abs(ra$estimate), abs(rb$estimate))
	rec_sub$CP[1] <- min(ra$p.value, rb$p.value)
	return(rec_sub)
})
stopCluster(cl)
rec <- do.call("rbind", results)
rec$p.adj <- p.adjust(rec$Pval, "BH")
rec$cp.adj <- p.adjust(rec$CP, "BH")
write.csv(rec, "cmp_vdv_raw2.csv")

bin_step <- ncol(osn_sct_raw[["trans"]])/(bins+0.25)
bin_size <- 1.25*bin_step
rec_ppp <- data.frame(Cell=colnames(osn_rna), Value=osn_rna$vap)
rec_ppp <- rec_ppp[order(rec_ppp$Value),]
rec_ppp$GS <- 1:nrow(rec_ppp)
rec_ppp$LH <- ""
rec_ppp$LH[1:round(nrow(rec_ppp)/4)] <- "L"
rec_ppp$LH[round(nrow(rec_ppp)*3/4):nrow(rec_ppp)] <- "H"
rec_ppp <- rec_ppp[match(colnames(osn_rna), rec_ppp$Cell),]
id_l <- match(rec_ppp$Cell[which(rec_ppp$LH == "L")], colnames(osn_sct_raw[["trans"]]))
id_h <- match(rec_ppp$Cell[which(rec_ppp$LH == "H")], colnames(osn_sct_raw[["trans"]]))
osn_sct_ap_lh <- data.frame(matrix(0, nrow=nrow(osn_sct_raw[["trans"]]), ncol=2, dimnames=list(rownames(osn_sct_raw[["trans"]]), c("L", "H"))))
osn_sct_ap_lh$L <- rowSums(osn_sct_raw[["trans"]][,which(rec_ppp$LH == "L")])
osn_sct_ap_lh$H <- rowSums(osn_sct_raw[["trans"]][,which(rec_ppp$LH == "H")])
rec_val_ap <- rep(0, bins)
osn_sct_ap <- data.frame(matrix(0, nrow=nrow(osn_sct_raw[["trans"]]), ncol=bins, dimnames=list(rownames(osn_sct_raw[["trans"]]), 1:bins)))
for (i in 1:bins)
{
	s <- round((i-1)*bin_step+1)
	e <- round((i-1)*bin_step+bin_size)
	if (i == bins) e <- ncol(osn_sct_raw[["trans"]])
	#print(paste(s, e, e - s + 1))
	ids <- which(rec_ppp$GS >= s & rec_ppp$GS <= e)
	rec_val_ap[i] <- mean(rec_ppp$Value[ids])
	osn_sct_ap[, i] <- rowSums(osn_sct_raw[["trans"]][,ids])
}
#x = "ENSMUSG00000028647"
#x = "ENSMUSG00000024158"
cl <- makeCluster(16, type="FORK")
results <- parLapply(cl, terms, function(x) {
	ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
	rec_sub <- data.frame(Gene=x, RNA=0, Num=0, Count=0, CM=0, Ctr=0, Pval=1, R2=0, Cor=0, CP=0, Diff=0, DS=0, TP=1, Trans="", Lv="", PV=0, SV="", SR=0)
	if (length(ids) < 2) return(rec_sub)
	mat <- as.matrix(osn_sct_ap[ids,])
	ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
	if (length(ts) < 2) return(rec_sub)
	for (i in ids) rec_sub$TP[1] <- min(rec_sub$TP[1], t.test(osn_sct_raw[["trans"]][i, id_l], osn_sct_raw[["trans"]][i, id_h])$p.value)
	rec_sub$RNA[1] <- sum(mat)/ncol(osn_rna)
	mat <- mat[ts,]
	rec_sub$Num[1] <- nrow(mat)
	rec_sub$Count[1] <- min(rowMaxs(mat))
	rec_sub$CM[1] <- max(mat)
	rec_sub$Ctr[1] <- min(rowSums(mat))*100/max(rowSums(mat))
	for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
	rk <- abs(mat[, 1] - mat[, bins])
	rec_sub$Diff[1] <- max(rk)
	rk <- rk[order(rk, decreasing=T)]
	rec_sub$Trans[1] <- paste(names(rk)[1:2], collapse=";")
	pp <- osn_sct_ap_lh[ids,]
	ppc <- pp
	for (i in 1:ncol(pp)) pp[, i] <- pp[, i]*100/sum(pp[, i])
	rec_sub$Lv[1] <- paste(round(pp[names(rk)[1:2],2]-pp[names(rk)[1:2],1]), collapse=";")
	rec_sub$PV[1] <- max(abs(pp[names(rk)[1:2],2]-pp[names(rk)[1:2],1]))
	ia <- round(pp[names(rk)[1:2],1])
	ib <- round(pp[names(rk)[1:2],2])
	ia <- as.numeric(round(ia[1]*100/sum(ia), 2))
	ib <- as.numeric(round(ib[1]*100/sum(ib), 2))
	rec_sub$SV[1] <- paste0(ia, ";", ib)
	rec_sub$SR[1] <- round(ia/ib, 2)
	rec_sub$Pval[1] <- fisher.test(ppc[names(rk)[1:2],])$p.value
	t <- rep(0, nrow(mat))
	for (k in 1:nrow(mat)) t[k] <- abs(summary(lm(mat[k,]~rec_val_ap))$adj.r.squared)
	rec_sub$R2[1] <- max(t)
	rec_sub$DS[1] <- 0
	ra <- cor.test(mat[names(rk)[1],], rec_val_ap, method="spearman")
	rb <- cor.test(mat[names(rk)[2],], rec_val_ap, method="spearman")
	rec_sub$Cor[1] <- max(abs(ra$estimate), abs(rb$estimate))
	rec_sub$CP[1] <- min(ra$p.value, rb$p.value)
	return(rec_sub)
})
stopCluster(cl)
rec <- do.call("rbind", results)
rec$p.adj <- p.adjust(rec$Pval, "BH")
rec$cp.adj <- p.adjust(rec$CP, "BH")
write.csv(rec, "cmp_vap_raw2.csv")

#osn_rna$vtt <- runif(ncol(osn_rna), -1, 1)
bin_step <- ncol(osn_sct_raw[["trans"]])/(bins+0.25)
bin_size <- 1.25*bin_step
rec_ppp <- data.frame(Cell=colnames(osn_rna), Value=osn_rna$vtt)
rec_ppp <- rec_ppp[order(rec_ppp$Value),]
rec_ppp$GS <- 1:nrow(rec_ppp)
rec_ppp$LH <- ""
rec_ppp$LH[1:round(nrow(rec_ppp)/4)] <- "L"
rec_ppp$LH[round(nrow(rec_ppp)*3/4):nrow(rec_ppp)] <- "H"
rec_ppp <- rec_ppp[match(colnames(osn_rna), rec_ppp$Cell),]
id_l <- match(rec_ppp$Cell[which(rec_ppp$LH == "L")], colnames(osn_sct_raw[["trans"]]))
id_h <- match(rec_ppp$Cell[which(rec_ppp$LH == "H")], colnames(osn_sct_raw[["trans"]]))
osn_sct_tt_lh <- data.frame(matrix(0, nrow=nrow(osn_sct_raw[["trans"]]), ncol=2, dimnames=list(rownames(osn_sct_raw[["trans"]]), c("L", "H"))))
osn_sct_tt_lh$L <- rowSums(osn_sct_raw[["trans"]][,which(rec_ppp$LH == "L")])
osn_sct_tt_lh$H <- rowSums(osn_sct_raw[["trans"]][,which(rec_ppp$LH == "H")])
rec_val_tt <- rep(0, bins)
osn_sct_tt <- data.frame(matrix(0, nrow=nrow(osn_sct_raw[["trans"]]), ncol=bins, dimnames=list(rownames(osn_sct_raw[["trans"]]), 1:bins)))
for (i in 1:bins)
{
	s <- round((i-1)*bin_step+1)
	e <- round((i-1)*bin_step+bin_size)
	if (i == bins) e <- ncol(osn_sct_raw[["trans"]])
	#print(paste(s, e, e - s + 1))
	ids <- which(rec_ppp$GS >= s & rec_ppp$GS <= e)
	rec_val_tt[i] <- mean(rec_ppp$Value[ids])
	osn_sct_tt[, i] <- rowSums(osn_sct_raw[["trans"]][,ids])
}
cl <- makeCluster(16, type="FORK")
results <- parLapply(cl, terms, function(x) {
	ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
	rec_sub <- data.frame(Gene=x, RNA=0, Num=0, Count=0, CM=0, Ctr=0, Pval=1, R2=0, Cor=0, CP=0, Diff=0, DS=0, TP=1, Trans="", Lv="", PV=0, SV="", SR=0)
	if (length(ids) < 2) return(rec_sub)
	mat <- as.matrix(osn_sct_tt[ids,])
	ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
	if (length(ts) < 2) return(rec_sub)
	for (i in ids) rec_sub$TP[1] <- min(rec_sub$TP[1], t.test(osn_sct_raw[["trans"]][i, id_l], osn_sct_raw[["trans"]][i, id_h])$p.value)
	rec_sub$RNA[1] <- sum(mat)/ncol(osn_rna)
	mat <- mat[ts,]
	rec_sub$Num[1] <- nrow(mat)
	rec_sub$Count[1] <- min(rowMaxs(mat))
	rec_sub$CM[1] <- max(mat)
	rec_sub$Ctr[1] <- min(rowSums(mat))*100/max(rowSums(mat))
	for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
	rk <- abs(mat[, 1] - mat[, bins])
	rec_sub$Diff[1] <- max(rk)
	rk <- rk[order(rk, decreasing=T)]
	rec_sub$Trans[1] <- paste(names(rk)[1:2], collapse=";")
	pp <- osn_sct_tt_lh[ids,]
	ppc <- pp
	for (i in 1:ncol(pp)) pp[, i] <- pp[, i]*100/sum(pp[, i])
	rec_sub$Lv[1] <- paste(round(pp[names(rk)[1:2],2]-pp[names(rk)[1:2],1]), collapse=";")
	rec_sub$PV[1] <- max(abs(pp[names(rk)[1:2],2]-pp[names(rk)[1:2],1]))
	ia <- round(pp[names(rk)[1:2],1])
	ib <- round(pp[names(rk)[1:2],2])
	ia <- as.numeric(round(ia[1]*100/sum(ia), 2))
	ib <- as.numeric(round(ib[1]*100/sum(ib), 2))
	rec_sub$SV[1] <- paste0(ia, ";", ib)
	rec_sub$SR[1] <- round(ia/ib, 2)
	rec_sub$Pval[1] <- fisher.test(pp[names(rk)[1:2],])$p.value
	t <- rep(0, nrow(mat))
	for (k in 1:nrow(mat)) t[k] <- abs(summary(lm(mat[k,]~rec_val_tt))$adj.r.squared)
	rec_sub$R2[1] <- max(t)
	rec_sub$DS[1] <- 0
	ra <- cor.test(mat[names(rk)[1],], rec_val_tt, method="spearman")
	rb <- cor.test(mat[names(rk)[2],], rec_val_tt, method="spearman")
	rec_sub$Cor[1] <- max(abs(ra$estimate), abs(rb$estimate))
	rec_sub$CP[1] <- min(ra$p.value, rb$p.value)
	return(rec_sub)
})
stopCluster(cl)
rec <- do.call("rbind", results)
rec$p.adj <- p.adjust(rec$Pval, "BH")
rec$cp.adj <- p.adjust(rec$CP, "BH")
write.csv(rec, "cmp_vtt_raw2.csv")

cmp_raw <- read.csv("cmp_vdv_raw2.csv", h=T, r=1)
cmp_raw <- cmp_raw[which(cmp_raw$Gene != "-" & cmp_raw$Num > 1 & cmp_raw$Count > 2 & cmp_raw$CM > 5),]
cmp_raw$TP[which(is.na(cmp_raw$TP))] <- 1
cmp_raw$Cor[which(is.na(cmp_raw$Cor))] <- 0
cmp_raw$p.adj <- cmp_raw$Pval
cmp_raw$cp.adj <- cmp_raw$CP
cmp_raw <- cmp_raw[order(cmp_raw$p.adj),]
cmp_raw$Symbol <- osn_sct_raw[["features"]]$Symbol[match(cmp_raw$Gene, osn_sct_raw[["features"]]$Gene)]
cmp_raw$Y <- -log10(cmp_raw$TP)
cmp_raw$Y[which(cmp_raw$Y > 10)] <- 10
cmp_raw$X <- cmp_raw$Cor
cmp_raw$Group <- "N.S."
cmp_raw$Group[which(cmp_raw$X > 0.3 & cmp_raw$Y > 1.3)] <- "IDG"
cmp_raw$Group <- factor(cmp_raw$Group, levels=c("IDG", "N.S."))
cmp_raw$Type <- "D > V"
cmp_raw$Val <- (cmp_raw$Y/max(cmp_raw$Y))^2 + (cmp_raw$X/max(cmp_raw$X))^2
cmp_raw <- cmp_raw[order(cmp_raw$Val, decreasing=T),]
cmp_raw$Anno <- ""
#cmp_raw$Anno[which(cmp_raw$Symbol == "Chd2" & cmp_raw$Pval < 0.05)] <- "Chd2"
#cmp_raw$Anno[which(cmp_raw$Symbol == "Psmd10" & cmp_raw$Pval < 0.05)] <- "Psmd10"
cmp_raw$ET <- "Diff"
for (i in 1:nrow(cmp_raw))
{
	if (cmp_raw$Group[i] != "IDG") next
	ids <- match(unlist(strsplit(cmp_raw$Trans[i], split=";")), osn_sct_raw[["features"]]$Name)
	strand <- osn_sct_raw[["features"]]$Strand[ids[1]]
	pos_list <- c()
	for (id in ids)
	{
		ts <- unlist(strsplit(osn_sct_raw[["features"]]$Exon[id], ";"))
		ts <- ts[min(2, length(ts))]
		if (strand == "+") pos_list <- c(pos_list, gsub("-.*", "", ts)) else pos_list <- c(pos_list, gsub(".*-", "", ts))
	}
	if (length(unique(pos_list)) == 1) cmp_raw$ET[i] <- "AS"
}
#write.csv(cmp_raw, "cmp_vdv_filter3.csv")
cmp_vdv_raw <- cmp_raw
cmp_vdv_raw <- read.csv("cmp_vdv_filter3.csv", h=T, r=1)
cmp_vdv_raw$Anno <- ""
aid <- which(cmp_vdv_raw$Group == "IDG")
cmp_vdv_raw$Anno[aid[1:min(10, length(aid))]] <- cmp_vdv_raw$Symbol[aid[1:min(10, length(aid))]]
cmp_vdv_raw$Annoc <- "black"
cmp_vdv_raw$Annoc[which(cmp_vdv_raw$Anno == "Fgf12")] <- "red"

cmp_raw <- read.csv("cmp_vap_raw2.csv", h=T, r=1)
cmp_raw <- cmp_raw[which(cmp_raw$Gene != "-" & cmp_raw$Num > 1 & cmp_raw$Count > 2 & cmp_raw$CM > 5),]
cmp_raw$TP[which(is.na(cmp_raw$TP))] <- 1
cmp_raw$Cor[which(is.na(cmp_raw$Cor))] <- 0
cmp_raw$p.adj <- cmp_raw$Pval
cmp_raw$cp.adj <- cmp_raw$CP
cmp_raw <- cmp_raw[order(cmp_raw$p.adj),]
cmp_raw$Symbol <- osn_sct_raw[["features"]]$Symbol[match(cmp_raw$Gene, osn_sct_raw[["features"]]$Gene)]
cmp_raw$Y <- -log10(cmp_raw$TP)
cmp_raw$Y[which(cmp_raw$Y > 10)] <- 10
cmp_raw$X <- cmp_raw$Cor
cmp_raw$Group <- "N.S."
cmp_raw$Group[which(cmp_raw$X > 0.3 & cmp_raw$Y > 1.3)] <- "IDG"
cmp_raw$Group <- factor(cmp_raw$Group, levels=c("IDG", "N.S."))
cmp_raw$Type <- "D > V"
cmp_raw$Val <- (cmp_raw$Y/max(cmp_raw$Y))^2 + (cmp_raw$X/max(cmp_raw$X))^2
cmp_raw <- cmp_raw[order(cmp_raw$Val, decreasing=T),]
cmp_raw$Anno <- ""
#cmp_raw$Anno[which(cmp_raw$Symbol == "Chd2" & cmp_raw$Pval < 0.05)] <- "Chd2"
#cmp_raw$Anno[which(cmp_raw$Symbol == "Psmd10" & cmp_raw$Pval < 0.05)] <- "Psmd10"
cmp_raw$ET <- "Diff"
for (i in 1:nrow(cmp_raw))
{
	if (cmp_raw$Group[i] != "IDG") next
	ids <- match(unlist(strsplit(cmp_raw$Trans[i], split=";")), osn_sct_raw[["features"]]$Name)
	strand <- osn_sct_raw[["features"]]$Strand[ids[1]]
	pos_list <- c()
	for (id in ids)
	{
		ts <- unlist(strsplit(osn_sct_raw[["features"]]$Exon[id], ";"))
		ts <- ts[min(2, length(ts))]
		if (strand == "+") pos_list <- c(pos_list, gsub("-.*", "", ts)) else pos_list <- c(pos_list, gsub(".*-", "", ts))
	}
	if (length(unique(pos_list)) == 1) cmp_raw$ET[i] <- "AS"
}
#write.csv(cmp_raw, "cmp_vap_filter3.csv")
cmp_vap_raw <- cmp_raw
cmp_vap_raw <- read.csv("cmp_vap_filter3.csv", h=T, r=1)
cmp_vap_raw$Anno <- ""
aid <- which(cmp_vap_raw$Group == "IDG")
cmp_vap_raw$Anno[aid[1:min(15, length(aid))]] <- cmp_vap_raw$Symbol[aid[1:min(15, length(aid))]]
cmp_vap_raw$Annoc <- "black"
cmp_vap_raw$Annoc[which(cmp_vap_raw$Anno == "Anapc16")] <- "red"

cmp_raw <- read.csv("cmp_vtt_raw2.csv", h=T, r=1)
cmp_raw <- cmp_raw[which(cmp_raw$Gene != "-" & cmp_raw$Num > 1 & cmp_raw$Count > 2 & cmp_raw$CM > 5),]
cmp_raw$TP[which(is.na(cmp_raw$TP))] <- 1
cmp_raw$Cor[which(is.na(cmp_raw$Cor))] <- 0
cmp_raw$p.adj <- cmp_raw$Pval
cmp_raw$cp.adj <- cmp_raw$CP
cmp_raw <- cmp_raw[order(cmp_raw$p.adj),]
cmp_raw$Symbol <- osn_sct_raw[["features"]]$Symbol[match(cmp_raw$Gene, osn_sct_raw[["features"]]$Gene)]
cmp_raw$Y <- -log10(cmp_raw$TP)
cmp_raw$Y[which(cmp_raw$Y > 10)] <- 10
cmp_raw$X <- cmp_raw$Cor
cmp_raw$Group <- "N.S."
cmp_raw$Group[which(cmp_raw$X > 0.3 & cmp_raw$Y > 1.3)] <- "IDG"
cmp_raw$Group <- factor(cmp_raw$Group, levels=c("IDG", "N.S."))
cmp_raw$Type <- "D > V"
cmp_raw$Val <- (cmp_raw$Y/max(cmp_raw$Y))^2 + (cmp_raw$X/max(cmp_raw$X))^2
cmp_raw <- cmp_raw[order(cmp_raw$Val, decreasing=T),]
cmp_raw$Anno <- ""
#cmp_raw$Anno[which(cmp_raw$Symbol == "Chd2" & cmp_raw$Pval < 0.05)] <- "Chd2"
#cmp_raw$Anno[which(cmp_raw$Symbol == "Psmd10" & cmp_raw$Pval < 0.05)] <- "Psmd10"
cmp_raw$ET <- "Diff"
for (i in 1:nrow(cmp_raw))
{
	if (cmp_raw$Group[i] != "IDG") next
	ids <- match(unlist(strsplit(cmp_raw$Trans[i], split=";")), osn_sct_raw[["features"]]$Name)
	strand <- osn_sct_raw[["features"]]$Strand[ids[1]]
	pos_list <- c()
	for (id in ids)
	{
		ts <- unlist(strsplit(osn_sct_raw[["features"]]$Exon[id], ";"))
		ts <- ts[min(2, length(ts))]
		if (strand == "+") pos_list <- c(pos_list, gsub("-.*", "", ts)) else pos_list <- c(pos_list, gsub(".*-", "", ts))
	}
	if (length(unique(pos_list)) == 1) cmp_raw$ET[i] <- "AS"
}
write.csv(cmp_raw, "cmp_vtt_filter2.csv")
cmp_vtt_raw <- cmp_raw
cmp_vtt_raw <- read.csv("cmp_vtt_filter2.csv", h=T, r=1)
cmp_vtt_raw$Anno <- ""
aid <- which(cmp_vtt_raw$Group == "IDG")
cmp_vtt_raw$Anno[aid[1:min(10, length(aid))]] <- cmp_vtt_raw$Symbol[aid[1:min(10, length(aid))]]

rec_c <- data.frame()
rec_p <- data.frame()
for (x in cmp_vdv_raw$Gene)
{
	gene <- osn_sct_raw[["features"]]$Symbol[match(x, osn_sct_raw[["features"]]$Gene)]
	ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
	mat <- as.matrix(osn_sct_dv[ids,])
	ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
	if (length(ts) < 2) next
	mat <- mat[ts,]
	bk <- ceiling(bins/3)-1
	#pp <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
	pp <- osn_sct_dv_lh[ids,]
	mat_count <- mat
	for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
	ppr <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
	rk <- abs(mat[, 1] - mat[, bins])
	rk <- rk[order(rk, decreasing=T)]
	rec_c <- rbind(rec_c, cbind(data.frame(Trans=names(rk)[1:2], Gene=x, Symbol=gene), mat_count[names(rk)[1:2],]))
	rec_p <- rbind(rec_p, cbind(data.frame(Trans=names(rk)[1:2], Gene=x, Symbol=gene), mat[names(rk)[1:2],]))
}
write.csv(rec_c, "cmp_vdv_filter2_rdc.csv", quote=F, row.names=F)
write.csv(rec_p, "cmp_vdv_filter2_rdp.csv", quote=F, row.names=F)
rec_c <- data.frame()
rec_p <- data.frame()
for (x in cmp_vap_raw$Gene)
{
	gene <- osn_sct_raw[["features"]]$Symbol[match(x, osn_sct_raw[["features"]]$Gene)]
	ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
	mat <- as.matrix(osn_sct_ap[ids,])
	ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
	if (length(ts) < 2) next
	mat <- mat[ts,]
	bk <- ceiling(bins/3)-1
	#pp <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
	pp <- osn_sct_ap_lh[ids,]
	mat_count <- mat
	for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
	ppr <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
	rk <- abs(mat[, 1] - mat[, bins])
	rk <- rk[order(rk, decreasing=T)]
	rec_c <- rbind(rec_c, cbind(data.frame(Trans=names(rk)[1:2], Gene=x, Symbol=gene), mat_count[names(rk)[1:2],]))
	rec_p <- rbind(rec_p, cbind(data.frame(Trans=names(rk)[1:2], Gene=x, Symbol=gene), mat[names(rk)[1:2],]))
}
write.csv(rec_c, "cmp_vap_filter2_rdc.csv", quote=F, row.names=F)
write.csv(rec_p, "cmp_vap_filter2_rdp.csv", quote=F, row.names=F)
rec_c <- data.frame()
rec_p <- data.frame()
for (x in cmp_vtt_raw$Gene)
{
	gene <- osn_sct_raw[["features"]]$Symbol[match(x, osn_sct_raw[["features"]]$Gene)]
	ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
	mat <- as.matrix(osn_sct_tt[ids,])
	ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
	if (length(ts) < 2) next
	mat <- mat[ts,]
	bk <- ceiling(bins/3)-1
	#pp <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
	pp <- osn_sct_tt_lh[ids,]
	mat_count <- mat
	for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
	ppr <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
	rk <- abs(mat[, 1] - mat[, bins])
	rk <- rk[order(rk, decreasing=T)]
	rec_c <- rbind(rec_c, cbind(data.frame(Trans=names(rk)[1:2], Gene=x, Symbol=gene), mat_count[names(rk)[1:2],]))
	rec_p <- rbind(rec_p, cbind(data.frame(Trans=names(rk)[1:2], Gene=x, Symbol=gene), mat[names(rk)[1:2],]))
}
write.csv(rec_c, "cmp_vtt_filter2_rdc.csv", quote=F, row.names=F)
write.csv(rec_p, "cmp_vtt_filter2_rdp.csv", quote=F, row.names=F)

#x_max <- max(max(max(cmp_vdv_raw$X), cmp_vap_raw$X), cmp_vtt_raw$X)
#y_max <- max(max(max(cmp_vdv_raw$Y), cmp_vap_raw$Y), cmp_vtt_raw$Y)
x_max <- max(max(cmp_vdv_raw$X), cmp_vap_raw$X)
y_max <- max(max(cmp_vdv_raw$Y), cmp_vap_raw$Y)
pea <- ggplot(cmp_vdv_raw, aes(x=X, y=Y, fill=Group, color=Group))+geom_point(shape=21, size=2)+
	labs(title=paste("D > V (#IDG: ", length(which(cmp_vdv_raw$Group == "IDG")), ")"), 
	x="Correlation", y="Pval (-log10)", color=NULL, fill=NULL)+
	#geom_text_repel(data=cmp_vdv_raw[which(cmp_vdv_raw$Anno != "" & cmp_vdv_raw$Annoc == "black"),,drop=F], aes(label=Anno), color="black", size=2, 
	#segment.size=0.5, direction="y", nudge_y=0.2, nudge_x=0, hjust=0)+
	geom_text_repel(data=cmp_vdv_raw[which(cmp_vdv_raw$Anno != "" & cmp_vdv_raw$Annoc == "red"),,drop=F], aes(label=Anno), 
	color="black", size=4, segment.size=0.5, direction="y", nudge_y=0.2, nudge_x=0, hjust=0)+
	scale_x_continuous(limits=c(0, x_max+0.05), expand=c(0, 0))+
	scale_y_continuous(limits=c(0, y_max+1), expand=c(0, 0))+
	scale_fill_manual(values=c("red4", "grey60"))+
	scale_color_manual(values=c("black", "grey60"))+
	guides(color=guide_legend(override.aes=list(size=4)))+
	geom_vline(xintercept=0.3, color="gray30", linetype="dashed", linewidth=0.6)+
	geom_hline(yintercept=1.3, color="gray30", linetype="dashed", linewidth=0.6)+
	#geom_text(aes(x=X, y=Y+0.1, label=Anno), size=0.3*text_size, vjust=0, color="black")+
	theme(axis.line=element_line(linetype=1, color='black'), panel.background=element_rect(0, linetype=0), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black", face="plain"), 
	axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	legend.box.background=element_rect(linewidth=0.4, color="black"),
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), legend.key.size=unit(15, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), legend.position=c(0.85, 0.83))+tag_thm

peb <- ggplot(cmp_vap_raw, aes(x=X, y=Y, fill=Group, color=Group))+geom_point(shape=21, size=2)+
	labs(title=paste("A > P (#IDG: ", length(which(cmp_vap_raw$Group == "IDG")), ")"), 
	x="Correlation", y="Pval (-log10)", color=NULL, fill=NULL)+
	#geom_text_repel(data=cmp_vap_raw[which(cmp_vap_raw$Anno != "" & cmp_vap_raw$Annoc == "black"),,drop=F], aes(label=Anno), color="black", size=2, 
	#segment.size=0.5, direction="y", nudge_y=0.2, nudge_x=0, hjust=0)+
	geom_text_repel(data=cmp_vap_raw[which(cmp_vap_raw$Anno != "" & cmp_vap_raw$Annoc == "red"),,drop=F], aes(label=Anno), 
	color="black", size=4, segment.size=0.5, direction="y", nudge_y=0.2, nudge_x=0, hjust=0)+
	scale_x_continuous(limits=c(0, x_max+0.05), expand=c(0, 0))+
	scale_y_continuous(limits=c(0, y_max+1), expand=c(0, 0))+
	scale_fill_manual(values=c("red4", "grey60"))+
	scale_color_manual(values=c("black", "grey60"))+
	guides(color=guide_legend(override.aes=list(size=4)))+
	geom_vline(xintercept=0.3, color="gray30", linetype="dashed", linewidth=0.6)+
	geom_hline(yintercept=1.3, color="gray30", linetype="dashed", linewidth=0.6)+
	#geom_text(aes(x=X, y=Y+0.1, label=Anno), size=0.3*text_size, vjust=0, color="black")+
	theme(axis.line=element_line(linetype=1, color='black'), panel.background=element_rect(0, linetype=0), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black", face="plain"), 
	axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.box.background=element_rect(linewidth=0.4, color="black"),
	legend.text=element_text(size=text_size, color="black"), legend.key.size=unit(15, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), legend.position=c(0.85, 0.83))+tag_thm

x <- "Nqo1"
id <- osn_sct_raw[["features"]]$Gene[match(x, osn_sct_raw[["features"]]$Symbol)]
ids <- which(osn_sct_raw[["features"]]$Gene == id & osn_sct_raw[["features"]]$Type == "protein_coding")
rec <- data.frame(Count=as.numeric(colSums(osn_sct_dv[ids,])), Score=rec_val_dv)
cr <- cor.test(rec$Count, rec$Score)
pca <- ggplot(rec, aes(x=Score, y=Count))+geom_point(color="grey40", size=3)+
	labs(title=x, x="Score of D > V", y="Counts")+
	annotate("text", x=0.18, y=max(rec$Count), label=paste0("R=", round(cr$estimate, 2), "\nPval=", round(cr$p.value, 2)), 
	color="black", size=5, hjust=0, vjust=1)+
	stat_smooth(method=lm, se=F, color=col_list[13], linetype="dashed", linewidth=1)+
	#scale_x_continuous(limits=c(min(rec$val)-0.01, max(rec$val)+0.01), expand=c(0, 0))+
	#scale_y_continuous(limits=c(min(rec$expr)-0.01, max(rec$expr)+0.01), expand=c(0, 0))+
	theme(axis.line=element_line(linetype=1, color='black'), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin())+tag_thm
x <- "Nfix"
id <- osn_sct_raw[["features"]]$Gene[match(x, osn_sct_raw[["features"]]$Symbol)]
ids <- which(osn_sct_raw[["features"]]$Gene == id & osn_sct_raw[["features"]]$Type == "protein_coding")
rec <- data.frame(Count=as.numeric(colSums(osn_sct_dv[ids,])), Score=rec_val_dv)
cr <- cor.test(rec$Count, rec$Score)
pcb <- ggplot(rec, aes(x=Score, y=Count))+geom_point(color="grey40", size=3)+
	labs(title=x, x="Score of D > V", y="Counts")+
	annotate("text", x=min(rec$Score), y=max(rec$Count), label=paste0("R=", round(cr$estimate, 2), "\nPval=", round(cr$p.value, 2)), 
	color="black", size=5, hjust=0, vjust=1)+
	stat_smooth(method=lm, se=F, color=col_list[13], linetype="dashed", linewidth=1)+
	#scale_x_continuous(limits=c(min(rec$val)-0.01, max(rec$val)+0.01), expand=c(0, 0))+
	#scale_y_continuous(limits=c(min(rec$expr)-0.01, max(rec$expr)+0.01), expand=c(0, 0))+
	theme(axis.line=element_line(linetype=1, color='black'), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin())+tag_thm
x <- "Tmsb4x"
id <- osn_sct_raw[["features"]]$Gene[match(x, osn_sct_raw[["features"]]$Symbol)]
ids <- which(osn_sct_raw[["features"]]$Gene == id & osn_sct_raw[["features"]]$Type == "protein_coding")
rec <- data.frame(Count=as.numeric(colSums(osn_sct_ap[ids,])), Score=rec_val_ap)
cr <- cor.test(rec$Count, rec$Score)
pda <- ggplot(rec, aes(x=Score, y=Count))+geom_point(color="grey40", size=3)+
	labs(title=x, x="Score of A > P", y="Counts")+
	annotate("text", x=0.2, y=max(rec$Count), label=paste0("R=", round(cr$estimate, 2), "\nPval=", round(cr$p.value, 2)), 
	color="black", size=5, hjust=0, vjust=1)+
	stat_smooth(method=lm, se=F, color=col_list[13], linetype="dashed", linewidth=1)+
	#scale_x_continuous(limits=c(min(rec$val)-0.01, max(rec$val)+0.01), expand=c(0, 0))+
	#scale_y_continuous(limits=c(min(rec$expr)-0.01, max(rec$expr)+0.01), expand=c(0, 0))+
	theme(axis.line=element_line(linetype=1, color='black'), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin())+tag_thm
x <- "Nrp1"
id <- osn_sct_raw[["features"]]$Gene[match(x, osn_sct_raw[["features"]]$Symbol)]
ids <- which(osn_sct_raw[["features"]]$Gene == id & osn_sct_raw[["features"]]$Type == "protein_coding")
rec <- data.frame(Count=as.numeric(colSums(osn_sct_ap[ids,])), Score=rec_val_ap)
cr <- cor.test(rec$Count, rec$Score)
pdb <- ggplot(rec, aes(x=Score, y=Count))+geom_point(color="grey40", size=3)+
	labs(title=x, x="Score of A > P", y="Counts")+
	annotate("text", x=min(rec$Score), y=max(rec$Count), label=paste0("R=", round(cr$estimate, 2), "\nPval=", round(cr$p.value, 2)), 
	color="black", size=5, hjust=0, vjust=1)+
	stat_smooth(method=lm, se=F, color=col_list[13], linetype="dashed", linewidth=1)+
	#scale_x_continuous(limits=c(min(rec$val)-0.01, max(rec$val)+0.01), expand=c(0, 0))+
	#scale_y_continuous(limits=c(min(rec$expr)-0.01, max(rec$expr)+0.01), expand=c(0, 0))+
	theme(axis.line=element_line(linetype=1, color='black'), 
	plot.title=element_text(size=title_size, hjust=0.5, color="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin())+tag_thm
#ggsave(plot=wrap_plots(list(pca, pcb, pda, pdb), nrow=2), width=10, height=8, dpi=200, filename="test.png", limitsize=F)

gene <- "Fgf12"
junction_info <- read.csv("/mnt/md0/oe_full_length/output/OSNfulllength/res_info_junction_fgf12_raw.csv", h=T)
density_info <- read.csv("/mnt/md0/oe_full_length/output/OSNfulllength/res_info_density_fgf12_raw.csv", h=T)
x <- osn_sct_raw[["features"]]$Gene[match(gene, osn_sct_raw[["features"]]$Symbol)]
ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
mat <- as.matrix(osn_sct_dv[ids,])
ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
mat <- mat[ts,]
bk <- ceiling(bins/3)-1
#pp <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
pp <- osn_sct_dv_lh[ids,]
mat_count <- mat
for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
ppr <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
rk <- abs(mat[, 1] - mat[, bins])
rk <- rk[order(rk, decreasing=T)]
rec <- data.frame(RNA=as.numeric(colSums(osn_sct_dv[ids,])), Count=mat_count[names(rk)[1],], Percent=mat[names(rk)[1],], Score=rec_val_dv)
cr <- cor.test(rec$Percent, rec$Score)
pib <- wrap_elements(ggplot(rec, aes(x=Score, y=Percent))+geom_point(color=col_list[6], size=3)+
	labs(title=names(rk)[1], x="Score of D > V", y="Percentage (%)")+
	annotate("text", x=min(rec$Score), y=max(rec$Percent), label=paste0("R=", round(cr$estimate, 2), "\nPval=", round(cr$p.value, 2)), 
	color="black", size=5, hjust=0, vjust=1)+
	stat_smooth(method=lm, se=F, color=col_list[13], linetype="dashed", linewidth=1)+
	scale_x_continuous(breaks=seq(-0.4, 0.4, 0.2))+
	#scale_y_continuous(limits=c(min(rec$expr)-0.01, max(rec$expr)+0.01), expand=c(0, 0))+
	theme(axis.line=element_line(linetype=1, color='black'), 
	plot.title=element_text(size=15, hjust=0.5, color=col_list[6]), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin()))+tag_thm
cell_info <- data.frame(cell=colnames(osn_rna), val=osn_rna$vdv)
cell_info <- cell_info[order(cell_info$val),]
cell_info$Group <- factor(rep(seq(1, ncol(osn_rna), round(ncol(osn_rna)/5))[1:5], each=round(ncol(osn_rna)/5))[1:ncol(osn_rna)], labels=1:5)
groups <- 1:5
osn_rna$Group <- cell_info$Group[match(colnames(osn_rna), rownames(cell_info))]
cell_num <- min(table(osn_rna$Group))
col_groups <- colorRampPalette(brewer.pal(9,"YlOrRd")[c(1,5,9)])(length(groups))
names(col_groups) <- groups
rec_ppb <- data.frame(Time=seq(min(osn_rna$vdv), max(osn_rna$vdv), 0.01))
piac <- ggplot(rec_ppb, aes(x=Time, y="", color=Time))+geom_point(size=2, alpha=0.5, shape=15)+
	labs(title=NULL, x=NULL, y=NULL)+scale_color_viridis()+
	scale_x_continuous(breaks=c(-0.8, 0.8), labels=c("D", "V"), expand=c(0,0))+
	scale_y_discrete(expand=c(0,0))+
	guides(color=guide_legend(override.aes=list(size=4, alpha=1)))+
	theme(panel.background=element_blank(), axis.line=element_blank(), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), legend.key.size=unit(15, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text.y=element_blank(), axis.text.x=element_text(size=text_size, color="black"), 
	axis.ticks=element_blank(), legend.position="none")+tag_thm
ids <- which(osn_sct_raw[["features"]]$Symbol == gene & osn_sct_raw[["features"]]$Type == "protein_coding")
rec_tt <- matrix(0, nrow=length(ids), ncol=length(groups), dimnames=list(osn_sct_raw[["features"]]$ID[ids], groups))
for (i in 1:length(groups))
{
	cells <- which(osn_rna$Group == i)
	rec_tt[, i] <- rowSums(osn_sct_raw[["trans"]][ids, cells])*cell_num/length(cells)
}
ids_sub <- which(rowMaxs(rec_tt) > 4)
if (length(ids_sub) < 2) next
rec_tt <- rec_tt[ids_sub,]
for (i in 1:ncol(rec_tt)) rec_tt[, i] <- rec_tt[, i]*100/sum(rec_tt[, i])
rec_ttt <- data.frame()
for (i in which(colSums(rec_tt) > 0)) rec_ttt <- rbind(rec_ttt, data.frame(Type=i, 
	Isoform=rownames(rec_tt), Count=rec_tt[, i]))
rec_ttt$Type <- factor(rec_ttt$Type, levels=1:length(groups))
rec_ttt$Isoform <- factor(rec_ttt$Isoform, levels=rev(rownames(rec_tt)))
piaa <- ggplot(rec_ttt, aes(x=Type, y=Count, stratum=Isoform, alluvium=Isoform, fill=Isoform))+
	labs(title=NULL, x=NULL, y="Percentage (%)")+
	geom_flow(alpha=1)+geom_stratum(alpha=1, color="white", linewidth=0)+
	scale_y_continuous(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
	scale_fill_manual(values=col_list[c(1:(nrow(rec_tt)))+5])+
	theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,color="black"), 
	axis.title=element_text(size=title_size, color="black"), axis.text.y=element_text(size=text_size, color="black"), 
	axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="none")+tag_thm
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
info_ranges <- paste0(min(info_exons$start), "-", max(info_exons$start), "\n", info_exons$seqnames[1], "(", strand, ")")
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
info_introns$transcript_name <- factor(info_introns$transcript_name, levels=rownames(rec_tt))
info_exons$transcript_name <- factor(info_exons$transcript_name, levels=rownames(rec_tt))
info_utrs$transcript_name <- factor(info_utrs$transcript_name, levels=rownames(rec_tt))
info_introns$rs <- coord_trans$r[match(info_introns$start, coord_trans$i)]
info_introns$re <- coord_trans$r[match(info_introns$end, coord_trans$i)]
info_utrs$rs <- coord_trans$r[match(info_utrs$start, coord_trans$i)]
info_utrs$re <- coord_trans$r[match(info_utrs$end, coord_trans$i)]
info_exons$rs <- coord_trans$r[match(info_exons$start, coord_trans$i)]
info_exons$re <- coord_trans$r[match(info_exons$end, coord_trans$i)]
info_col <- rev(col_list[c(1:(nrow(rec_tt)))+5])
names(info_col) <- rownames(rec_tt)
piab <- ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
	geom_intron(data=info_introns, aes(color=transcript_name, strand=strand), linewidth=0.5, arrow.min.intron.length=300)+
	geom_range(aes(fill=transcript_name, color=transcript_name), linewidth=0, height=0.3)+
	geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name, fill=transcript_name, 
	color=transcript_name), linewidth=0.5, height=0.1)+
	labs(title=NULL, x=NULL, y=NULL)+
	scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, labels=info_ranges)+
	scale_y_discrete(expand=expansion(mult=c(0.1,0.23)))+
	geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), 
	x=min(info_utrs$start, info_exons$start), hjust=0, vjust=0, nudge_y=0.3, size=4)+
	scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
	axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
	axis.text.x=element_text(size=text_size, color="black"), legend.position="none")+tag_thm
pia <- wrap_elements(wrap_plots(A=piaa, B=piab, C=piac, design=c("AB\nCB"), heights=c(40, 1))+
	plot_annotation(gene, theme=theme(plot.title=element_text(size=title_size, colour="black", hjust=0.5), 
	plot.margin=margin(0,10,0,10), panel.spacing=unit(0, "pt"))))+tag_thm
junction_info$is <- coord_trans$i[match(junction_info$xa-1, coord_trans$r)]
junction_info$ie <- coord_trans$i[match(junction_info$xb, coord_trans$r)]
density_info <- density_info[order(density_info$x),]
density_info$x <- density_info$x + 1
density_info <- rbind(density_info, density_info[nrow(density_info),])
density_info$x[nrow(density_info)] <- density_info$x[nrow(density_info)] + 1
density_info$type <- "intron"
for (i in 1:nrow(info_utrs)) density_info$type[which(density_info$x >= info_utrs$rs[i] & density_info$x <= info_utrs$re[i])] <- "exon"
for (i in 1:nrow(info_exons)) density_info$type[which(density_info$x >= info_exons$rs[i] & density_info$x <= info_exons$re[i])] <- "exon"
t <- 1
density_regions <- data.frame()
for (i in 2:nrow(density_info)) if (density_info$type[i] != density_info$type[t] | i == nrow(density_info))
{
	if (density_info$type[t] == "exon") density_regions <- rbind(density_regions, 
		data.frame(start=density_info$x[t], end=density_info$x[i] - 1, type=density_info$type[t]))
	else if (density_info$type[t] == "intron") density_regions <- rbind(density_regions, 
		data.frame(start=density_info$x[t] - 1, end=density_info$x[i], type=density_info$type[t]))
	t <- i
}
density_regions$is <- coord_trans$i[match(density_regions$start, coord_trans$r)]
density_regions$ie <- coord_trans$i[match(density_regions$end, coord_trans$r)]
density_regions$ie[nrow(density_regions)] <- density_regions$ie[nrow(density_regions)] + 1
if (density_regions$type[nrow(density_regions)] == "intron") density_regions <- density_regions[1:nrow(density_regions)-1,]
density_regions$il <- density_regions$ie - density_regions$is
density_regions$sc <- (density_regions$end - density_regions$start) / density_regions$il
density_info <- density_info[1:(nrow(density_info)-1),]
gene_list <- gene
junc_list <- osn_sct_raw[["features"]]$Gene[match(gene, osn_sct_raw[["features"]]$Symbol)]
sample_list <- unique(density_info$s)
color_list <- setNames(rev(col_list[c(1, 2)]), sample_list)
pls <- lapply(sample_list, function(id)
{
	density <- density_info[which(density_info$g == gene & density_info$s == id),]
	density_fix <- data.frame()
	for (i in 1:nrow(density_regions))
	{
		for (j in 1:density_regions$il[i]) density_fix <- rbind(density_fix, data.frame(x=j+density_regions$is[i]-1, 
			y=mean(density$y[which(density$x >= round((j-1)*density_regions$sc[i])+density_regions$start[i] & 
			density$x <= round(j*density_regions$sc[i])+density_regions$start[i])])))
	}
	density <- density_fix
	junction = junction_info[which(junction_info$g == gene & junction_info$s == id),]
	junction = junction[order(junction$xa, junction$xb),]
	junction$xa <- junction$is
	junction$xb <- junction$ie
	junction$term = paste0(junction$g, ":", junction$ao, "-", junction$bo)
	junction$xmid <- 0
	gp = ggplot(density)+geom_bar(aes(x, y), width=1, position='identity', stat='identity', fill=color_list[[id]])+
		scale_x_continuous(expand=c(0.01, 0.01))+labs(y=id)+
		theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,color="black"), 
		axis.title.y=element_text(size=title_size, color="black", angle=0, vjust=0.5), 
		axis.title.x=element_blank(), axis.line.x=element_blank(), 
		axis.text=element_blank(), axis.ticks=element_blank(), plot.margin=margin(b=4))
	ymax = max(density$y)*1.1
	ymin = -max(density$y)*0.1
	for (i in 1:nrow(junction))
	{
		j = as.numeric(junction[i,1:5])
		xmid = mean(j[1:2])
		junction$xmid[i] <- xmid
		curve_par = gpar(lwd=1.5, col=color_list[[id]])
		#if (i > (nrow(junction) - 2)) curve_par = gpar(lwd=1.5, col="red")
		pcol <- "black"
		if (i > (nrow(junction) - 2)) pcol <- "red"
		if (i%%2 == 0) {
			ymid = -runif(1, 0.1, 0.3)*max(density$y)
			ymin = min(ymin, ymid*1.1)
			gp = gp+annotation_custom(xsplineGrob(c(0, 0, 1, 1), c(1, 0, 0, 0), shape=1, gp=curve_par), j[1], xmid, 0, ymid)+
				annotation_custom(xsplineGrob(c(1, 1, 0, 0), c(1, 0, 0, 0), shape=1, gp=curve_par), xmid, j[2], 0, ymid)
		} else {
			ymid = runif(1, 1.2, 1.4)*max(j[3:4])
			ymax = max(ymax, ymid*1.1)
			gp = gp+annotation_custom(xsplineGrob(c(0, 0, 1, 1), c(0, 1, 1, 1), shape=1, gp=curve_par), j[1], xmid, j[3], ymid)+
				annotation_custom(xsplineGrob(c(1, 1, 0, 0), c(0, 1, 1, 1), shape=1, gp=curve_par), xmid, j[2], j[4], ymid)
		}
		gp = gp+annotate("label", x=xmid, y=ymid, label=as.character(j[5]), color=pcol, 
			vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"), size=5, border.color="white")
	}
	gp <- gp+scale_y_continuous(breaks=c(0, max(density$y)), limits=c(ymin, ymax))+coord_cartesian(clip="off")
	#ggsave(plot=gp, width=10, height=4, dpi=200, "test.png", limitsize=F)
	return(gp)
})
picb <- ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
	geom_intron(data=info_introns, aes(color=transcript_name, strand=strand), linewidth=0.5, arrow.min.intron.length=300)+
	geom_range(aes(fill=transcript_name, color=transcript_name), linewidth=0, height=0.3)+
	geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name, fill=transcript_name, 
	color=transcript_name), linewidth=0.5, height=0.1)+
	labs(title=NULL, x=NULL, y=NULL)+
	scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, labels=info_ranges, expand=c(0.01, 0.01))+
	#geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), 
	#x=x_min, hjust=0, vjust=0, nudge_y=0.3, size=0.3*text_size)+
	scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
	axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
	axis.text.x=element_text(size=text_size, color="black"), legend.position="none")+tag_thm
pic <- wrap_elements(wrap_plots(c(pls, list(picb)), ncol=1, heights=rep(c(3, 3, 2), length(gene_list)))+
	plot_annotation(gene, theme=theme(plot.title=element_text(size=title_size, colour="black", hjust=0.5), 
	plot.margin=margin(0,0,0,0), panel.spacing=unit(0, "pt"))))+tag_thm

gene <- "Anapc16"
junction_info <- read.csv("/mnt/md0/oe_full_length/output/OSNfulllength/res_info_junction_anapc16_raw_ap.csv", h=T)
density_info <- read.csv("/mnt/md0/oe_full_length/output/OSNfulllength/res_info_density_anapc16_raw_ap.csv", h=T)
x <- osn_sct_raw[["features"]]$Gene[match(gene, osn_sct_raw[["features"]]$Symbol)]
ids <- which(osn_sct_raw[["features"]]$Gene == x & osn_sct_raw[["features"]]$Type == "protein_coding")
mat <- as.matrix(osn_sct_ap[ids,])
ts <- which(rowMaxs(mat) > 2 & rowSums(mat) > bins)
mat <- mat[ts,]
bk <- ceiling(bins/3)-1
#pp <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
pp <- osn_sct_ap_lh[ids,]
mat_count <- mat
for (i in 1:ncol(mat)) mat[, i] <- mat[, i]*100/sum(mat[, i])
ppr <- cbind(rowSums(mat[, 1:bk]), rowSums(mat[, (bins-bk+1):bins]))
rk <- abs(mat[, 1] - mat[, bins])
rk <- rk[order(rk, decreasing=T)]
rec <- data.frame(RNA=as.numeric(colSums(osn_sct_ap[ids,])), Count=mat_count[names(rk)[1],], Percent=mat[names(rk)[1],], Score=rec_val_ap)
cr <- cor.test(rec$Percent, rec$Score)
pjb <- wrap_elements(ggplot(rec, aes(x=Score, y=Percent))+geom_point(color=col_list[7], size=3)+
	labs(title=names(rk)[1], x="Score of A > P", y="Percentage (%)")+
	annotate("text", x=min(rec$Score), y=min(rec$Percent), label=paste0("R=", round(cr$estimate, 2), "\nPval=", round(cr$p.value, 2)), 
	color="black", size=5, hjust=0, vjust=0)+
	stat_smooth(method=lm, se=F, color=col_list[13], linetype="dashed", linewidth=1)+
	#scale_x_continuous(limits=c(min(rec$val)-0.01, max(rec$val)+0.01), expand=c(0, 0))+
	#scale_y_continuous(limits=c(min(rec$expr)-0.01, max(rec$expr)+0.01), expand=c(0, 0))+
	theme(axis.line=element_line(linetype=1, color='black'), 
	plot.title=element_text(size=15, hjust=0.5, color=col_list[7]), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	panel.background=element_rect(0, linetype=0), plot.margin=margin()))+tag_thm
cell_info <- data.frame(cell=colnames(osn_rna), val=osn_rna$vap)
cell_info <- cell_info[order(cell_info$val),]
cell_info$Group <- factor(rep(seq(1, ncol(osn_rna), round(ncol(osn_rna)/5))[1:5], each=round(ncol(osn_rna)/5))[1:ncol(osn_rna)], labels=1:5)
groups <- 1:5
osn_rna$Group <- cell_info$Group[match(colnames(osn_rna), rownames(cell_info))]
cell_num <- min(table(osn_rna$Group))
col_groups <- colorRampPalette(brewer.pal(9,"YlOrRd")[c(1,5,9)])(length(groups))
names(col_groups) <- groups
rec_ppa <- data.frame(Time=seq(min(osn_rna$vap), max(osn_rna$vap), 0.01))
pjac <- ggplot(rec_ppa, aes(x=Time, y="", color=Time))+geom_point(size=2, alpha=0.5, shape=15)+
	labs(title=NULL, x=NULL, y=NULL)+scale_color_viridis()+
	scale_x_continuous(breaks=c(-0.8, 0.8), labels=c("A", "P"), expand=c(0,0))+
	scale_y_discrete(expand=c(0,0))+
	guides(color=guide_legend(override.aes=list(size=4, alpha=1)))+
	theme(panel.background=element_blank(), axis.line=element_blank(), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), legend.key.size=unit(15, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text.y=element_blank(), axis.text.x=element_text(size=text_size, color="black"), 
	axis.ticks=element_blank(), legend.position="none")+tag_thm
ids <- which(osn_sct_raw[["features"]]$Symbol == gene & osn_sct_raw[["features"]]$Type == "protein_coding")
rec_tt <- matrix(0, nrow=length(ids), ncol=length(groups), dimnames=list(osn_sct_raw[["features"]]$ID[ids], groups))
for (i in 1:length(groups))
{
	cells <- which(osn_rna$Group == i)
	rec_tt[, i] <- rowSums(osn_sct_raw[["trans"]][ids, cells])*cell_num/length(cells)
}
ids_sub <- which(rowMaxs(rec_tt) > 4)
if (length(ids_sub) < 2) next
rec_tt <- rec_tt[ids_sub,]
for (i in 1:ncol(rec_tt)) rec_tt[, i] <- rec_tt[, i]*100/sum(rec_tt[, i])
rec_ttt <- data.frame()
for (i in which(colSums(rec_tt) > 0)) rec_ttt <- rbind(rec_ttt, data.frame(Type=i, 
	Isoform=rownames(rec_tt), Count=rec_tt[, i]))
rec_ttt$Type <- factor(rec_ttt$Type, levels=1:length(groups))
rec_ttt$Isoform <- factor(rec_ttt$Isoform, levels=rev(rownames(rec_tt)))
pjaa <- ggplot(rec_ttt, aes(x=Type, y=Count, stratum=Isoform, alluvium=Isoform, fill=Isoform))+
	labs(title=NULL, x=NULL, y="Percentage (%)")+
	geom_flow(alpha=1)+geom_stratum(alpha=1, color="white", linewidth=0)+
	scale_y_continuous(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
	scale_fill_manual(values=col_list[c(1:(nrow(rec_tt)))+5])+
	theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,color="black"), 
	axis.title=element_text(size=title_size, color="black"), axis.text.y=element_text(size=text_size, color="black"), 
	axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="none")+tag_thm
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
info_ranges <- paste0(min(info_exons$start), "-", max(info_exons$start), "\n", info_exons$seqnames[1], "(", strand, ")")
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
info_introns$transcript_name <- factor(info_introns$transcript_name, levels=rownames(rec_tt))
info_exons$transcript_name <- factor(info_exons$transcript_name, levels=rownames(rec_tt))
info_utrs$transcript_name <- factor(info_utrs$transcript_name, levels=rownames(rec_tt))
info_introns$rs <- coord_trans$r[match(info_introns$start, coord_trans$i)]
info_introns$re <- coord_trans$r[match(info_introns$end, coord_trans$i)]
info_utrs$rs <- coord_trans$r[match(info_utrs$start, coord_trans$i)]
info_utrs$re <- coord_trans$r[match(info_utrs$end, coord_trans$i)]
info_exons$rs <- coord_trans$r[match(info_exons$start, coord_trans$i)]
info_exons$re <- coord_trans$r[match(info_exons$end, coord_trans$i)]
info_col <- rev(col_list[c(1:(nrow(rec_tt)))+5])
names(info_col) <- rownames(rec_tt)
pjab <- ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
	geom_intron(data=info_introns, aes(color=transcript_name, strand=strand), linewidth=0.5, arrow.min.intron.length=300)+
	geom_range(aes(fill=transcript_name, color=transcript_name), linewidth=0, height=0.3)+
	geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name, fill=transcript_name, 
	color=transcript_name), linewidth=0.5, height=0.1)+
	labs(title=NULL, x=NULL, y=NULL)+
	scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, labels=info_ranges)+
	scale_y_discrete(expand=expansion(mult=c(0.1,0.23)))+
	geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), 
	x=min(info_utrs$start, info_exons$start), hjust=0, vjust=0, nudge_y=0.3, size=4)+
	scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
	axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
	axis.text.x=element_text(size=text_size, color="black"), legend.position="none")+tag_thm
pja <- wrap_elements(wrap_plots(A=pjaa, B=pjab, C=pjac, design=c("AB\nCB"), heights=c(40, 1))+
	plot_annotation(gene, theme=theme(plot.title=element_text(size=title_size, colour="black", hjust=0.5), 
	plot.margin=margin(0,10,0,10), panel.spacing=unit(5, "pt"))))+tag_thm
junction_info$is <- coord_trans$i[match(junction_info$xa-1, coord_trans$r)]
junction_info$ie <- coord_trans$i[match(junction_info$xb, coord_trans$r)]
density_info <- density_info[order(density_info$x),]
density_info$x <- density_info$x + 1
density_info <- rbind(density_info, density_info[nrow(density_info),])
density_info$x[nrow(density_info)] <- density_info$x[nrow(density_info)] + 1
density_info$type <- "intron"
for (i in 1:nrow(info_utrs)) density_info$type[which(density_info$x >= info_utrs$rs[i] & density_info$x <= info_utrs$re[i])] <- "exon"
for (i in 1:nrow(info_exons)) density_info$type[which(density_info$x >= info_exons$rs[i] & density_info$x <= info_exons$re[i])] <- "exon"
t <- 1
density_regions <- data.frame()
for (i in 2:nrow(density_info)) if (density_info$type[i] != density_info$type[t] | i == nrow(density_info))
{
	if (density_info$type[t] == "exon") density_regions <- rbind(density_regions, 
		data.frame(start=density_info$x[t], end=density_info$x[i] - 1, type=density_info$type[t]))
	else if (density_info$type[t] == "intron") density_regions <- rbind(density_regions, 
		data.frame(start=density_info$x[t] - 1, end=density_info$x[i], type=density_info$type[t]))
	t <- i
}
density_regions$is <- coord_trans$i[match(density_regions$start, coord_trans$r)]
density_regions$ie <- coord_trans$i[match(density_regions$end, coord_trans$r)]
density_regions$ie[nrow(density_regions)] <- density_regions$ie[nrow(density_regions)] + 1
if (density_regions$type[nrow(density_regions)] == "intron") density_regions <- density_regions[1:nrow(density_regions)-1,]
density_regions$il <- density_regions$ie - density_regions$is
density_regions$sc <- (density_regions$end - density_regions$start) / density_regions$il
density_info <- density_info[1:(nrow(density_info)-1),]
gene_list <- gene
junc_list <- osn_sct_raw[["features"]]$Gene[match(gene, osn_sct_raw[["features"]]$Symbol)]
sample_list <- unique(density_info$s)
color_list <- setNames(rev(col_list[c(1, 2)]), sample_list)
pls <- lapply(sample_list, function(id)
{
	density <- density_info[which(density_info$g == gene & density_info$s == id),]
	density_fix <- data.frame()
	for (i in 1:nrow(density_regions))
	{
		for (j in 1:density_regions$il[i]) density_fix <- rbind(density_fix, data.frame(x=j+density_regions$is[i]-1, 
			y=mean(density$y[which(density$x >= round((j-1)*density_regions$sc[i])+density_regions$start[i] & 
			density$x <= round(j*density_regions$sc[i])+density_regions$start[i])])))
	}
	density <- density_fix
	junction = junction_info[which(junction_info$g == gene & junction_info$s == id),]
	junction = junction[order(junction$xa, junction$xb),]
	junction$xa <- junction$is
	junction$xb <- junction$ie
	junction$term = paste0(junction$g, ":", junction$ao, "-", junction$bo)
	junction$xmid <- 0
	gp = ggplot(density)+geom_bar(aes(x, y), width=1, position='identity', stat='identity', fill=color_list[[id]])+
		scale_x_continuous(expand=c(0.01, 0.01))+labs(y=id)+
		theme(panel.background=element_blank(), axis.line.y=element_line(linetype=1,color="black"), 
		axis.title.y=element_text(size=title_size, color="black", angle=0, vjust=0.5), 
		axis.title.x=element_blank(), axis.line.x=element_blank(), 
		axis.text=element_blank(), axis.ticks=element_blank(), plot.margin=margin(b=4))
	ymax = max(density$y)*1.1
	ymin = -max(density$y)*0.1
	for (i in 1:nrow(junction))
	{
		j = as.numeric(junction[i,1:5])
		xmid = mean(j[1:2])
		junction$xmid[i] <- xmid
		curve_par = gpar(lwd=1.5, col=color_list[[id]])
		#if (i > (nrow(junction) - 2)) curve_par = gpar(lwd=1.5, col="red")
		pcol <- "black"
		if (i > (nrow(junction) - 2)) pcol <- "red"
		if (i%%2 == 0) {
			ymid = -runif(1, 0.1, 0.3)*max(density$y)
			ymin = min(ymin, ymid*1.1)
			gp = gp+annotation_custom(xsplineGrob(c(0, 0, 1, 1), c(1, 0, 0, 0), shape=1, gp=curve_par), j[1], xmid, 0, ymid)+
				annotation_custom(xsplineGrob(c(1, 1, 0, 0), c(1, 0, 0, 0), shape=1, gp=curve_par), xmid, j[2], 0, ymid)
		} else {
			ymid = runif(1, 1.2, 1.4)*max(j[3:4])
			ymax = max(ymax, ymid*1.1)
			gp = gp+annotation_custom(xsplineGrob(c(0, 0, 1, 1), c(0, 1, 1, 1), shape=1, gp=curve_par), j[1], xmid, j[3], ymid)+
				annotation_custom(xsplineGrob(c(1, 1, 0, 0), c(0, 1, 1, 1), shape=1, gp=curve_par), xmid, j[2], j[4], ymid)
		}
		gp = gp+annotate("label", x=xmid, y=ymid, label=as.character(j[5]), color=pcol, 
			vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"), size=5, border.color="white")
	}
	gp <- gp+scale_y_continuous(breaks=c(0, max(density$y)), limits=c(ymin, ymax))+coord_cartesian(clip="off")
	#ggsave(plot=gp, width=10, height=4, dpi=200, "test.png", limitsize=F)
	return(gp)
})
pjcb <- ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
	geom_intron(data=info_introns, aes(color=transcript_name, strand=strand), linewidth=0.5, arrow.min.intron.length=300)+
	geom_range(aes(fill=transcript_name, color=transcript_name), linewidth=0, height=0.3)+
	geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name, fill=transcript_name, 
	color=transcript_name), linewidth=0.5, height=0.1)+
	labs(title=NULL, x=NULL, y=NULL)+
	scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, labels=info_ranges, expand=c(0.01, 0.01))+
	#geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), 
	#x=x_min, hjust=0, vjust=0, nudge_y=0.3, size=0.3*text_size)+
	scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
	axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
	axis.text.x=element_text(size=text_size, color="black"), legend.position="none")+tag_thm
pjc <- wrap_elements(wrap_plots(c(pls, list(pjcb)), ncol=1, heights=rep(c(3, 3, 2), length(gene_list)))+
	plot_annotation(gene, theme=theme(plot.title=element_text(size=title_size, colour="black", hjust=0.5), 
	plot.margin=margin(0,0,0,0), panel.spacing=unit(0, "pt"))))+tag_thm

rec_mm <- data.frame(Group=rep(c("D>V", "A>P"), each=2), Type=rep(c("AS", "ATSS"), 2), 
	Number=c(table(cmp_vdv_raw$ET[which(cmp_vdv_raw$Group == "IDG")]),
	table(cmp_vap_raw$ET[which(cmp_vap_raw$Group == "IDG")])))
rec_mm$Group <- factor(rec_mm$Group, levels=c("D>V", "A>P"))
pg <- ggplot(rec_mm, aes(x=Group, y=Number, fill=Type))+
	geom_bar(stat="identity", width=0.6, position=position_dodge(0.8))+
	labs(title=NULL, x=NULL, y="Number", fill="Isoform")+
	scale_fill_manual(values=col_list[c(1,3)])+
	scale_y_continuous(expand=c(0, 0))+
	guides(fill=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, color='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(color="black"), 
	axis.text.y=element_text(size=text_size, color="black"), 
	axis.text.x=element_text(size=title_size, color="black"), 
	axis.title=element_text(size=title_size, color="black"), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())+tag_thm
ids <- which(cmp_vdv_raw$Group == "IDG" & cmp_vdv_raw$ET != "Diff")
rec_mk <- data.frame(Group="D>V", Total=length(ids), table(cmp_vdv_raw$ETS[ids]))
ids <- which(cmp_vap_raw$Group == "IDG" & cmp_vap_raw$ET != "Diff")
rec_mk <- rbind(rec_mk, data.frame(Group="A>P", Total=length(ids), table(cmp_vap_raw$ETS[ids])))
colnames(rec_mk) <- c("Group", "Total", "Type", "Number")
rec_mk$Group <- factor(rec_mk$Group, levels=c("D>V", "A>P"))
rec_mk$Type <- factor(rec_mk$Type, levels=c("SE", "A3SS", "A5SS"))
rec_mk$Rate <- rec_mk$Number*100/rec_mk$Total
ph <- ggplot(rec_mk, aes(x=Group, y=Rate, fill=Type))+geom_bar(stat="identity", width=0.6)+
	labs(title=NULL, x=NULL, y="Percentage (%)", fill="AS")+
	scale_fill_manual(values=col_list[c(2,5,6)])+
	scale_y_continuous(expand=c(0, 0))+
	guides(fill=guide_legend(byrow=T))+
	theme(axis.line=element_line(linetype=1, color='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(color="black"), 
	axis.text.y=element_text(size=text_size, color="black"), 
	axis.text.x=element_text(size=title_size, color="black"), 
	axis.title=element_text(size=title_size, color="black"), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank())+tag_thm
pec <- wrap_plots(list(pg, ph), nrow=1)+plot_layout(guides="collect")+tag_thm

pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pc, pca, pcb), nrow=1, widths=c(1.11, 1, 1, 1))+plot_annotation(tag_levels=list(c("A", "C", "D", "")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pb, pd, pda, pdb), nrow=1, widths=c(1, 1, 1, 1))+plot_annotation(tag_levels=list(c("B", "E", "F", "")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pea, peb, pec), nrow=1, widths=c(1, 1, 1.2))+
	plot_annotation(tag_levels=list(c("G", "H", "I", "")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pia, pib, pic), nrow=1, widths=c(3, 2.2, 2.3))+
	plot_annotation(tag_levels=list(c("J", "", "")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pja, pjb, pjc), nrow=1, widths=c(3, 2.2, 2.3))+
	plot_annotation(tag_levels=list(c("K", "", "")), theme=tag_thm))+tag_thm), 
	ncol=1, heights=c(1, 1, 1, 1, 1)), width=13, height=16, dpi=200, filename="oe_fl_fig_03.png", limitsize=F)



