library(dplyr)
library(future)
library(scater)
library(Seurat)
library(Matrix)
library(edgeR)
library(ggplot2)
library(patchwork)
library(BiocParallel)
library(scDblFinder)
library(ggbeeswarm)
library(RColorBrewer)
library(Biostrings)
library(DropletUtils)
library(viridis)
library(ggpointdensity)
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
font_scale <- 2
text_size <- 6*font_scale
title_size <- 7.5*font_scale
tag_thm <- theme(plot.tag=element_text(size=10*font_scale, colour="black", face="bold"), plot.margin=margin())


#/md01/yangjw28/data/OSN/jointF18w230419/jointF18w0509/outs
sce <- CreateSeuratObject(Read10X_h5("filtered_feature_bc_matrix.h5")[1], project="OSNnano")
sce[["cell.sample"]] <- Idents(sce)
sce[["cell.mt"]] <- PercentageFeatureSet(sce, features=rownames(sce)[grep("^mt-", rownames(sce))])
sce <- as.SingleCellExperiment(sce)
sce <- scDblFinder(sce, samples="cell.sample", BPPARAM=MulticoreParam(16))

cs_info <- data.frame(gene=colSums(counts(sce) > 0), gi="F", count=colSums(counts(sce)), ci="F", mt=sce$cell.mt, mti="F", dbi=sce$scDblFinder.class)
cs_info$mti[which(cs_info$mt < 10)] <- "P"
pa <- ggplot(cs_info, aes(x="Term", y=mt, colour=mti))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Mito. (%)", y=NULL, color="Type")+scale_y_continuous(limits=c(0, 100), breaks=seq(10, 90, 10))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=10)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
cs_limit <- c(500, 5000)
cs_info$gi[which(cs_info$gene > cs_limit[1] & cs_info$gene < cs_limit[2])] <- "P"
pb <- ggplot(cs_info, aes(x="Term", y=gene, colour=gi))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Genes", y=NULL, color="Type")+scale_y_continuous(breaks=c(cs_limit[1], seq(1000, max(cs_info$gene), 1000)))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=cs_limit)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
cs_limit <- c(500, 15000)
cs_info$ci[which(cs_info$count > cs_limit[1] & cs_info$count < cs_limit[2])] <- "P"
pc <- ggplot(cs_info, aes(x="Term", y=count, color=ci))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+#guides(colour="none")+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	labs(title=NULL, x="Counts", y=NULL, color="Type")+scale_y_continuous(breaks=c(cs_limit[1], seq(2000, max(cs_info$count), 2000)))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=cs_limit)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())

ggsave(plot=wrap_plots(list(pa, pb, pc), nrow=1)+plot_annotation(
	title=paste0("#cells after screening: ", length(which(cs_info$mti == "P" & cs_info$gi == "P" & cs_info$dbi == "singlet"))),
	theme=theme(plot.title=element_text(size=title_size, hjust=0.5))), width=7, height=6, dpi=200, "qc_test_osn_rna.png")

sce <- sce[, which(cs_info$mti == "P" & cs_info$gi == "P" & cs_info$dbi == "singlet")]
sce <- addPerCellQC(sce, percent_top=c(20, 50, 100, 200))
sce$total_features <- sce$detected
sce$log10_total_features <- log10(sce$detected)
sce$total_counts <- sce$sum
sce$log10_total_counts <- log10(sce$sum)
sce$featcount_ratio <- sce$log10_total_counts/sce$log10_total_features
mod <- loess(colData(sce)$log10_total_features~colData(sce)$log10_total_counts)
pred <- predict(mod, newdata=data.frame(log10_total_counts=colData(sce)$log10_total_counts))
sce$featcount_dist <- colData(sce)$log10_total_features - pred
sce$pct_counts_top_20_features <- colData(sce)[[intersect(c("percent_top_20","pct_counts_in_top_20_features","percent.top_20"), colnames(colData(sce)))[[1]]]]
sce$pct_counts_top_50_features <- colData(sce)[[intersect(c("percent_top_50","pct_counts_in_top_50_features","percent.top_50"), colnames(colData(sce)))[[1]]]]
vars <- c("log10_total_counts:both:5", "log10_total_features:both:5", "pct_counts_top_20_features:both:5", "featcount_dist:both:5")
out <- table(unlist(lapply(strsplit(vars,":"), function(f){which(isOutlier(sce[[f[1]]], log=FALSE, nmads=as.numeric(f[3]), type=f[2]))})))
out <- as.numeric(names(out)[which(out > 0)])
fc <- data.frame(Feature=colData(sce)$log10_total_features, Count=colData(sce)$log10_total_counts, Group="Pass")
rownames(fc) <- rownames(colData(sce))
if (length(out) > 0) fc$Group[out] <- "Filter"
fc <- fc[order(fc$Group),]
pc <- ggplot(fc, aes(x=Feature, y=Count, colour=Group))+geom_point(stroke=0, shape=16, size=2, alpha=0.5)+
	stat_smooth(method=loess, se=FALSE, colour="black")+
	labs(title=NULL, x="Feature (log10)", y="Count (log10)")+
	scale_colour_manual(values=col_list[c(3, 1)])+
	guides(colour=guide_legend(override.aes=list(size=5)))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))
ggsave(plot=wrap_plots(list(pb, pc), nrow=1, widths=c(1, 6))+
	plot_annotation(title=paste0("The QC results (", length(fc$Group == "Pass"), "/", nrow(cs_info), ")"), 
	theme=theme(plot.title=element_text(size=16, hjust=0.5))), 
	width=8, height=6, dpi=200, "qc_1_osn_rna.png")

if (length(out) > 0) sce <- sce[, -out]
sce <- sce[rowSums(counts(sce)) > 0,]
sce <- as.Seurat(sce, data=NULL)
names(sce@assays) <- "RNA"
DefaultAssay(sce) <- "RNA"
sce <- SCTransform(sce, method="glmGamPoi")
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims=1:20)
#sce@reductions[["umap"]]@cell.embeddings <- -sce@reductions[["umap"]]@cell.embeddings
ggsave(plot=UMAPPlot(sce, group.by="cell.sample", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL), width=12, height=10, dpi=200, "qc_2.png")
sce <- FindNeighbors(sce, reduction="umap", dims=1:2)
sce <- FindClusters(sce, algorithm=4, method="igraph", resolution=2)
sce[["cell.cls"]] <- Idents(sce)
ggsave(plot=UMAPPlot(sce, group.by="cell.cls", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL), width=12, height=10, dpi=200, "qc_2_osn_rna.png")

cell_marker <- read.delim("cellmarker_sel.txt")
marker_list <- intersect(unique(unlist(strsplit(cell_marker$Marker, split=","))), rownames(sce))
marker_list <- data.frame(Term=marker_list, ID=marker_list, Type="")
for (i in 1:nrow(marker_list)) marker_list$Type[i] <- cell_marker$Type[grep(marker_list$Term[i], cell_marker$Marker)[1]]
pd_list <- lapply(1:nrow(marker_list), function(i)
{
	p <- FeaturePlot(sce, slot="data", features=marker_list$ID[i], min.cutoff=0, cols=c("gray98", "red2"))+
		labs(title=paste0(marker_list$Term[i], " (", marker_list$Type[i], ")"), x="UMAP1", y="UMAP2", colour="RNA")+
		theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
		axis.line=element_line(colour="black"), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
		legend.text=element_text(colour="black", size=12))
	p$data <- p$data[order(p$data[, 4]),]
	return(p)
})
ggsave(plot=wrap_plots(pd_list, ncol=6), width=36, height=25, dpi=200, filename="sce_umap_marker_osn_rna.png", limitsize=F)
bcs <- read.delim("osn_bc_info.tsv", h=F)
bcs <- intersect(bcs[, 1], colnames(sce))
sce <- subset(sce, cells=bcs)
sce <- SCTransform(sce, method="glmGamPoi")
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims=1:20)
saveRDS(sce, "osn_rna2.rds")

sce <- readRDS("osn_rna2.rds")
groups <- names(table(sce[["cell.cls"]][,1]))
groups <- lapply(groups, function(x){FindMarkers(sce, ident.1=x, group.by="cell.cls")})
names(groups) <- names(table(sce[["cell.cls"]][,1]))
for (i in 1:length(groups)) groups[[i]] <- groups[[i]][which(groups[[i]]$p_val_adj < 0.01 & groups[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(groups)) groups[[i]] <- groups[[i]][order(groups[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(groups)) groups[[i]] <- groups[[i]][order(groups[[i]]$p_val_adj),]
saveRDS(groups, "sce_markers_cls_osn_rna.rds")


sce <- subset(sce, cells=colnames(sce)[which(as.numeric(sce$cell.cls) < 39 & as.numeric(sce$cell.cls) != 36)])
sce <- SCTransform(sce, method="glmGamPoi")
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims=1:17)
ggsave(plot=UMAPPlot(sce, group.by="cell.cls", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL), width=12, height=10, dpi=200, paste0("qc_2_osn_rna_", 17, ".png"))
saveRDS(sce, "osn_rna2_fix.rds")


for (i in 16:30)
{
	sce <- RunUMAP(sce, dims=1:i)
	ggsave(plot=UMAPPlot(sce, group.by="cell.cls", pt.size=1, label=T, label.size=5)+
		labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL), width=12, height=10, dpi=200, paste0("qc_2_osn_rna_", i, ".png"))
}



#############################################################################################
library(dplyr)
library(future)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
font_scale <- 2
text_size <- 6*font_scale
title_size <- 7.5*font_scale
tag_thm <- theme(plot.tag=element_text(size=10*font_scale, colour="black", face="bold"), plot.margin=margin())

sce <- readRDS("osn_rna2.rds")
cell_marker <- read.delim("cellmarker_sel.txt")
type_list <- c(cell_marker$Type, "UNK")

ggsave(plot=UMAPPlot(sce, group.by="cell.cls", pt.size=1.2, label=T, label.size=5)+
	labs(title="Round 1", x="UMAP1", y="UMAP2", colour="Group"),
	width=10, height=8, dpi=200, "tttt.png")

sce_sub <- sce
types <- names(table(sce_sub[["cell.cls"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce_sub, ident.1=x, group.by="cell.cls")})
names(types) <- names(table(sce_sub[["cell.cls"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]

sce_sub$Type <- "Mature"
sce_sub$Type[which(sce_sub$cell.cls == 36)] <- "Immature"
sce_sub$Type[which(sce_sub$cell.cls == 39)] <- "Immature"
sce_sub$Type[which(sce_sub$cell.cls == 40)] <- "Basophil"
sce_sub$Type[which(sce_sub$cell.cls == 41)] <- "Basophil"
sce_sub$Type[which(sce_sub$cell.cls == 42)] <- "Basophil"

ggsave(plot=UMAPPlot(sce_sub, group.by="Type", pt.size=1.2, label=T, label.size=5)+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour="Type"),
	width=10, height=8, dpi=200, "ttttt.png")



for (i in 1:nrow(cell_marker)) 
{
	a = intersect(rownames(types[[42]]), unlist(strsplit(cell_marker$Marker[i], ",")))
	if (length(a) > 0) print(cell_marker$Type[i])
}


######################################################################
# Round 1
######################################################################
marker_list <- strsplit(cell_marker$Marker, ",")
names(marker_list) <- cell_marker$Type
for (i in 1:length(marker_list)) marker_list[[i]] <- match(marker_list[[i]], rownames(sce[["SCT"]]@data))

rec <- matrix(0, nrow=ncol(sce), ncol=length(marker_list), dimnames=list(colnames(sce), names(marker_list)))
for (i in 1:length(marker_list)) rec[, i] <- apply(sce[["SCT"]]@data[marker_list[[i]],, drop=F], 2, min)
res <- rep(0, nrow(rec))
for (i in 1:nrow(rec))
{
	id <- which(rec[i,] > 0)
	if (length(id) == 1) res[i] <- id
}
sce[["cell.subtype_fix"]] <- c("Unknow", names(marker_list))[res+1]
cell_filter <- data.frame()
for (i in 1:ncol(rec))
{
	terms <- which(res == i)
	if (length(terms) < 3) next
	terms <- data.frame(ID=terms, Val=rec[terms, i])
	terms <- terms[order(terms$Val, decreasing=T),]
	terms <- terms[2:min(50, nrow(terms)),]
	cell_filter <- rbind(cell_filter, data.frame(Cell=rownames(terms), Type=names(marker_list)[i]))
}
sce_sub <- subset(sce, cells=cell_filter$Cell)
filterList <- ""
groups <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
for (group in groups)
{
	umap_sub <- as.data.frame(sce_sub@reductions[["umap"]]@cell.embeddings[which(sce_sub[["cell.subtype_fix"]][, 1] == group),])
	umap_sub$Count <- 0
	umap_sub_dist <- as.matrix(dist(umap_sub))
	group_lim <- median(umap_sub_dist)*2
	for (i in 1:nrow(umap_sub_dist)) umap_sub$Count[i] <- length(which(umap_sub_dist[i,] < group_lim))
	if (length(filter) > 0) filterList <- c(filterList, rownames(umap_sub)[which(umap_sub$Count > nrow(umap_sub)*0.1)])
}
sce_sub <- subset(sce_sub, cells=filterList)

types <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce_sub, ident.1=x, group.by="cell.subtype_fix")})
names(types) <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "sce_markers_type_1.rds")

print("Round 1:")
print(table(sce_sub[["cell.subtype_fix"]]))
sce_sub[["cell.subtype_fix"]] <- factor(sce_sub[["cell.subtype_fix"]][, 1], levels=intersect(cell_marker$Type, sce_sub[["cell.subtype_fix"]][, 1]))
ggsave(plot=UMAPPlot(sce_sub, group.by="cell.subtype_fix", pt.size=1.2, label=T, label.size=5)+
	labs(title="Round 1", x="UMAP1", y="UMAP2", colour="Group"),
	width=10, height=8, dpi=200, "sce_umap_type_1.png")

######################################################################
# Round 2
######################################################################
types <- readRDS("sce_markers_type_1.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- match(rownames(types[[i]])[1:3], rownames(sce[["SCT"]]@data))
marker_list <- types

rec <- matrix(0, nrow=ncol(sce), ncol=length(marker_list), dimnames=list(colnames(sce), names(marker_list)))
for (i in 1:length(marker_list)) rec[, i] <- apply(sce[["SCT"]]@data[marker_list[[i]],, drop=F], 2, min)
res <- rep(0, nrow(rec))
for (i in 1:nrow(rec))
{
	id <- which(rec[i,] > 0)
	if (length(id) == 1) res[i] <- id
}
sce[["cell.subtype_fix"]] <- c("Unknow", names(marker_list))[res+1]
cell_filter <- data.frame()
for (i in 1:ncol(rec))
{
	terms <- which(res == i)
	if (length(terms) < 3) next
	terms <- data.frame(ID=terms, Val=rec[terms, i])
	terms <- terms[order(terms$Val, decreasing=T),]
	terms <- terms[2:min(50, nrow(terms)),]
	cell_filter <- rbind(cell_filter, data.frame(Cell=rownames(terms), Type=names(marker_list)[i]))
}
sce_sub <- subset(sce, cells=cell_filter$Cell)
filterList <- ""
groups <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
for (group in groups)
{
	umap_sub <- as.data.frame(sce_sub@reductions[["umap"]]@cell.embeddings[which(sce_sub[["cell.subtype_fix"]][, 1] == group),])
	umap_sub$Count <- 0
	umap_sub_dist <- as.matrix(dist(umap_sub))
	group_lim <- median(umap_sub_dist)*2
	for (i in 1:nrow(umap_sub_dist)) umap_sub$Count[i] <- length(which(umap_sub_dist[i,] < group_lim))
	if (length(filter) > 0) filterList <- c(filterList, rownames(umap_sub)[which(umap_sub$Count > nrow(umap_sub)*0.1)])
}
sce_sub <- subset(sce_sub, cells=filterList)

types <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce_sub, ident.1=x, group.by="cell.subtype_fix")})
names(types) <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "sce_markers_type_2.rds")

print("Round 2:")
print(table(sce_sub[["cell.subtype_fix"]]))
sce_sub[["cell.subtype_fix"]] <- factor(sce_sub[["cell.subtype_fix"]][, 1], levels=intersect(cell_marker$Type, sce_sub[["cell.subtype_fix"]][, 1]))
ggsave(plot=UMAPPlot(sce_sub, group.by="cell.subtype_fix", pt.size=1.2, label=T, label.size=5)+
	labs(title="Round 2", x="UMAP1", y="UMAP2", colour="Group"),
	width=10, height=8, dpi=200, "sce_umap_type_2.png")

######################################################################
# Round 3
######################################################################
term_count <- 50
types <- readRDS("sce_markers_type_1.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- types[[i]][1:3,]
marker_list_ori <- types
types <- readRDS("sce_markers_type_2.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.001),]
for (i in 1:length(types)) types[[i]] <- rbind(marker_list_ori[[names(types)[i]]], types[[i]][setdiff(rownames(types[[i]]), rownames(marker_list_ori[[names(types)[i]]])),])
for (i in 1:length(types)) term_count <- min(term_count, nrow(types[[i]]))
ref_info <- tibble(type=names(types), gene=list(""))
for (i in 1:length(types)) ref_info$gene[i] <- list(rownames(types[[i]])[1:min(term_count, nrow(types[[i]]))])
markers <- unique(unlist(ref_info$gene))

types <- readRDS("sce_markers_cls.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.001),]
pred <- data.frame()
for (type in types)
{
	terms <- intersect(rownames(type)[1:min(term_count, nrow(type))], markers)
	res <- enricher(terms, universe=markers, TERM2GENE=ref_info, minGSSize=1)
	term <- data.frame(Type="NA",  Count="0", FDR=1)
	if (length(res$ID) > 0)
	{
		term <- data.frame(Type=res$ID[1],  Count=res$GeneRatio[1], FDR=res$p.adjust[1])
		hyper <- cell_marker$Hyper[match(res$ID, cell_marker$Type)]
		ratio <- as.numeric(gsub("/.*", "", res$GeneRatio))/as.numeric(gsub(".*/", "", res$GeneRatio))
		for (i in 1:length(res$ID)) if (hyper[i] != hyper[1] & ratio[i]*2 > ratio[1]) term <- data.frame(Type="NA",  Count="0", FDR=1)
	}
	pred <- rbind(pred, term)
}
sce[["cell.subtype_fix"]] <- ""
for (i in 1:nrow(pred)) sce[["cell.subtype_fix"]][which(sce[["cell.cls"]][, 1] == rownames(pred)[i]), 1] <- pred$Type[i]
sce[["cell.hypertype"]] <- cell_marker$Hyper[match(sce[["cell.subtype_fix"]][, 1], cell_marker$Type)]
sce_sub <- subset(sce, cells=colnames(sce)[which(sce[["cell.subtype_fix"]][, 1] != "NA")])
types <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce_sub, ident.1=x, group.by="cell.subtype_fix")})
names(types) <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "sce_markers_type_3.rds")

term_count <- 50
types <- readRDS("sce_markers_type_1.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- types[[i]][1:3,]
marker_list_ori <- types
types <- readRDS("sce_markers_type_3.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.001),]
for (i in 1:length(types)) types[[i]] <- rbind(marker_list_ori[[names(types)[i]]], types[[i]][setdiff(rownames(types[[i]]), rownames(marker_list_ori[[names(types)[i]]])),])
for (i in 1:length(types)) term_count <- min(term_count, nrow(types[[i]]))
ref_info <- tibble(type=names(types), gene=list(""))
for (i in 1:length(types)) ref_info$gene[i] <- list(rownames(types[[i]])[1:min(term_count, nrow(types[[i]]))])
markers <- unique(unlist(ref_info$gene))

types <- readRDS("sce_markers_cls.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.001),]
pred <- data.frame(Type="", Count=0, FDR=0, Other="", Count.o="", FDR.o="")
for (type in types)
{
	terms <- intersect(rownames(type)[1:min(term_count, nrow(type))], markers)
	res <- enricher(terms, universe=markers, TERM2GENE=ref_info, minGSSize=1)
	if (length(res$ID) == 0)
	{
		pred <- rbind(pred, c("NA", "0", 1, "NA", "NA", "NA"))
		next
	}
	pred <- rbind(pred, c(res$ID[1], res$GeneRatio[1], res$p.adjust[1], paste(res$ID[-1], collapse=","), 
		paste(res$GeneRatio[-1], collapse=","), paste(res$p.adjust[-1], collapse=",")))
}
pred <- pred[-1,]
rownames(pred) <- 1:nrow(pred)
sce[["cell.subtype_fix"]] <- ""
for (i in 1:nrow(pred)) sce[["cell.subtype_fix"]][which(sce[["cell.cls"]][, 1] == rownames(pred)[i]), 1] <- pred$Type[i]
sce[["cell.hypertype"]] <- cell_marker$Hyper[match(sce[["cell.subtype_fix"]][, 1], cell_marker$Type)]
sce_sub <- subset(sce, cells=colnames(sce)[which(sce[["cell.subtype_fix"]][, 1] != "NA")])
saveRDS(pred, "sce_cls_type.rds")

print("Round 3:")
print(table(sce_sub[["cell.subtype_fix"]]))
sce_sub[["cell.subtype_fix"]] <- factor(sce_sub[["cell.subtype_fix"]][, 1], levels=intersect(cell_marker$Type, sce_sub[["cell.subtype_fix"]][, 1]))
ggsave(plot=UMAPPlot(sce_sub, group.by="cell.subtype_fix", pt.size=1.2, label=T, label.size=5)+
	labs(title="Round 3", x="UMAP1", y="UMAP2", colour="Group"),
	width=10, height=8, dpi=200, "sce_umap_type_3.png")

######################################################################
# Round 4
######################################################################
term_count <- 50
types <- readRDS("sce_markers_type_1.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@scale.data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- types[[i]][1:3,]
marker_list_ori <- types
types <- readRDS("sce_markers_type_2.rds")
for (i in 1:length(types)) types[[i]] <- types[[i]][intersect(rownames(types[[i]]), rownames(sce[["SCT"]]@scale.data)),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.001),]
for (i in 1:length(types)) types[[i]] <- rbind(marker_list_ori[[names(types)[i]]], types[[i]][setdiff(rownames(types[[i]]), rownames(marker_list_ori[[names(types)[i]]])),])
for (i in 1:length(types)) term_count <- min(term_count, nrow(types[[i]]))
pred <- readRDS("sce_cls_type.rds")

sce[["cell.subtype_ori"]] <- pred$Type[sce[["cell.cls"]][, 1]]
sce[["cell.subtype_fix"]] <- sce[["cell.subtype_ori"]][, 1]
sce[["cell.hypertype"]] <- cell_marker$Hyper[match(sce[["cell.subtype_fix"]][, 1], cell_marker$Type)]
sce_sub <- subset(sce, cells=colnames(sce)[which(sce[["cell.subtype_fix"]][, 1] != "NA")])

hyper_list <- table(cell_marker$Hyper)
hyper_list <- names(hyper_list)[which(hyper_list > 1)]
for (hyper in hyper_list)
{
	ci <- which(sce_sub[["cell.hypertype"]][, 1] == hyper)
	if (length(ci) < 3) next
	sce_hyper <- subset(sce_sub, cells=colnames(sce_sub)[ci])
	type_list <- intersect(cell_marker$Type[which(cell_marker$Hyper == hyper)], names(types))
	cls_list <- unique(sce_hyper[["cell.cls"]][, 1])
	
	tc <- 50
	for (type in type_list) tc <- min(term_count, nrow(types[[type]]))
	marker_list <- list()
	for (type in type_list) marker_list <- c(marker_list, list(match(rownames(types[[type]])[1:tc], rownames(sce_hyper[["SCT"]]@scale.data))))
	names(marker_list) <- type_list
	rec <- matrix(0, nrow=ncol(sce_hyper), ncol=length(marker_list), dimnames=list(colnames(sce_hyper), names(marker_list)))
	for (i in 1:length(marker_list)) rec[, i] <- apply(sce_hyper[["SCT"]]@scale.data[marker_list[[i]],], 2, sum)
	type_dect <- apply(rec, 1, function (x) 
	{
		tr <- which.max(x)
		if (length(which(x == x[tr])) > 1) return(0)
		return(tr)
	})
	sce_hyper[["cell.subtype_fix"]] <- c("Unknow", colnames(rec))[type_dect+1]
	sce_hyper[["cell.subtype_fix"]][which(sce_hyper[["cell.subtype_fix"]][, 1] == "Unknow"), 1] <- sce_hyper[["cell.subtype_ori"]][which(sce_hyper[["cell.subtype_fix"]][, 1] == "Unknow"), 1]
	for (cls in cls_list)
	{
		tr <- table(sce_hyper[["cell.subtype_fix"]][which(sce_hyper[["cell.cls"]][, 1] == cls), 1])
		td <- names(tr)[which.max(tr)]
		tr <- setdiff(names(tr), names(tr)[which(tr > max(tr)*0.05)])
		for (t in tr) sce_hyper[["cell.subtype_fix"]][which(sce_hyper[["cell.cls"]][, 1] == cls & sce_hyper[["cell.subtype_fix"]][, 1] == t), 1] <- td
	}
	sce_sub[["cell.subtype_fix"]][ci, 1] <- sce_hyper[["cell.subtype_fix"]][, 1]
}

types <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce_sub, ident.1=x, group.by="cell.subtype_fix")})
names(types) <- names(table(sce_sub[["cell.subtype_fix"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "sce_markers_type_4.rds")

print("Round 4:")
print(table(sce_sub[["cell.subtype_fix"]]))
sce_sub[["cell.subtype_fix"]] <- factor(sce_sub[["cell.subtype_fix"]][, 1], levels=intersect(cell_marker$Type, sce_sub[["cell.subtype_fix"]][, 1]))
ggsave(plot=UMAPPlot(sce_sub, group.by="cell.subtype_fix", pt.size=1.2, label=T, label.size=5)+
	labs(title="Round 4", x="UMAP1", y="UMAP2", colour="Group"),
	width=10, height=8, dpi=200, "sce_umap_type_4.png")


