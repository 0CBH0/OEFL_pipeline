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
library(org.Mm.eg.db)
library(GO.db)
options(stringsAsFactors=FALSE)

col_list <- c(c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8"), brewer.pal(12,"Set3")[-c(2, 9)])
text_size <- 13
title_size <- 15
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, color="black"), plot.margin=margin(-3,-3,-3,-3), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))

osn_rna <- readRDS("osn_rna2.rds")
#types <- names(table(osn_rna[["cell.cls"]][,1]))
#types <- lapply(types, function(x){FindMarkers(osn_rna, ident.1=x, group.by="cell.cls")})
#names(types) <- names(table(osn_rna[["cell.cls"]][,1]))
#for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
#for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
#for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
osn_rna$cell.subtype_fix <- "mOSN"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 36)] <- "iOSN"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 39)] <- "iOSN"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 40)] <- "Basophil"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 41)] <- "Basophil"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 42)] <- "Basophil"
#ggsave(plot=UMAPPlot(osn_rna, group.by="cell.subtype_fix", pt.size=1.2, label=T, label.size=5)+
#	labs(title=NULL, x="UMAP1", y="UMAP2", color="Type"), width=10, height=8, dpi=200, "test.png")

types <- c("Basophil", "iOSN", "mOSN")
cell_col <- col_list[1:length(types)]
names(cell_col) <- types

osn_nano <- readRDS("osn_nano.rds")
#types <- names(table(osn_nano[["cell.cls"]][,1]))
#types <- lapply(types, function(x){FindMarkers(osn_nano, ident.1=x, group.by="cell.cls")})
#names(types) <- names(table(osn_nano[["cell.cls"]][,1]))
#for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
#for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
#for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
#table(osn_rna$cell.cls[match(colnames(osn_nano)[which(osn_nano$cell.cls == 30)], colnames(osn_rna))]) # 36
#table(osn_rna$cell.cls[match(colnames(osn_nano)[which(osn_nano$cell.cls == 33)], colnames(osn_rna))]) # xxx
#table(osn_rna$cell.cls[match(colnames(osn_nano)[which(osn_nano$cell.cls == 36)], colnames(osn_rna))]) # 39
#table(osn_rna$cell.cls[match(colnames(osn_nano)[which(osn_nano$cell.cls == 38)], colnames(osn_rna))]) # 40 42
#table(osn_rna$cell.cls[match(colnames(osn_nano)[which(osn_nano$cell.cls == 39)], colnames(osn_rna))]) # 41
osn_nano$cell.subtype_fix <- "mOSN"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 30)] <- "iOSN"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 33)] <- "iOSN"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 36)] <- "iOSN"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 38)] <- "Basophil"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 39)] <- "Basophil"
#ggsave(plot=UMAPPlot(osn_nano, group.by="cell.subtype_fix", pt.size=1.2, label=T, label.size=5)+
#	labs(title=NULL, x="UMAP1", y="UMAP2", color="Type"), width=10, height=8, dpi=200, "test.png")

cs_info <- data.frame(osn_rna@meta.data)
pa <- wrap_elements(ggplot(cs_info, aes(x="Term", y=nFeature_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Genes", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	axis.line.x=element_blank(), axis.line.y=element_line(color="black"), panel.background=element_blank()))+tag_thm
pb <- wrap_elements(ggplot(cs_info, aes(x="Term", y=nCount_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Counts", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	axis.line.x=element_blank(), axis.line.y=element_line(color="black"), panel.background=element_blank()))+tag_thm
pc <- wrap_elements(ggplot(cs_info, aes(x="Term", y=cell.mt))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Mito. (%)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	axis.line.x=element_blank(), axis.line.y=element_line(color="black"), panel.background=element_blank()))+tag_thm

res_info <- read.delim("ass_type_info_osn.tsv", h=T)
res_info$sample <- "ONT_OSN"
res_info$count[1] <- res_info$count[1] + res_info$count[2]
res_info <- res_info[-2,]
res_info$rate <- res_info$count*100/sum(res_info$count)
res_info$type <- factor(res_info$type, levels=rev(res_info$type), labels=rev(c("FSM", "ISM3'", "ISM5'", "ISM internal", "Mono exonic")))
pd <- wrap_elements(ggplot(res_info, aes(y=type, x=rate))+
	geom_bar(stat="identity", position=position_dodge(0.8), fill=col_list[1])+
	labs(title=NULL, y=NULL, x="Percentage of reads (%)")+
	scale_x_continuous(expand=c(0, 0))+
	theme(legend.position="none", axis.line=element_line(linetype=1, color='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.y=element_blank(), axis.ticks.x=element_line(color="black"), 
	axis.text=element_text(size=text_size, color="black"), 
	axis.title=element_text(size=title_size, color="black"), 
	legend.title=element_text(size=title_size, color="black", face="bold"), 
	legend.text=element_text(size=text_size, color="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

terms <- read.csv("sce_nano_len_sub.csv", r=1)[, 1]
cell_info <- read.csv("sce_nano_info.csv", r=1, h=T)
pe <- ggplot(, aes(x=terms, fill=type, color=type))+geom_density(fill=col_list[1], color=col_list[1], alpha=0.7, linewidth=0.8)+
	scale_y_continuous(expand=c(0, 0))+scale_x_continuous(breaks=seq(0, 3000, 500), limits=c(0, 2000), expand=c(0, 0))+
	labs(title=NULL, x="Length of tagged reads", y="Density\n")+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), axis.title=element_text(size=title_size, color="black"), 
	axis.ticks.y=element_blank(), axis.ticks.x=element_line(color="black"), 
	axis.text.y=element_blank(), axis.text.x=element_text(size=text_size, color="black"), axis.line=element_line(color="black"))+tag_thm
pf <- ggplot(cell_info, aes(x="Term", y=Len))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Length (mean)", y=NULL)+scale_x_discrete(breaks=NULL)+theme(
	axis.title=element_text(size=title_size, color="black"), axis.text=element_text(size=text_size, color="black"), 
	axis.line.x=element_blank(), axis.line.y=element_line(color="black"), panel.background=element_blank())+tag_thm
pg <- ggplot(cell_info, aes(x="Term", y=Isoform))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Isoforms", y=NULL)+scale_x_discrete(breaks=NULL)+theme(
	axis.title=element_text(size=title_size, color="black"), axis.text=element_text(size=text_size, color="black"), 
	axis.line.x=element_blank(), axis.line.y=element_line(color="black"), panel.background=element_blank())+tag_thm

osn_rna$cell.subtype_fix <- factor(osn_rna$cell.subtype_fix, levels=types)
rec_da <- data.frame(X=osn_rna@reductions$umap@cell.embeddings[,2], 
	Y=-osn_rna@reductions$umap@cell.embeddings[,1], Type=osn_rna$cell.subtype_fix)
pda <- ggplot(rec_da, aes(x=X, y=Y, color=Type))+geom_point(size=1)+
	labs(title="Identified with illumina", x=NULL, y="UMAP2", color=NULL)+
	scale_color_manual(values=cell_col, drop=F)+
	guides(color=guide_legend(nrow=1, override.aes=list(size=4), reverse=T))+
	theme(legend.position="none", axis.line=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=1), 
	axis.text.x=element_text(size=text_size, color="black"), 
	axis.text.y=element_text(size=text_size, color="black", angle=90, hjust=0.5, vjust=0.5), 
	axis.title=element_text(size=title_size, color="black"), 
	plot.title=element_text(size=title_size, hjust=0.8, color="black"), 
	legend.box.spacing=unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.text=element_text(size=text_size, color="black"), 
	legend.title=element_text(size=title_size, color="black"), 
	panel.background=element_rect(0, linetype=0))+tag_thm
osn_nano$cell.subtype_fix <- factor(osn_nano$cell.subtype_fix, levels=types)
rec_db <- data.frame(Y=osn_nano@reductions$umap@cell.embeddings[,2], 
	X=osn_nano@reductions$umap@cell.embeddings[,1], Type=osn_nano$cell.subtype_fix)
pdb <- ggplot(rec_db, aes(x=X, y=Y, color=Type))+geom_point(size=1)+
	labs(title="Identified with nanopore", x="UMAP1", y="UMAP2", color=NULL)+
	scale_color_manual(values=cell_col, drop=F)+
	guides(color=guide_legend(nrow=1, override.aes=list(size=4), reverse=T))+
	theme(legend.position="bottom", panel.border=element_rect(color="black", fill=NA, linewidth=1), 
	axis.text.x=element_text(size=text_size, color="black"), 
	axis.text.y=element_text(size=text_size, color="black", angle=90, hjust=0.5, vjust=0.5), 
	axis.title=element_text(size=title_size, color="black"), 
	plot.title=element_text(size=title_size, hjust=0.8, color="black"), 
	legend.box.spacing=unit(0, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))+tag_thm
pdx <- wrap_elements(wrap_plots(list(pda, pdb), ncol=1))+tag_thm

rec_c <- data.frame()
for (i in 1:length(types))
{
	for (j in 1:length(types))
	{
		cells_rna <- colnames(osn_rna)[which(osn_rna$cell.subtype_fix == types[i])]
		cells_nano <- colnames(osn_nano)[which(osn_nano$cell.subtype_fix == types[j])]
		if (length(cells_rna) == 0 | length(cells_nano) == 0) next
		rec_c <- rbind(rec_c, data.frame(X=types[i], Y=types[j], 
			Val=length(intersect(cells_rna, cells_nano))*100/length(union(cells_rna, cells_nano))))
	}
}
rec_c$X <- factor(rec_c$X, levels=types)
rec_c$Y <- factor(rec_c$Y, levels=types)
pca <- ggplot(rec_c, aes(x=X, y=Y, fill=Val))+geom_tile()+
	labs(title="Consistency of cell types", x="Identified with illumina", y="Identified with nanopore", fill="%")+
	scale_fill_viridis()+scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.ticks=element_blank(), 
	strip.text=element_blank(), strip.background=element_blank(), 
	axis.text.x=element_text(size=text_size, color="black"), axis.text.y=element_text(size=text_size, color="black", angle=90, hjust=0.5), 
	axis.title=element_text(size=title_size, color="black"), 
	plot.title=element_text(size=title_size, color="black", hjust=0.5), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"))+tag_thm

#cmp_vdv_raw <- read.csv("cmp_vdv_filter3_oe.csv", h=T, r=1)
#cmp_vap_raw <- read.csv("cmp_vap_filter3_oe.csv", h=T, r=1)
#for (type in c("A > P", "D > V"))
#{
#	terms <- cmp_vap_raw$Symbol[which(cmp_vap_raw$Group == "IDG")]
#	if (type == "D > V") terms <- cmp_vdv_raw$Symbol[which(cmp_vdv_raw$Group == "IDG")]
#	eid <- mapIds(org.Mm.eg.db, keys=terms, column="ENTREZID", keytype="SYMBOL", multiVals="first")
#	ego <- enrichGO(gene=eid, keyType="ENTREZID", OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH", readable=T)
#	ego <- ego@result[which((ego@result$pvalue < 0.01 | ego@result$qvalue < 0.05) & ego@result$Count > 1),]
#	if (nrow(ego) == 0) next
#	ego_id <- data.frame(ID=ego$ID, level=0)
#	id_filtered <- c()
#	for (i in 1:nrow(ego_id))
#	{
#		if (!is.na(match(ego_id$ID[i], id_filtered))) next
#		rec <- as.character(unlist(mget(ego_id$ID[i], GOBPPARENTS, ifnotfound=NA)))
#		parents <- rec
#		level <- 0
#		while(length(parents) > 0)
#		{
#			parents_ori <- parents
#			parents <- c()
#			for (p in parents_ori)
#			{
#				if (p == "all" & ego_id$level[i] == 0) ego_id$level[i] <- level
#				if (is.na(p) | p == "all") next
#				rec <- c(rec, p)
#				parents <- c(parents, as.character(unlist(mget(p, GOBPPARENTS, ifnotfound=NA))))
#			}
#			level <- level + 1
#			parents <- unique(parents)
#		}
#		rec <- unique(rec)
#		rec <- which(!is.na(match(ego_id$ID, rec)))
#		if (length(rec) > 0) id_filtered <- c(id_filtered, ego_id$ID[rec])
#	}
#	ego_id$level[match(unique(id_filtered), ego_id$ID)] <- 0
#	if (length(which(ego_id$level > 2)) == 0) next
#	ego_id <- ego_id[which(ego_id$level > 2),, drop=F]
#	ego <- ego[match(ego_id$ID, ego$ID),, drop=F]
#	ego$Level <- ego_id$level
#	ego$Type <- type
#	ego <- ego[order(ego$pvalue),]
#	ego <- ego[order(ego$Count, decreasing=T),]
#	if (type == "D > V") write.csv(ego, "cmp_dv_raw_go2.csv")
#	if (type == "A > P") write.csv(ego, "cmp_ap_raw_go2.csv")
#}
ego_ap <- read.csv("cmp_ap_raw_go2.csv", r=1, h=T)
ego_dv <- read.csv("cmp_dv_raw_go2.csv", r=1, h=T)
ego_total <- rbind(ego_dv[1:min(nrow(ego_dv), 10),], ego_ap[1:min(nrow(ego_ap), 10),])
#ego_total <- rbind(ego_dv[1:min(nrow(ego_dv), 10),], ego_ap[c(2:7, 10:13),])
ego_total$Rank <- factor(rev(1:nrow(ego_total)))
ego_total$Type <- factor(ego_total$Type, levels=c("D > V", "A > P"))
pga <- wrap_elements(ggplot(ego_total, aes(x=Count, y=Rank, colour=pvalue, size=Count))+geom_point()+
	scale_color_gradient(low="#440255", high="#FFFFBF", guide=guide_colorbar(reverse=T))+
	#scale_color_continuous(low=brewer.pal(11,"Spectral")[1], high=brewer.pal(11,"Spectral")[6], 
	#guide=guide_colorbar(reverse=TRUE))+
	scale_y_discrete(breaks=ego_total$Rank, labels=ego_total$Description, position="right")+
	facet_grid(Type~., scales="free_y", space="free_y", switch="y")+
	labs(title=NULL, x=NULL, y=NULL, color="p.val")+
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(colour="#800026", fill="#800026"), 
	strip.text=element_text(size=title_size*0.8, colour="white", face="bold"), 
	panel.spacing=unit(5, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	axis.line=element_blank(), panel.background=element_rect(fill='gray98')))+tag_thm

ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pb, pc, pd), nrow=1, widths=c(1,1,1,2))+
	plot_annotation(tag_levels=list(c("A", "B", "C", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pe, pf, pg, pca), nrow=1, widths=c(1.4,0.8,0.8,1.5))+
	plot_annotation(tag_levels=list(c("E", "F", "G", "H")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pdx, pga), widths=c(1,3.5))+
	plot_annotation(tag_levels=list(c("I", "", "J"))))+tag_thm), 
	ncol=1, heights=c(1,1,1.5))+tag_thm, width=13, height=14, dpi=200, filename="oe_fl_fig_S03.png", limitsize=F)













tt <- names(table(osn_rna[["Type"]][,1]))
tt <- lapply(tt, function(x){FindMarkers(osn_rna, ident.1=x, group.by="Type")})
names(tt) <- names(table(osn_rna[["Type"]][,1]))
for (i in 1:length(tt)) tt[[i]] <- tt[[i]][which(tt[[i]]$p_val_adj < 0.01 & tt[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(tt)) tt[[i]] <- tt[[i]][order(tt[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(tt)) tt[[i]] <- tt[[i]][order(tt[[i]]$p_val_adj),]



osn_nano <- readRDS("osn_nano.rds")

tt <- names(table(osn_nano[["cell.cls"]][,1]))
tt <- lapply(types, function(x){FindMarkers(osn_nano, ident.1=x, group.by="cell.cls")})
#names(types) <- names(table(osn_nano[["cell.cls"]][,1]))


ggsave(plot=UMAPPlot(osn_nano, group.by="cell.cls", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", color=NULL), width=12, height=10, dpi=200, "test.png")

tt <- data.frame()
for (c in names(types))
{
	ts <- rownames(types[[c]])[which(types[[c]]$p_val_adj < 0.001)[1:20]]
	rec <- data.frame()
	for (m in names(markers)) rec <- rbind(rec, data.frame(Type=m, 
		Count=length(intersect(ts, rownames(markers[[m]])[which(markers[[m]]$p_val_adj < 0.001)[1:20]]))))
	tt <- rbind(tt, data.frame(rec[which.max(rec$Count),], ID=c))
}

head(types[["30"]])

30 iOSN
33 iOSN
36
38
39

sce_rna <- readRDS("OER_fix.rds")
markers <- names(table(sce_rna[["cell.subtype_fix"]][,1]))
markers <- lapply(markers, function(x){FindMarkers(sce_rna, ident.1=x, group.by="cell.subtype_fix")})
names(markers) <- names(table(sce_rna[["cell.subtype_fix"]][,1]))

sce_rna

osn_rna

terms <- intersect(rownames(sce_rna[["SCT"]]@scale.data), rownames(osn_rna[["SCT"]]@scale.data))

cc <- cor(sce_rna[["SCT"]]@scale.data[terms,], osn_rna[["SCT"]]@scale.data[terms,])

AverageExpression



for (name)

osn_rna <- readRDS("osn_rna2.rds")
sce_rna <- readRDS("OER_fix.rds")


sce_terms <- AverageExpression(sce_rna, assays="SCT", group.by="cell.subtype_fix")$SCT
osn_terms <- AverageExpression(osn_rna, assays="SCT", group.by="cell.cls")$SCT
terms <- intersect(rownames(sce_terms)[which(rowSums(sce_terms) > 0)], rownames(osn_terms)[which(rowSums(osn_terms) > 0)])
cc <- cor(as.matrix(sce_terms[terms,]), as.matrix(osn_terms[terms,]))
rec <- data.frame()
for (i in 1:ncol(cc)) rec <- rbind(rec, data.frame(ID=colnames(cc)[i], Type=rownames(cc)[which.max(cc[, i])], Cor=cc[which.max(cc[, i]), i]))


markers <- readRDS("sce_markers_type_4.rds")

types <- c("Monocyte", "Macrophage", "BC", "RC", "Microvillar", "Neutrophils", "Basophil", 
	"Osteogenic", "Pericytes", "Ensheathing", "Sustentacular", "Bowman", 
	"HBC", "GBC", "INP", "Immature", "Mature")
cell_col <- col_list[1:length(types)]
names(cell_col) <- types
types_osn <- c("HBC", "GBC", "INP", "Immature", "Mature")
sce_rna <- readRDS("osn_rna2_fix.rds")
sce_nano <- readRDS("oe_nano_fix.rds")

cs_info <- data.frame(sce_rna@meta.data)
pa <- wrap_elements(ggplot(cs_info, aes(x="Term", y=nFeature_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Genes", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(color="black"), axis.line=element_line(color="black"), panel.background=element_blank()))+tag_thm
pb <- wrap_elements(ggplot(cs_info, aes(x="Term", y=nCount_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Counts", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(color="black"), axis.line=element_line(color="black"), panel.background=element_blank()))+tag_thm
pc <- wrap_elements(ggplot(cs_info, aes(x="Term", y=cell.mt))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Mito. (%)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(color="black"), axis.line=element_line(color="black"), panel.background=element_blank()))+tag_thm

res_info <- read.delim("ass_type_info_osn.tsv", h=T)
res_info$sample <- "ONT_OSN"
res_info$count[1] <- res_info$count[1] + res_info$count[2]
res_info <- res_info[-2,]
res_info$rate <- res_info$count*100/sum(res_info$count)
res_info$type <- factor(res_info$type, levels=rev(res_info$type), labels=rev(c("FSM", "ISM3'", "ISM5'", "ISM internal", "Mono exonic")))
pd <- wrap_elements(ggplot(res_info, aes(x=type, y=rate))+
	geom_bar(stat="identity", position=position_dodge(0.8), fill=col_list[1])+
	labs(title=NULL, x=NULL, y="Reads percentage", fill="Type")+
	#scale_fill_manual(values=rev(c("#db5f56", "#8ddd3e", "#fdb54e", "#7dcac1", "#5684da")), guide=guide_legend(reverse=T), drop=F)+
	scale_y_continuous(expand=c(0, 0))+
	theme(legend.position="none", axis.line=element_line(linetype=1, color='black'), panel.background=element_rect(0, linetype=0), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(color="black"), 
	axis.text.y=element_text(size=text_size, color="black"), 
	#axis.text.x=element_text(size=text_size, color="black", angle=-45, vjust=1, hjust=0), 
	axis.text.x=element_text(size=text_size*0.8, color="black"), 
	axis.title=element_text(size=title_size, color="black"), 
	legend.title=element_text(size=title_size, color="black", face="bold"), 
	legend.text=element_text(size=text_size, color="black"), 
	legend.key.size=unit(20, "pt"), legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, design="ABCD", widths=c(1,1,1,2)), 
	width=10, height=4, dpi=200, filename="oe_fl_fig_S032.png", limitsize=F)



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
pe <- wrap_elements(ggplot(, aes(x=terms, fill=type, color=type))+geom_density(fill=col_list[1], color=col_list[1], alpha=0.7, linewidth=0.8)+
	scale_y_continuous(expand=c(0, 0))+scale_x_continuous(breaks=seq(0, 3000, 500), limits=c(0, 2000), expand=c(0, 0))+
	labs(title=NULL, x="Length of tagged reads", y="Density\n")+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	axis.text=element_text(size=text_size, color="black"), axis.line=element_line(color="black")))+tag_thm
pf <- wrap_elements(ggplot(cell_info, aes(x="Term", y=Isoform))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Isoforms", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(color="black"), axis.line=element_line(color="black"), panel.background=element_blank()))+tag_thm
pg <- wrap_elements(ggplot(cell_info, aes(x="Term", y=Len))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Length (mean)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(color="black"), axis.line=element_line(color="black"), panel.background=element_blank()))+tag_thm

ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, D=pd, design="ABCD", widths=c(2,2,2,3)), 
	width=8, height=4, dpi=200, filename="oe_fl_fig_S032.png", limitsize=F)

