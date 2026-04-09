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
source("GeomBarSignif.R")
options(stringsAsFactors=FALSE)

col_list <- c(c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8"), brewer.pal(12,"Set3")[-c(2, 9)])
text_size <- 8
title_size <- 9
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, face="bold", color="black"), plot.margin=margin(0,-3,0,-3), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", color=NA),  plot.background=element_rect(fill="transparent", color=NA), 
	legend.box.spacing=unit(0, "pt"))

osn_rna <- readRDS("osn_rna2.rds")
osn_rna$cell.subtype_fix <- "mOSN"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 36)] <- "iOSN"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 39)] <- "iOSN"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 40)] <- "Basophil"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 41)] <- "Basophil"
osn_rna$cell.subtype_fix[which(osn_rna$cell.cls == 42)] <- "Basophil"

types <- c("Basophil", "iOSN", "mOSN")
cell_col <- col_list[1:length(types)]
names(cell_col) <- types

osn_nano <- readRDS("osn_nano.rds")
osn_nano$cell.subtype_fix <- "mOSN"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 30)] <- "iOSN"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 33)] <- "iOSN"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 36)] <- "iOSN"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 38)] <- "Basophil"
osn_nano$cell.subtype_fix[which(osn_nano$cell.cls == 39)] <- "Basophil"

cs_info <- data.frame(osn_rna@meta.data)
pa <- wrap_elements(ggplot(cs_info, aes(x="Term", y=nFeature_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="# of Genes", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	panel.background=element_blank()))+tag_thm
pb <- wrap_elements(ggplot(cs_info, aes(x="Term", y=nCount_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="# of counts", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	panel.background=element_blank()))+tag_thm
pc <- wrap_elements(ggplot(cs_info, aes(x="Term", y=cell.mt))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="% of mitochondria", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(size=text_size, color="black"), axis.title=element_text(size=title_size, color="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	panel.background=element_blank()))+tag_thm

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
	theme(legend.position="none", panel.background=element_rect(0, linetype=0), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks.x=element_line(linewidth=0.35, color="black"), axis.ticks.y=element_blank(), 
	axis.text=element_text(size=text_size, color="black"), 
	axis.title=element_text(size=title_size, color="black"), 
	legend.title=element_text(size=title_size, color="black", face="bold"), 
	legend.text=element_text(size=text_size, color="black"), 
	legend.key.size=unit(12, "pt"), #legend.box.spacing = unit(2, "pt"), 
	legend.key=element_blank(), legend.background=element_blank()))+tag_thm

terms <- read.csv("sce_nano_len_sub.csv", r=1)[, 1]
cell_info <- read.csv("sce_nano_info.csv", r=1, h=T)
pe <- ggplot(, aes(x=terms, fill=type, color=type))+geom_density(fill=col_list[1], color=col_list[1], alpha=0.7, linewidth=0.8)+
	scale_y_continuous(expand=c(0, 0))+scale_x_continuous(breaks=seq(0, 3000, 500), limits=c(0, 2000), expand=c(0, 0))+
	labs(title=NULL, x="Tagged reads length (bp)", y="Density\n")+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), axis.title=element_text(size=title_size, color="black"), 
	axis.line=element_line(linewidth=0.35, color="black"), axis.ticks.x=element_line(linewidth=0.35, color="black"), axis.ticks.y=element_blank(), 
	axis.text.y=element_blank(), axis.text.x=element_text(size=text_size, color="black"))+tag_thm
pf <- ggplot(cell_info, aes(x="Term", y=Len))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="Length of reads \nper cell (mean)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	axis.title=element_text(size=title_size, color="black"), axis.text=element_text(size=text_size, color="black"), 
	panel.background=element_blank())+tag_thm
pg <- ggplot(cell_info, aes(x="Term", y=Isoform))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.2)+
	labs(title=NULL, x="# of isoforms", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.line=element_line(linewidth=0.35, color="black"), axis.ticks=element_line(linewidth=0.35, color="black"), 
	axis.title=element_text(size=title_size, color="black"), axis.text=element_text(size=text_size, color="black"), 
	panel.background=element_blank())+tag_thm

osn_rna$cell.subtype_fix <- factor(osn_rna$cell.subtype_fix, levels=types)
rec_da <- data.frame(X=osn_rna@reductions$umap@cell.embeddings[,2], 
	Y=-osn_rna@reductions$umap@cell.embeddings[,1], Type=osn_rna$cell.subtype_fix)
pda <- ggplot(rec_da, aes(x=X, y=Y, color=Type))+geom_point(size=0.1)+
	labs(title="Identified with illumina", x=NULL, y="UMAP2", color=NULL)+
	scale_color_manual(values=cell_col, drop=F)+
	guides(color=guide_legend(nrow=1, override.aes=list(size=3), reverse=T))+
	theme(legend.position="none", axis.line=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=0.35), 
	axis.text.x=element_text(size=text_size, color="black"), 
	axis.text.y=element_text(size=text_size, color="black", angle=90, hjust=0.5, vjust=0.5), 
	axis.title=element_text(size=title_size, color="black"), 
	plot.title=element_text(size=title_size, hjust=0.8, color="black"), 
	legend.box.spacing=unit(0, "pt"), legend.key.size=unit(12, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.text=element_text(size=text_size, color="black"), 
	legend.title=element_text(size=title_size, color="black"), 
	panel.background=element_rect(0, linetype=0))+tag_thm
osn_nano$cell.subtype_fix <- factor(osn_nano$cell.subtype_fix, levels=types)
rec_db <- data.frame(Y=osn_nano@reductions$umap@cell.embeddings[,2], 
	X=osn_nano@reductions$umap@cell.embeddings[,1], Type=osn_nano$cell.subtype_fix)
pdb <- ggplot(rec_db, aes(x=X, y=Y, color=Type))+geom_point(size=0.1)+
	labs(title="Identified with nanopore", x="UMAP1", y="UMAP2", color=NULL)+
	scale_color_manual(values=cell_col, drop=F)+
	guides(color=guide_legend(nrow=1, override.aes=list(size=3), reverse=T))+
	theme(legend.position="bottom", panel.border=element_rect(color="black", fill=NA, linewidth=0.35), 
	axis.text.x=element_text(size=text_size, color="black"), 
	axis.text.y=element_text(size=text_size, color="black", angle=90, hjust=0.5, vjust=0.5), 
	axis.title=element_text(size=title_size, color="black"), 
	plot.title=element_text(size=title_size, hjust=0.8, color="black"), 
	legend.box.spacing=unit(0, "pt"), legend.key.size=unit(12, "pt"), 
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
	labs(title="Consistency of cell types", x="Identified with illumina", y="Identified with nanopore", fill="PCT. (%)")+
	scale_fill_viridis()+scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))+
	theme(panel.background=element_blank(), axis.ticks=element_blank(), legend.key.size=unit(12, "pt"), 
	strip.text=element_blank(), strip.background=element_blank(), 
	axis.text.x=element_text(size=text_size, color="black"), axis.text.y=element_text(size=text_size, color="black", angle=90, hjust=0.5), 
	axis.title=element_text(size=title_size, color="black"), 
	plot.title=element_text(size=title_size, color="black", hjust=0.5), 
	legend.title=element_text(size=title_size, color="black"), 
	legend.text=element_text(size=text_size, color="black"))+tag_thm

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
	theme(axis.text.x=element_blank(), axis.ticks=element_blank(), legend.key.size=unit(12, "pt"), 
	axis.text.y=element_text(size=text_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(colour="#800026", fill="#800026"), 
	strip.text=element_text(size=title_size, colour="white", face="bold"), 
	panel.spacing=unit(5, "pt"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
	axis.line=element_blank(), panel.background=element_rect(fill='gray98')))+tag_thm

cmp_raw_diff <- read.csv("cmp_res_raw_filter.csv", r=1, h=T)
res_l <- data.frame(Group="Differentiation", Diff=cmp_raw_diff$Diff*100)
cmp_vdv_raw <- read.csv("cmp_vdv_filter3.csv", h=T, r=1)
res_l <- rbind(res_l, data.frame(Group="D-V projection", Diff=cmp_vdv_raw$Diff[which(cmp_vdv_raw$Group == "IDG")]))
cmp_vap_raw <- read.csv("cmp_vap_filter3.csv", h=T, r=1)
res_l <- rbind(res_l, data.frame(Group="A-P projection", Diff=cmp_vap_raw$Diff[which(cmp_vap_raw$Group == "IDG")]))
res_l <- res_l[which(!is.na(res_l$Diff)),]
res_l <- res_l[which(res_l$Diff != 100),]
res_l$Group <- factor(res_l$Group, levels=c("Differentiation", "D-V projection", "A-P projection"), 
	labels=c(paste0("Differentiation\n(n=", nrow(cmp_raw_diff), ")"), 
	paste0("D-V projection\n(n=", length(which(cmp_vdv_raw$Group == "IDG")),")"), 
	paste0("A-P projection\n(n=", length(which(cmp_vap_raw$Group == "IDG")), ")")))
pll <- wrap_elements(ggplot(res_l, aes(x=Group, y=Diff, fill=Group, color=Group))+geom_violin()+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white")+scale_y_continuous()+
	scale_x_discrete(drop=F)+scale_fill_manual(values=col_list[c(3,4,5)], drop=F)+
	scale_color_manual(values=col_list[c(3,4,5)], drop=F)+
	labs(title=NULL, x=NULL, y="Difference in proportions")+
	theme(panel.background=element_rect(0, linetype = 0), panel.border=element_rect(fill='transparent', linewidth=0.35, color='black'), 
	axis.title=element_text(size=title_size, colour="black"), axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black", angle=30, hjust=1), 
	axis.ticks=element_line(linewidth=0.35, color="black"), 
	legend.position="none", plot.margin=margin(-5,0,-5,0)))+tag_thm

targets <- c("Nqo1", "Nfix", "Fgf12-24", "Fgf12-99")
groups <- c("D", "V")
data_raw <- read.csv("qpcr_info.csv", h=T)
res <- data.frame()
for (target in targets) for (group in groups)
{
	ids <- which(data_raw$Target == target & data_raw$Group == group)
	res <- rbind(res, data.frame(Target=target, Group=group, Rep=9, 
		RQ.mean=mean(data_raw$RQ[ids]), RQ.sd=sd(data_raw$RQ[ids])))
}
res$RQ.sem <- res$RQ.sd/sqrt(res$Rep)
res$Target <- factor(res$Target, levels=targets)
#t.test(data_raw$RQ[which(data_raw$Target == "Nqo1" & data_raw$Group == "V")], data_raw$RQ[which(data_raw$Target == "Nqo1" & data_raw$Group == "D")])
#t.test(data_raw$RQ[which(data_raw$Target == "Nfix" & data_raw$Group == "V")], data_raw$RQ[which(data_raw$Target == "Nfix" & data_raw$Group == "D")])
#t.test(data_raw$RQ[which(data_raw$Target == "Fgf12-24" & data_raw$Group == "V")], data_raw$RQ[which(data_raw$Target == "Fgf12-24" & data_raw$Group == "D")])
#t.test(data_raw$RQ[which(data_raw$Target == "Fgf12-99" & data_raw$Group == "V")], data_raw$RQ[which(data_raw$Target == "Fgf12-99" & data_raw$Group == "D")])
res$sig <- c("***", "***", "***", "***", "NS", "NS", "***", "***")
res_sub <- res[which(res$Target == "Nqo1" | res$Target == "Nfix"),]
data_sub <- data_raw[which(data_raw$Target == "Nqo1" | data_raw$Target == "Nfix"),]
pia <- wrap_elements(ggplot(res_sub, aes(x=Group, y=RQ.mean))+
	geom_bar(stat="identity", color="black", position=position_dodge(0.7), width=0.6, linewidth=0.3)+
	geom_errorbar(aes(ymin=RQ.mean, ymax=RQ.mean+RQ.sem), size=0.3, width=0.2, position=position_dodge(0.7))+
	geom_barsignif(aes(y=RQ.mean+RQ.sem, signif=sig), control_index=1, vjust=0.04, position=position_dodge(0.7))+
	geom_quasirandom(aes(x=Group, y=RQ), data_sub, pch=21, size=1, fill="white", groupOnX=T)+
	facet_wrap(~Target, nrow=1, scales="free")+
	labs(title=NULL, x=NULL, y="Relative Expression")+
	theme(panel.background=element_rect(0, linetype = 0), panel.border=element_rect(fill='transparent',colour='black'), 
	strip.background=element_blank(), axis.ticks=element_line(colour="black"), 
	plot.title=element_text(size=title_size, colour="black", hjust=0.5), 
	strip.text=element_text(size=title_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), plot.margin=margin(-10,0,0,0)))+tag_thm
res_sub <- res[which(res$Target == "Fgf12-24" | res$Target == "Fgf12-99"),]
data_sub <- data_raw[which(data_raw$Target == "Fgf12-24" | data_raw$Target == "Fgf12-99"),]
res_line <- data.frame(Target=c("Fgf12-24", "Fgf12-99"), 
	a=res_sub$RQ.mean[c(2, 4)] - res_sub$RQ.mean[c(1, 3)], 
	b=res_sub$RQ.mean[c(1, 3)] * 2 - res_sub$RQ.mean[c(2, 4)])
pib <- wrap_elements(ggplot(res_sub, aes(x=Group, y=RQ.mean))+
	geom_bar(stat="identity", color="black", position=position_dodge(0.7), width=0.6, linewidth=0.3)+
	geom_errorbar(aes(ymin=RQ.mean, ymax=RQ.mean+RQ.sem), size=0.3, width=0.2, position=position_dodge(0.7))+
	geom_barsignif(aes(y=RQ.mean+RQ.sem, signif=sig), control_index=1, vjust=0.04, position=position_dodge(0.7))+
	geom_abline(aes(slope=a, intercept=b), res_line, linetype="dashed", color=col_list[3], size=1)+
	geom_quasirandom(aes(x=Group, y=RQ), data_sub, pch=21, size=1, fill="white", groupOnX=T)+
	facet_wrap(~Target, nrow=1, scales="free")+
	labs(title=NULL, x=NULL, y="Relative Expression")+
	theme(panel.background=element_rect(0, linetype = 0), panel.border=element_rect(fill='transparent',colour='black'), 
	strip.background=element_blank(), axis.ticks=element_line(colour="black"), 
	plot.title=element_text(size=title_size, colour="black", hjust=0.5), 
	strip.text=element_text(size=title_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), plot.margin=margin(-10,0,0,0)))+tag_thm
pii <- wrap_plots(list(pia, pib), nrow=1)+tag_thm

pblank <- wrap_elements(ggplot()+geom_blank()+theme(panel.background=element_blank()))+tag_thm
ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pb, pc, pd), nrow=1, widths=c(1,1,1,2))+
	plot_annotation(tag_levels=list(c("A", "B", "C", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pe, pf, pg, pca), nrow=1, widths=c(1.4,0.8,0.8,1.5))+
	plot_annotation(tag_levels=list(c("E", "F", "G", "H")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pdx, pga), widths=c(1,3))+
	plot_annotation(tag_levels=list(c("I", "J"))))+tag_thm, 
	wrap_elements(wrap_plots(list(pii, pll), widths=c(3, 1))+
	plot_annotation(tag_levels=list(c("K", "", "L"))))+tag_thm), 
	ncol=1, heights=c(1,1,1.6,1))+tag_thm, 
	width=210, height=297, dpi=300, units="mm", filename="oe_fl_fig_S03.png", limitsize=F)

ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(pa, pb, pc, pd), nrow=1, widths=c(1,1,1,2))+
	plot_annotation(tag_levels=list(c("A", "B", "C", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pe, pf, pg, pca), nrow=1, widths=c(1.4,0.8,0.8,1.5))+
	plot_annotation(tag_levels=list(c("E", "F", "G", "H")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pdx, pga), widths=c(1,3))+
	plot_annotation(tag_levels=list(c("I", "J"))))+tag_thm, 
	wrap_elements(wrap_plots(list(pii, pll), widths=c(3, 1))+
	plot_annotation(tag_levels=list(c("K", "", "L"))))+tag_thm), 
	ncol=1, heights=c(1,1,1.6,1))+tag_thm, 
	width=210, height=297, dpi=300, units="mm", filename="oe_fl_fig_S03.pdf", limitsize=F)
