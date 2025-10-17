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
options(stringsAsFactors=FALSE)

col_list <- c(c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8"), brewer.pal(12,"Set3")[-c(2, 9)])
text_size <- 13
title_size <- 15
choose_font("Arial")
tag_thm <- theme(plot.tag=element_text(size=title_size, colour="black"), plot.margin=margin(0,-3,2,-3), panel.spacing=unit(0, "pt"), 
	panel.background=element_rect(fill="transparent", colour=NA),  plot.background=element_rect(fill="transparent", colour=NA), 
	legend.box.spacing=unit(0, "pt"))
tag_thm2 <- theme(plot.tag=element_text(size=title_size, colour="black"), plot.margin=margin(-10,-3,2,2), panel.spacing=unit(0, "pt"), 
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

rec_aa <- data.frame(X=osn_rna@reductions$umap@cell.embeddings[,1], 
	Y=osn_rna@reductions$umap@cell.embeddings[,2], Type=osn_rna$cell.subtype_fix)
rec_aa$Type <- factor(rec_aa$Type, levels=rev(types))
paa <- ggplot(rec_aa, aes(x=X, y=Y, color=Type))+geom_point(size=1)+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour=NULL)+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(9,8,7,6,5)], drop=F)+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	annotate("text", x=-4, y=4, label="HBC → mOSN", color="black", size=4, hjust=0, vjust=1)+
	theme(axis.line=element_blank(), 
	panel.border=element_rect(color="black", fill=NA, linewidth=0.8), 
	plot.title=element_text(size=title_size, hjust=0.5, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"),  axis.title=element_text(size=title_size, colour="black"), 
	legend.position=c(0.2, 0.75))+tag_thm

rec_ab <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_iso.tsv", h=T)
rec_ab$type[which(rec_ab$type == "Mature")] <- "mOSN"
rec_ab$type[which(rec_ab$type == "Immature")] <- "iOSN"
rec_ab$type <- factor(rec_ab$type, levels=types)
r = summary(aov(count~type, data=rec_ab))
pab <- ggplot(rec_ab, aes(x=count, color=type)) + stat_ecdf(linewidth=0.8)+
	labs(title=NULL, x="Isoforms per gene", y="ECDF", color=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)))+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	scale_x_continuous(breaks=seq(1, 10, 2), limits=c(0.9, 6.1))+
	#scale_x_continuous(expand=c(0, 0))+scale_y_continuous(expand=c(0, 0))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	legend.title=element_text(size=text_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), axis.line=element_line(colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.position = c(0.76, 0.33))+tag_thm

genes <- unique(rec_ab$gene)
res <- matrix(0, nrow=length(genes), ncol=length(types), dimnames=list(genes, types))
for (i in 1:length(types))
{
	rec_sub <- rec_ab[which(rec_ab$type == types[i]),]
	res[match(rec_sub$gene, genes), i] <- rec_sub$count
}
res <- res[intersect(which(rowMaxs(res) > 1), which(rowSums(res > 0) > 3)),]
rec <- data.frame(Term=rownames(res), Cor=cor(t(res), 1:5))
rec <- rec[order(rec$Cor, decreasing=T),]
rec$Rank <- 1:nrow(rec)
rec$Group <- "N.S."
rec$Group[which(rec$Cor > 0.3)] <- "Increased"
rec$Group[which(rec$Cor < -0.3)] <- "Decreased"
rec$Group <- factor(rec$Group, levels=c("Increased", "Decreased", "N.S."))
rec$Symbol <- gtf$gene_name[match(gsub(":.*", "", rec$Term), gtf$gene_id)]
rec$Anno <- ""
rec$Anno[which(rec$Symbol == "Usp48")] <- "Usp48"
rec_ac <- rec
#terms <- rec_ac$Symbol[which(rec_ac$Group == "Increased")]
#eid <- mapIds(org.Mm.eg.db, keys=terms, column="ENTREZID", keytype="SYMBOL", multiVals="first")
#ego <- enrichGO(gene=eid, keyType="ENTREZID", OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH", readable=T)
#ego <- ego@result[which((ego@result$p.adjust < 0.05 | ego@result$qvalue < 0.05) & ego@result$Count > 1),]
#ego_id <- data.frame(ID=ego$ID, level=0)
#id_filtered <- c()
#for (i in 1:nrow(ego_id))
#{
#	if (!is.na(match(ego_id$ID[i], id_filtered))) next
#	rec <- as.character(unlist(mget(ego_id$ID[i], GOBPPARENTS, ifnotfound=NA)))
#	parents <- rec
#	level <- 0
#	while(length(parents) > 0)
#	{
#		parents_ori <- parents
#		parents <- c()
#		for (p in parents_ori)
#		{
#			if (p == "all" & ego_id$level[i] == 0) ego_id$level[i] <- level
#			if (is.na(p) | p == "all") next
#			rec <- c(rec, p)
#			parents <- c(parents, as.character(unlist(mget(p, GOBPPARENTS, ifnotfound=NA))))
#		}
#		level <- level + 1
#		parents <- unique(parents)
#	}
#	rec <- unique(rec)
#	rec <- which(!is.na(match(ego_id$ID, rec)))
#	if (length(rec) > 0) id_filtered <- c(id_filtered, ego_id$ID[rec])
#}
#ego_id$level[match(unique(id_filtered), ego_id$ID)] <- 0
#ego_id <- ego_id[which(ego_id$level > 2),, drop=F]
#ego <- ego[match(ego_id$ID, ego$ID),, drop=F]
#ego$Level <- ego_id$level
#ego <- ego[order(ego$p.adjust),]
#ego <- ego[order(ego$Count, decreasing=T),]
#write.csv(ego, "num_go.csv")
pac <- ggplot(rec_ac, aes(x=Rank, y=Cor, color=Group)) + geom_point(size=0.8)+
	labs(title=NULL, x="Rank", y="Incrementality", color="")+
	geom_point(data=rec_ac[which(rec_ac$Anno != ""),,drop=F], aes(x=Rank, y=Cor, color=Group), color="red2", size=0.8)+
	geom_text_repel(data=rec_ac[which(rec_ac$Anno != ""),,drop=F], aes(label=Anno), color="black", size=5, 
	segment.size=0.5, direction="y", nudge_y=0.1, nudge_x=20, hjust=0)+
	geom_hline(yintercept=c(0.3, -0.3), color="gray30", linetype="dashed", linewidth=0.6)+
	scale_y_continuous(breaks=c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9))+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	scale_color_manual(values=c(col_list[3], col_list[4], "gray75"))+
	#scale_x_continuous(expand=c(0, 0))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position = c(0.37, 0.22))+tag_thm

gene <- "Usp48"
cell_num <- min(table(osn_rna$cell.subtype_fix))
col_types <- brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)]
names(col_types) <- types
ids <- which(osn_sct_raw[["features"]]$Symbol == gene & osn_sct_raw[["features"]]$Type == "protein_coding")
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
	labs(title=NULL, x=NULL, y="Expression")+
	#geom_flow(alpha=1)+geom_stratum(alpha=1, color="white", linewidth=0)+
	geom_flow(alpha=0.5)+geom_stratum(alpha=1, color="white")+guides(fill="none")+
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
#info_ranges <- paste0(info_exons$seqnames[1], strand, ":", min(info_exons$start), "-", max(info_exons$start))
info_ranges <- paste0(info_exons$seqnames[1], strand, ":", min(info_exons$start), "-", max(info_exons$start))
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
	scale_y_discrete(expand=expansion(mult=c(0.1,0.23)))+
	scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, labels=gene)+
	geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), x=0, hjust=0, vjust=0, nudge_y=0.3, size=4)+
	scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
	theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
	axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
	axis.text.x=element_text(size=text_size, colour="black"), legend.position="none")+tag_thm
rec_ppp <- data.frame(Type=c("HBC", "GBC", "INP", "iOSN", "   mOSN"), Info="Line")
rec_ppp$Type <- factor(rec_ppp$Type, levels=c("HBC", "GBC", "INP", "iOSN", "   mOSN"))
ppp <- ggplot(rec_ppp, aes(y=1, x=Type, fill=Type))+geom_bar(stat="identity", width=0.32)+
	theme_minimal()+labs(title=NULL, y=NULL, x=info_ranges, fill=NULL)+
	scale_fill_manual(values=as.character(col_types[c("HBC", "GBC", "INP", "iOSN", "mOSN")]))+
	scale_x_discrete(expand=c(0, 0))+scale_y_continuous(expand=c(0, 0))+
	theme(legend.position="none", panel.background=element_blank(), axis.title.y=element_blank(), 
	axis.title.x=element_text(size=title_size, colour="black", hjust=-0.4), 
	axis.text.y=element_blank(), axis.text.x=element_text(size=text_size, colour="black"), 
	axis.ticks=element_blank())+tag_thm
pad <- wrap_elements(wrap_plots(A=ppa, B=ppb, C=ppp, design=c("AB\nCB"), heights=c(40, 1))+
	plot_annotation(NULL, theme=theme(plot.title=element_text(size=title_size, colour="black", hjust=0.5), 
	plot.margin=margin(0,0,-32,0), panel.spacing=unit(0, "pt"))))+tag_thm

ics <- which(osn_sct_raw[["features"]]$Type == "NOVEL_IC")
ncs <- which(osn_sct_raw[["features"]]$Type == "NOVEL_NC")
#terms <- unique(osn_sct_raw[["features"]]$Gene[c(ics, ncs)])
#rec <- data.frame()
#for (x in terms)
#{
#	ids <- which(osn_sct_raw[["features"]]$Gene == x)
#	cs <- colSums(osn_sct_raw[["trans"]][ids,, drop=F])
#	rs <- rowSums(osn_sct_raw[["trans"]][ids,, drop=F])
#	ids <- grep("NOVEL", names(rs))
#	rec <- rbind(rec, data.frame(Term=names(rs)[ids], Count=rs[ids], Rate=rs[ids]*100/sum(rs), Cell=length(which(cs > 0))))
#}
#terms <- unique(osn_sct_raw[["features"]]$Symbol[match(rec$Term[which(rec$Count > 20 & rec$Rate > 40)], osn_sct_raw[["features"]]$ID)])
#eid <- mapIds(org.Mm.eg.db, keys=terms, column="ENTREZID", keytype="SYMBOL", multiVals="first")
#ego <- enrichGO(gene=eid, keyType="ENTREZID", OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH", readable=T)
#ego <- ego@result[which((ego@result$p.adjust < 0.05 | ego@result$qvalue < 0.05) & ego@result$Count > 1),]
#ego_id <- data.frame(ID=ego$ID, level=0)
#id_filtered <- c()
#for (i in 1:nrow(ego_id))
#{
#	if (!is.na(match(ego_id$ID[i], id_filtered))) next
#	rec <- as.character(unlist(mget(ego_id$ID[i], GOBPPARENTS, ifnotfound=NA)))
#	parents <- rec
#	level <- 0
#	while(length(parents) > 0)
#	{
#		parents_ori <- parents
#		parents <- c()
#		for (p in parents_ori)
#		{
#			if (p == "all" & ego_id$level[i] == 0) ego_id$level[i] <- level
#			if (is.na(p) | p == "all") next
#			rec <- c(rec, p)
#			parents <- c(parents, as.character(unlist(mget(p, GOBPPARENTS, ifnotfound=NA))))
#		}
#		level <- level + 1
#		parents <- unique(parents)
#	}
#	rec <- unique(rec)
#	rec <- which(!is.na(match(ego_id$ID, rec)))
#	if (length(rec) > 0) id_filtered <- c(id_filtered, ego_id$ID[rec])
#}
#ego_id$level[match(unique(id_filtered), ego_id$ID)] <- 0
#ego_id <- ego_id[which(ego_id$level > 2),, drop=F]
#ego <- ego[match(ego_id$ID, ego$ID),, drop=F]
#ego$Level <- ego_id$level
#ego <- ego[order(ego$pvalue),]
#ego <- ego[order(ego$Count, decreasing=T),]
#write.csv(ego, "test_go4.csv")
rec_baa <- data.frame(Type=c("IC", "NC"), Number=c(length(ics), length(ncs)))
rec_baa$Type <- factor(rec_baa$Type, levels=c("NC", "IC"))
pbaa <- ggplot(rec_baa, aes(x=4, y=Number, fill=Type))+geom_col()+
	labs(title="Total", x=NULL, y=NULL)+
	geom_text(aes(label=Number), size=4, position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+scale_x_continuous(limits=c(2.5, 4.5), expand=c(0, 0)) +
	scale_fill_manual(values=col_list[c(6,2)])+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	plot.background=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), 
	axis.text=element_blank(), axis.title=element_blank(), legend.title=element_text(colour="black", size=title_size), 
	legend.text=element_text(colour="black", size=text_size), legend.position="none")+tag_thm
rec_bab <- data.frame()
for (type in types)
{
	cs <- which(osn_rna$cell.subtype_fix == type)
	rec_bab <- rbind(rec_bab, data.frame(Type=type, Number=length(which(rowSums(osn_sct_raw[["trans"]][c(ics, ncs), cs] > 1) > 1))))
}
rec_bab$Type <- factor(rec_bab$Type, levels=types)
rec_bab$Group <- "In catalog"
rec_bac <- data.frame()
for (type in types)
{
	cs <- which(osn_rna$cell.subtype_fix == type)
	rec_bac <- rbind(rec_bac, data.frame(Type=type, Number=length(which(rowSums(osn_sct_raw[["trans"]][ncs, cs] > 1) > 1))))
}
rec_bac$Type <- factor(rec_bac$Type, levels=types)
rec_bac$Group <- "Novel junctions"
rec_bab <- rbind(rec_bab, rec_bac)
pbab <- ggplot(rec_bab, aes(x=Type, y=Number, fill=Group))+geom_bar(stat="identity", width=0.6)+
	labs(title=NULL, x=NULL, y="Novel isoforms", fill=NULL)+
	scale_fill_manual(values=col_list[c(2, 6)], drop=F)+
	scale_y_continuous(expand=c(0, 0))+guides(fill=guide_legend(nrow=1))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position=c(0.5, -0.15))+tag_thm
pba <- pbab+inset_element(pbaa, 0.02, 0.2, 0.55, 1)+tag_thm

#terms <- unique(osn_sct_raw[["features"]]$Symbol[match(rec$Term[which(rec$Count > 20 & rec$Rate > 40)], osn_sct_raw[["features"]]$ID)])
#eid <- mapIds(org.Mm.eg.db, keys=terms, column="ENTREZID", keytype="SYMBOL", multiVals="first")
#ego <- enrichGO(gene=eid, keyType="ENTREZID", OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH", readable=T)
#ego <- ego@result[which((ego@result$p.adjust < 0.05 | ego@result$qvalue < 0.05) & ego@result$Count > 1),]
#ego_id <- data.frame(ID=ego$ID, level=0)
#id_filtered <- c()
#for (i in 1:nrow(ego_id))
#{
#	if (!is.na(match(ego_id$ID[i], id_filtered))) next
#	rec <- as.character(unlist(mget(ego_id$ID[i], GOBPPARENTS, ifnotfound=NA)))
#	parents <- rec
#	level <- 0
#	while(length(parents) > 0)
#	{
#		parents_ori <- parents
#		parents <- c()
#		for (p in parents_ori)
#		{
#			if (p == "all" & ego_id$level[i] == 0) ego_id$level[i] <- level
#			if (is.na(p) | p == "all") next
#			rec <- c(rec, p)
#			parents <- c(parents, as.character(unlist(mget(p, GOBPPARENTS, ifnotfound=NA))))
#		}
#		level <- level + 1
#		parents <- unique(parents)
#	}
#	rec <- unique(rec)
#	rec <- which(!is.na(match(ego_id$ID, rec)))
#	if (length(rec) > 0) id_filtered <- c(id_filtered, ego_id$ID[rec])
#}
#ego_id$level[match(unique(id_filtered), ego_id$ID)] <- 0
#ego_id <- ego_id[which(ego_id$level > 2),, drop=F]
#ego <- ego[match(ego_id$ID, ego$ID),, drop=F]
#ego$Level <- ego_id$level
#ego <- ego[order(ego$pvalue),]
#ego <- ego[order(ego$Count, decreasing=T),]
#write.csv(ego, "test_go.csv")

reads_info <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_len.tsv", h=T)
reads_info$type <- factor(reads_info$type, levels=types)
reads_info_bb <- reads_info[which(reads_info$len > 10),]
r = summary(aov(len~type, data=reads_info_bb))
pbb <- ggplot(reads_info_bb, aes(x=len, color=type)) + stat_ecdf(linewidth=0.8)+
	labs(title=NULL, x="Length of isoforms", y="ECDF", color=NULL)+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	scale_x_continuous(limits=c(0, 3000), expand=c(0, 0))+
	annotate("text", x=100, y=1, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	color="black", size=4, hjust=0, vjust=1)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position = c(0.82, 0.3))+tag_thm

reads_info <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_len.tsv", h=T)
reads_info$type <- factor(reads_info$type, levels=types)
reads_info_bc <- reads_info[which(reads_info$utr3 > 10 & reads_info$body > 10),]
r = summary(aov(len~type, data=reads_info_bc))
pbc <- ggplot(reads_info_bc, aes(x=utr3, color=type)) + stat_ecdf(linewidth=0.8)+
	labs(title=NULL, x="Length of 3'UTRs", y="ECDF", color=NULL)+
	scale_x_continuous(limits=c(0, 1500), expand=c(0, 0))+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,8,9)])+
	annotate("text", x=50, y=1, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	color="black", size=4, hjust=0, vjust=1)+
	theme(plot.title=element_text(size=title_size, hjust=0.5), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	legend.key=element_blank(), legend.background=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position = c(0.82, 0.3))+tag_thm

sce_trans <- CreateSeuratObject(scg[["matrix"]], project="OEnano")
sce_trans$cell.subtype_fix <- as.character(osn_rna$cell.subtype_fix[match(colnames(sce_trans), colnames(osn_rna))])
sce_trans <- SCTransform(sce_trans, method="glmGamPoi")
reads_info <- read.delim("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_len.tsv", h=T)
reads_info$type[which(reads_info$type == "Mature")] <- "mOSN"
reads_info$type[which(reads_info$type == "Immature")] <- "iOSN"
reads_info$type <- factor(reads_info$type, levels=types)
rec <- data.frame()
for (type in types[-4])
{
	terms <- FindMarkers(sce_trans, ident.1=type, group.by="cell.subtype_fix")
	terms <- rownames(terms)[which(terms$p_val_adj < 0.01 & terms$avg_log2FC > 2)]
	rec <- rbind(rec, reads_info[match(paste0(type, ",", terms), paste0(reads_info$type, ",", gsub("_", "-", reads_info$trans))),])
}
rec_bd <- rec[which(rec$utr3 > 10 & rec$body > 10 & rec$utr3 < 2000),]
rec_bd$type <- factor(rec_bd$type, levels=types[-4], labels=c(types[1:3], "OSN"))
r = summary(aov(utr3~type, data=rec_bd))
pbd <- ggplot(rec_bd, aes(x=type, y=utr3, fill=type, color=type))+geom_violin()+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white")+scale_y_continuous()+
	scale_x_discrete(drop=F)+scale_fill_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,9)], drop=F)+
	scale_color_manual(values=brewer.pal(9,"YlGnBu")[c(5,6,7,9)], drop=F)+
	annotate("text", x=0.5, y=1500, label=paste0("F=", round(r[[1]][["F value"]][1], 2), "\nPval=", round(r[[1]][["Pr(>F)"]][1], 2)), 
	color="black", size=4, hjust=0, vjust=1)+
	labs(title=NULL, x="3'UTR length of\nmarker isoforms", y="Length of 3'UTR")+
	theme(axis.title=element_text(size=title_size, colour="black"), axis.text=element_text(size=text_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position="none")+tag_thm

sce_trans <- CreateSeuratObject(osn_sct[["genes"]], project="OEnano")
sce_trans$cell.subtype_fix <- as.character(osn_rna$cell.subtype_fix[match(colnames(sce_trans), colnames(osn_rna))])
sce_trans <- SCTransform(sce_trans, method="glmGamPoi", variable.features.n=nrow(sce_trans))
mat_utr3 <- round(scg[["utr3"]]*scg[["body"]]/(scg[["body"]]+0.1))
mat_utr3 <- scg[["utr3"]][, match(colnames(osn_rna), colnames(scg[["utr3"]]))]
mat_utr3 <- mat_utr3[which(rowSums(mat_utr3 > 10) > 10),]
mat_utr5 <- round(scg[["utr5"]]*scg[["body"]]/(scg[["body"]]+0.1))
mat_utr5 <- scg[["utr5"]][, match(colnames(osn_rna), colnames(scg[["utr5"]]))]
mat_utr5 <- mat_utr5[which(rowSums(mat_utr5 > 10) > 10),]
rec <- data.frame(Term=rownames(mat_utr3), Gene=scg[["features"]]$Gene[match(rownames(mat_utr3), scg[["features"]]$Name)], 
	Trans=scg[["features"]]$Trans[match(rownames(mat_utr3), scg[["features"]]$Name)], 
	Symbol=scg[["features"]]$Symbol[match(rownames(mat_utr3), scg[["features"]]$Name)], Pval.e=1, FC.e=0, Pval=1, FC=0, Len=0, Diff=0)
rec <- rec[which(!is.na(match(rec$Gene, rownames(sce_trans)))),]
term_ids <- match(rec$Term, rownames(mat_utr3))
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
	tta <- mat_utr3[term_ids[i], ids_a]
	ttb <- mat_utr3[term_ids[i], ids_b]
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
	#ggsave(plot=ptt, width=5, height=4, dpi=200, paste0("utr3_", rec_sub$Symbol, "_dist.png"), limitsize=F)
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
deg_info$Anno[which(deg_info$Symbol == "Calm1")] <- "Calm1"
#deg_info$Anno <- deg_info$Symbol
#deg_info$Anno[which(deg_info$Group == "N.S.")] <- ""
#deg_info$Anno[which(deg_info$Anno != "")[-c(1:10)]] <- ""
deg_info$Group <- factor(deg_info$Group, levels=c("Increased", "Decreased", "N.S."))
deg_info <- deg_info[which(abs(deg_info$Diff) > 10),]
pca <- ggplot(deg_info, aes(x=X, y=Y, color=Group))+geom_point(size=2)+
	labs(title=NULL, x="Difference of 3'UTR length", y="P.adj (-log10)", color=NULL)+
	scale_color_manual(values=c("#FC8D62", "#8DA0CB", "gray80"))+scale_y_continuous(expand=c(0, 0))+
	guides(colour=guide_legend(override.aes=list(size=4)))+
	#geom_text(aes(y=Y+1, label=Anno), size=3, col="black", vjust=0)+
	geom_text_repel(data=deg_info[which(deg_info$Anno != ""),,drop=F], aes(label=Anno), color="black", size=5, 
	segment.size=0.5, direction="y", , nudge_y=2, nudge_x=0.05, hjust=0)+
	geom_hline(yintercept=-log10(0.05), colour="gray30", linetype="dashed", linewidth=0.5)+
	geom_vline(xintercept=c(-30, 30), colour="gray30", linetype="dashed", linewidth=0.5)+
	theme(plot.title=element_text(size=title_size, hjust=0.5, colour="black"), panel.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), legend.text=element_text(size=text_size, colour="black"), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(colour="black"), legend.position=c(0.8, 0.83))+tag_thm

#terms <- rec$Term[which(rec$Padj < 0.05 & rec$FC > 0)]
term <- deg_info$Term[which(deg_info$Symbol == "Calm1")]
trans <- scg[["features"]]$Trans[match(term, scg[["features"]]$Name)]
gene <- scg[["features"]]$Symbol[match(term, scg[["features"]]$Name)]
term_data <- data.frame(Cell=colnames(mat_utr3), Type=osn_rna$cell.subtype_fix, Length=mat_utr3[term,])
#term_data <- term_data[which(term_data$Length > max(mean(term_data$Length)*2/3, mean(term_data$Length)-100)),]
term_data <- term_data[which(term_data$Length > 0),]
term_data$Group <- "Early"
term_data$Group[which(term_data$Type == "mOSN" | term_data$Type == "iOSN")] <- "Later"
term_data$Type <- factor(term_data$Type, levels=types)
term_data$Group <- factor(term_data$Group, levels=c("Early", "Later"))
term_data$LG <- "M"
term_data$LG[which(term_data$Length < 250)] <- "S"
term_data$LG[which(term_data$Length > 600)] <- "L"

rec <- data.frame()
for (type in c("Early", "Later")) rec <- rbind(rec, data.frame(Type=type, 
	S=length(which(term_data$Length < 250 & term_data$Group == type)), 
	M=length(which(term_data$Length >= 250 & term_data$Length <= 600 & term_data$Group == type)), 
	L=length(which(term_data$Length > 600 & term_data$Group == type))))
rec$RS <- rec$S*100/(rec$S+rec$M+rec$L)
rec$RM <- rec$M*100/(rec$S+rec$M+rec$L)
rec$RL <- rec$L*100/(rec$S+rec$M+rec$L)
rec_cba <- rbind(data.frame(Type=rec$Type, Group="RS", Count=rec$RS), 
	data.frame(Type=rec$Type, Group="RM", Count=rec$RM), 
	data.frame(Type=rec$Type, Group="RL", Count=rec$RL))
rec_cba$Type <- factor(rec_cba$Type, levels=c("Early", "Later"))
rec_cba$Group <- factor(rec_cba$Group, levels=c("RL", "RM", "RS"))
pcba <- ggplot(rec_cba, aes(x=Type, y=Count, fill=Group))+geom_bar(stat="identity", width=0.6)+
	labs(title=NULL, x="Calm1", y="Percentage (%)", fill=NULL)+
	scale_fill_manual(values=col_list[c(5,4,3)], drop=F)+
	scale_y_continuous(breaks=c(25, 50, 75), expand=c(0, 0))+guides(fill=guide_legend(nrow=1))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	legend.key=element_blank(), legend.background=element_blank(), legend.position="none", axis.ticks.x=element_blank(), 
	axis.text=element_text(size=text_size, colour="black"), axis.title.y=element_text(size=title_size, colour="black"), 
	axis.title.x=element_text(size=title_size, colour="black", hjust=2), 
	axis.line=element_line(colour="black"))+tag_thm
dens <- data.frame(density(term_data$Length)[c("x", "y")])
dens <- dens[which(dens$x > 0 & dens$x < 1100),]
dens$g <- "M"
dens$g[which(dens$x < 250)] <- "S" 
dens$g[which(dens$x > 600)] <- "L" 
#ss <- which(dens$g == "M")
#ss <- dens[c(ss[1], ss[length(ss)]),]
#ss$g <- c("S", "L")
#dens <- rbind(dens, ss)
dens$g <- factor(dens$g, levels=c("S", "M", "L"), labels=c("Short", "Medium", "Long")) 
pcbb <- ggplot(dens, aes(x=x, y=y, color=g, fill=g))+geom_area(alpha=1)+
	labs(title=NULL, x="Length of 3'UTR", y=NULL, color=NULL, fill=NULL)+
	scale_color_manual(values=col_list[c(3,4,5)], drop=F)+
	scale_fill_manual(values=col_list[c(3,4,5)], drop=F)+
	guides(color=guide_legend(reverse=T), fill=guide_legend(reverse=T))+
	scale_x_continuous(breaks=c(0, 250, 500, 750, 1000))+
	scale_y_continuous(breaks=c((max(dens$y)+min(dens$y))/2), labels=c("Distribution"), expand=c(0, 0))+coord_flip()+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	legend.key=element_blank(), legend.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), legend.text=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black", angle=90, hjust=0.5, vjust=0.5), 
	axis.title=element_text(size=title_size, colour="black"), 
	axis.ticks.y=element_line(colour="black"), axis.ticks.x=element_blank(), 
	axis.line.y=element_line(colour="black"), axis.line.x=element_blank(), 
	legend.position=c(0.81, 0.83))+tag_thm
pcb <- wrap_plots(list(pcba, pcbb), nrow=1)+tag_thm

terms <- names(which(table(osn_sct_raw[["features"]]$Gene) > 1))
cells <- lapply(types, function(x) which(osn_rna$cell.subtype_fix == x))
cl <- makeCluster(16, type="FORK")
results <- parLapply(cl, terms, function(x) {
	ids <- which(osn_sct_raw[["features"]]$Gene == x)
	rec_sub <- data.frame(Gene=x, Type=1:length(types), RNA=0, Num=0, Count=0, Ctr=0, Pval=1, Cor=0, Diff=0, Trans="")
	for (j in 1:nrow(rec_sub))
	{
		aid <- j
		bid <- j + 1
		if (j == 5) bid <- 1
		cell_a <- cells[[aid]]
		cell_b <- cells[[bid]]
		cell_num <- min(length(cell_a), length(cell_b))
		mat <- rbind(a=round(rowSums(osn_sct_raw[["trans"]][ids, cell_a])*cell_num/length(cell_a)), b=round(rowSums(osn_sct_raw[["trans"]][ids, cell_b])*cell_num/length(cell_b)))
		cs <- colSums(mat)
		ts <- which(cs > 2)
		if (length(ts) < 2) next
		mat <- mat[, ts]
		cs <- cs[ts]
		rs <- rowSums(mat)
		if (min(rs) < 2) next
		rec_sub$RNA[j] <- sum(osn_sct_raw[["trans"]][ids, c(cell_a, cell_b)])/(length(cell_a)+length(cell_b))
		rec_sub$Num[j] <- ncol(mat)
		rec_sub$Count[j] <- min(rowMaxs(mat))
		rec_sub$Ctr[j] <- min(rs)*100/max(rs)
		mat_df <- mat[1,]/sum(mat[1,])-mat[2,]/sum(mat[2,])
		mat_sub <- mat[, c(which.max(mat_df), which.min(mat_df))]
		rec_sub$Trans[j] <- paste(colnames(mat_sub), collapse=";")
		rec_sub$Pval[j] <- fisher.test(mat_sub)$p.value
		rec_sub$Cor[j] <- cor(mat[1,], mat[2,])
		rec_sub$Diff[j] <- max(abs(mat[1,]/sum(mat[1,])-mat[2,]/sum(mat[2,])))
	}
	return(rec_sub)
})
stopCluster(cl)
rec <- do.call("rbind", results)
write.csv(rec, "cmp_res_raw.csv")
cmp_raw <- read.csv("cmp_res_raw.csv", h=T, r=1)
cmp_raw <- cmp_raw[which(cmp_raw$Type != 5 & cmp_raw$Gene != "-" & cmp_raw$Num > 1 & cmp_raw$Count > 5),]
cmp_raw$Cor[which(is.na(cmp_raw$Cor))] <- 0
cmp_raw$p.adj <- p.adjust(cmp_raw$Pval, method="BH")
#cmp_raw$p.adj <- cmp_raw$Pval
#write.csv(cmp_raw, "cmp_res_filter.csv")
#for (j in 1:5) cmp_raw$p.adj[which(cmp_raw$Type == j)] <- p.adjust(cmp_raw$Pval[which(cmp_raw$Type == j)], "BH")
cmp_raw <- cmp_raw[order(cmp_raw$p.adj),]
cmp_raw$Symbol <- gtf$gene_name[match(cmp_raw$Gene, gtf$gene_id)]
cmp_raw$Y <- -log10(cmp_raw$p.adj)
#cmp_raw$Y[which(cmp_raw$Y > 10)] <- 10
cmp_raw$X <- cmp_raw$Diff
#cmp_raw$X[which(cmp_raw$X > 0.8)] <- 0.78
cmp_raw$Group <- cmp_raw$Type
cmp_raw$Group[which(cmp_raw$p.adj >= 0.2 | cmp_raw$Diff < 0.1)] <- 5
cmp_raw$Group <- factor(cmp_raw$Group, levels=1:5, labels=c("HBC>GBC", "GBC>INP", "INP>iOSN", "iOSN>mOSN", "N.S."))
cmp_raw$Val <- (cmp_raw$Y/max(cmp_raw$Y))^2 + (cmp_raw$X/max(cmp_raw$X))^2
cmp_raw <- cmp_raw[order(cmp_raw$Val, decreasing=T),]
cmp_raw$Anno <- ""
cmp_raw$Anno[which(cmp_raw$Symbol == "Cxadr" & cmp_raw$Pval < 0.05)] <- "Cxadr"
#cmp_raw$Anno[which(cmp_raw$Group != "N.S.")[1:10]] <- cmp_raw$Gene[which(cmp_raw$Group != "N.S.")[1:10]]
write.csv(cmp_raw[which(cmp_raw$Group != "N.S."),], "cmp_res_raw_filter.csv")
rec_da <- cmp_raw
pda <- ggplot(rec_da, aes(x=X, y=Y, color=Group, fill=Group))+geom_point(shape=21, size=3)+
	labs(title=NULL, x="Difference of proportions", y="P.adj (-log10)", fill=NULL, color=NULL)+
	geom_text_repel(data=rec_da[which(rec_da$Anno != ""),,drop=F], aes(label=Anno), color="black", size=5, 
	segment.size=0.5, direction="y", nudge_y=-0.1, nudge_x=0.05, hjust=0)+
	scale_x_continuous(breaks=seq(0.2, 1, 0.2), limits=c(0, max(rec_da$X)+0.1), expand=c(0, 0))+
	scale_y_continuous(limits=c(0, 11), breaks=seq(0, 10, 2), expand=c(0, 0))+
	scale_fill_manual(values=c(colorRampPalette(brewer.pal(9,"YlOrRd")[c(2,5,9)])(4), "grey60"))+
	scale_color_manual(values=c("black", "black", "black", "black", "grey60"))+
	guides(color=guide_legend(override.aes=list(size=4)), fill=guide_legend(override.aes=list(size=4)))+
	geom_vline(xintercept=0.1, color="gray30", linetype="dashed", linewidth=0.6)+
	geom_hline(yintercept=-log10(0.2), color="gray30", linetype="dashed", linewidth=0.6)+
	#geom_text(aes(x=X, y=Y+1, label=Anno), size=0.5*text_size, vjust=0,position=position_jitter(height=1), color="black")+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_rect(0, linetype=0), 
	axis.line=element_line(linetype=1, colour='black'), 
	axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), legend.key.size=unit(15, "pt"), 
	legend.key=element_blank(), legend.background=element_blank(), legend.position=c(0.8, 0.8))+tag_thm

rec_db <- data.frame(table(cmp_raw$Type[which(cmp_raw$Group != "N.S.")]))
rec_db$Var1 <- factor(rec_db$Var1, levels=1:4, 
	labels=c("HBC>GBC", "GBC>INP", "INP>iOSN", "iOSN>mOSN"))
pdb <- ggplot(rec_db, aes(x=Var1, y=Freq, fill=Var1))+geom_bar(stat="identity", width=0.7)+
	labs(title=NULL, x="Group", y="# IDGs", fill=NULL)+
	scale_y_continuous(limits=c(0, max(rec_db$Freq)+4), expand=c(0, 0))+
	geom_text(aes(x=Var1, y=Freq+0.5, label=Freq), size=5, vjust=0)+
	#scale_fill_manual(values=brewer.pal(9,"YlOrRd")[c(3,5,7,9)])+
	scale_fill_manual(values=colorRampPalette(brewer.pal(9,"YlOrRd")[c(2,5,9)])(4))+
	theme(plot.title=element_text(size=title_size, hjust=0.5), panel.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_line(color="black"), 
	axis.text.x=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	axis.title=element_text(size=title_size, colour="black"), 
	axis.line=element_line(color="black"), legend.position=c(0.3, 0.83))+tag_thm

genes <- c("Eef1d", "Cxadr")
plds <- lapply(genes, function (gene) {
	ids <- which(osn_sct_raw[["features"]]$Symbol == gene & osn_sct_raw[["features"]]$Type == "protein_coding")
	rec_tt <- matrix(0, nrow=length(ids), ncol=5, dimnames=list(osn_sct_raw[["features"]]$ID[ids], types))
	for (i in 1:5)
	{
		#cells <- sample(which(osn_rna$cell.subtype_fix == types[i]), cell_num)
		cells <- which(osn_rna$cell.subtype_fix == types[i])
		rec_tt[, i] <- rowSums(osn_sct_raw[["trans"]][ids, cells])*cell_num/length(cells)
	}
	rec_tt <- rec_tt[which(rowMaxs(rec_tt) > 1),]
	isoforms <- rownames(rec_tt)
	#for (i in 1:5) rec_tt[, i] <- rec_tt[, i]*100/sum(rec_tt[, i])
	rec_ttt <- data.frame()
	for (i in which(colSums(rec_tt) > 0)) rec_ttt <- rbind(rec_ttt, data.frame(Type=types[i], 
		Isoform=rownames(rec_tt), Count=rec_tt[, i]))
	rec_ttt$Type <- factor(rec_ttt$Type, levels=types)
	rec_ttt$Isoform <- factor(rec_ttt$Isoform, levels=rev(rownames(rec_tt)))
	pcca <- ggplot(rec_ttt, aes(x=Type, y=Count, stratum=Isoform, alluvium=Isoform, fill=Isoform))+
		labs(title=NULL, x=NULL, y="Expression")+
		#geom_flow(alpha=1)+geom_stratum(alpha=1, color="white", linewidth=0)+
		geom_flow(alpha=0.5)+geom_stratum(alpha=1, color="white")+guides(fill="none")+
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
	info_ranges <- paste0(info_exons$seqnames[1], strand, ":", min(info_exons$start), "-", max(info_exons$start))
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
	junction_info <- read.csv(paste0("/mnt/md0/oe_full_length/output/OEfulllength/res_info_junction_", gene, "_raw.csv"), h=T)
	density_info <- read.csv(paste0("/mnt/md0/oe_full_length/output/OEfulllength/res_info_density_", gene, "_raw.csv"), h=T)
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
	if (density_regions$type[1] == "intron") density_regions <- density_regions[2:nrow(density_regions),]
	if (density_regions$type[nrow(density_regions)] == "intron") density_regions <- density_regions[1:nrow(density_regions)-1,]
	density_regions$il <- density_regions$ie - density_regions$is
	density_regions$sc <- (density_regions$end - density_regions$start) / density_regions$il
	density_info <- density_info[1:(nrow(density_info)-1),]
	gene_list <- gene
	junc_list <- osn_sct_raw[["features"]]$Gene[match(gene, osn_sct_raw[["features"]]$Symbol)]
	sample_list <- unique(density_info$s)
	#xlim <- c(min(density_info$x[which(density_info$y > max(density_info$y)/10)]), max(density_regions$end+1))
	plsc <- lapply(sample_list, function(id)
	{
		density <- data.frame()
		for (t in isoforms)
		{
			density_sub <- density_info[which(density_info$g == gene & density_info$s == id & density_info$t == t),]
			for (i in 1:nrow(density_regions)) for (j in 1:density_regions$il[i]) density <- rbind(density, data.frame(x=j+density_regions$is[i]-1, 
				y=mean(density_sub$y[which(density_sub$x >= round((j-1)*density_regions$sc[i])+density_regions$start[i] & 
				density_sub$x <= round(j*density_regions$sc[i])+density_regions$start[i])]), t=t))
		}
		density$t <- factor(density$t, levels=isoforms)
		gp <- ggplot(density, aes(x, y, color=t, fill=t))+geom_area(position="identity", alpha=0.5, linewidth=0.1)+
			scale_fill_manual(values=col_list[c(2:1)+5])+scale_color_manual(values=col_list[c(2:1)+5])+
			scale_y_continuous(breaks=c((max(density$y)+min(density$y))/2), labels=c(id))+scale_x_continuous(expand=c(0.01, 0.01))+
			theme(panel.background=element_blank(), axis.line=element_blank(), 
			axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), plot.margin=margin(b=4), legend.position="none")
		return(gp)
	})
	pccc <- ggplot(info_exons, aes(xstart=start, xend=end, y=transcript_name))+
		geom_intron(data=info_introns, aes(color=transcript_name, strand=strand), linewidth=0.5, arrow.min.intron.length=200)+
		geom_range(aes(fill=transcript_name, color=transcript_name), linewidth=0, height=0.3)+
		geom_range(data=info_utrs, mapping=aes(xstart=start, xend=end, y=transcript_name, fill=transcript_name, 
		color=transcript_name), linewidth=0.5, height=0.1)+
		labs(title=NULL, x=NULL, y=NULL)+scale_y_discrete(expand=expansion(mult=c(0.5,1)))+
		scale_x_continuous(breaks=(max(c(info_exons$end, info_utrs$end))+min(c(info_exons$start, info_utrs$start)))/2, 
		labels=gene, expand=c(0.01, 0.01))+
		geom_text(aes(y=transcript_name, label=transcript_name, color=transcript_name), 
		x=min(density_regions$is), hjust=0, vjust=0, nudge_y=0.3, size=4)+
		scale_fill_manual(values=info_col, drop=F)+scale_color_manual(values=info_col, drop=F)+
		theme(panel.border=element_rect(color="black", fill=NA, linewidth=0.8), panel.background=element_blank(), 
		axis.line=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), 
		axis.text.x=element_text(size=text_size, colour="black"), legend.position="none")+tag_thm
	rec_ppp <- data.frame(Type=unique(rec_ttt$Type), Info="Line")
	rec_ppp$Type <- factor(rec_ppp$Type, levels=unique(rec_ttt$Type))
	pccb <- ggplot(rec_ppp, aes(y=1, x=Type, fill=Type))+geom_bar(stat="identity", width=0.32)+
		theme_minimal()+labs(title=NULL, y=NULL, x=info_ranges, fill=NULL)+
		scale_fill_manual(values=col_types[unique(rec_ttt$Type)])+
		scale_x_discrete(expand=c(0, 0))+scale_y_continuous(expand=c(0, 0))+
		theme(legend.position="none", panel.background=element_blank(), axis.title.y=element_blank(), 
		axis.title.x=element_text(size=title_size, color="black", hjust=-0.6), 
		axis.text.y=element_blank(), axis.text.x=element_text(size=text_size, color="black"), 
		axis.ticks=element_blank())+tag_thm
	pccd <- ggplot(rec_ppp, aes(x=1, y=Type, fill=Type))+geom_tile()+
		theme_minimal()+labs(title=NULL, y=NULL, x=info_ranges, fill=NULL)+
		scale_fill_manual(values=rev(as.character(col_types[unique(rec_ttt$Type)])))+
		scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
		theme(legend.position="none", panel.background=element_blank(), 
		axis.title=element_blank(), axis.text=element_blank(), 
		axis.ticks=element_blank())+tag_thm
	pccm <- wrap_plots(A=plsc[[1]], B=plsc[[2]], C=plsc[[3]], D=plsc[[4]], E=plsc[[5]], F=pccc, G=pccd, 
		design="GA\nGB\nGC\nGD\nGE\n#F", heights=c(1,1,1,1,1,2), widths=c(1,30))
	return(wrap_elements(wrap_plots(A=pcca, B=pccb, C=pccm, design=c("AC\nBC"), 
		heights=c(40, 1))+plot_annotation(NULL, theme=theme(plot.margin=margin(0,0,-32,10), panel.spacing=unit(0, "pt"))))+tag_thm)
})

ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(paa,pab,pac,pad), nrow=1, widths=c(0.8,0.6,0.6,1.5))+
	plot_annotation(tag_levels=list(c("A", "B", "C", "D")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pba,pbb,pbc,pbd), nrow=1, widths=c(1.5,0.8,0.8,1))+
	plot_annotation(tag_levels=list(c("E", "", "F", "G", "H")), theme=tag_thm))+tag_thm, 
	wrap_elements(wrap_plots(list(pca, pcb, plds[[1]]), nrow=1, widths=c(1, 1, 2))+
	plot_annotation(tag_levels=list(c("I", "J", "", "M")), theme=tag_thm2))+tag_thm, 
	wrap_elements(wrap_plots(list(pda,pdb, plds[[2]]), nrow=1, widths=c(1, 1, 2))+
	plot_annotation(tag_levels=list(c("K", "L", "N")), theme=tag_thm))+tag_thm), 
	ncol=1, heights=c(1,1,1,1)), width=13, height=14, dpi=200, filename="oe_fl_fig_02.png", limitsize=F)


#cmp_raw_filter <- rec_ca[which(rec_ca$Group != "N.S."),]
#ego_total <- data.frame()
#for (type in 1:4)
#{
#	terms <- cmp_raw_filter$Symbol[which(cmp_raw_filter$Type == type)]
#	eid <- mapIds(org.Mm.eg.db, keys=terms, column="ENTREZID", keytype="SYMBOL", multiVals="first")
#	ego <- enrichGO(gene=eid, keyType="ENTREZID", OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH", readable=T)
#	ego <- ego@result[which((ego@result$pvalue < 0.05 | ego@result$qvalue < 0.05) & ego@result$Count > 1),]
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
#	ego <- ego[order(ego$p.adjust),]
#	ego <- ego[order(ego$Count, decreasing=T),]
#	#if (nrow(ego) > 5) ego <- ego[1:5,]
#	if (nrow(ego_total) == 0) ego_total <- ego else ego_total <- rbind(ego_total, ego)
#}
#write.csv(ego_total, "cmp_raw_go2.csv")

