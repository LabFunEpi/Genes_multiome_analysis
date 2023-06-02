library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
keepBSgenomeSequences <- function(genome, seqnames)
{
    stopifnot(all(seqnames %in% seqnames(genome)))
    genome@user_seqnames <- setNames(seqnames, seqnames)
    genome@seqinfo <- genome@seqinfo[seqnames]
    genome
}
sequences_to_keep <- paste0("chr", c(1:22, "X", "Y"))
genome <- keepBSgenomeSequences(genome, sequences_to_keep)

setwd("/COMPARISON")
colors.mn <- c(`CD4 T` = "#99CCCC", `Treg` = "#E794EA", `CD8 T` = "#A9E0AB", `gdT` = "#CE5B5B", `MAIT` = "#EDA6A6", `NK` = "#E2E8A6", `B` = "#D68645", `Plasmablast` = "#48B758", `CD14 Mono` = "#B4AEE5", `CD16 Mono` = "#F9BD95", `Int Mono` = "#B81254", `cDC` = "#67A5E5", `pDC` = "#D8E1A7", `HSPC` = "#C683ED", `Neutrophils` = "#3B69B7", `Basophils` = "black")

GEX <- readRDS("/FINAL/scRNA/filtered.rds")
pred.GEX.monaco <- readRDS("/FINAL/scRNA/SingleR.monaco.fine.rds")
GEX$monaco <- pred.GEX.monaco$new.fine

MULT <- readRDS("/FINAL/scMULT/filtered.rds")

sub_samples <- GEX$sample %>% table() %>% data.frame() %>% set_colnames(c("sample", "GEX")) %>% full_join(MULT$sample %>% table() %>% data.frame() %>% set_colnames(c("sample", "MULT"))) %>% drop_na() %>% select(sample) %>% unlist() %>% as.character()


###################################################################################

print(paste0(Sys.time(), ":: Cell-type identification using SingleR for MULT ... "))

library(celldex)
monaco <- celldex::MonacoImmuneData()
colors.monaco <- c(`CD4 T` = "#81a66e", `Treg` = "#32584e", `CD8 T` = "#f5df87", `gdT` = "#b4c29f", `MAIT` = "#fdbd92", `NK` = "#546046", `B` = "#7c6289", `Plasmablast` = "#5991ba", `CD14 Mono` = "#d18c5a", `CD16 Mono` = "#f18e8c", `Int Mono` = "#8e402b", `cDC` = "#65a84c", `pDC` = "#b5d75b", `HSPC` = "#ffcbe7", `Neutrophils` = "#a2bcc2", `Basophils` = "#66a5a4")
fine_map <- data.frame(colData(monaco)) %>% select(-label.ont) %>% distinct() %>% 
    mutate(new.fine = case_when(
        label.fine == "MAIT cells" ~ "MAIT", 
        (label.fine == "Vd2 gd T cells" | label.fine == "Non-Vd2 gd T cells") ~ "gdT",
        label.fine == "T regulatory cells" ~ "Treg", 
        label.fine == "Classical monocytes" ~ "CD14 Mono",
        label.fine == "Non classical monocytes" ~ "CD16 Mono",
        label.fine == "Intermediate monocytes" ~ "Int Mono",
        label.fine == "Plasmacytoid dendritic cells" ~ "pDC",
        label.fine == "Myeloid dendritic cells" ~ "cDC",
        label.fine == "Plasmablasts" ~ "Plasmablast",
        label.fine == "Progenitor cells" ~ "HSPC",
        label.fine == "Natural killer cells" ~ "NK",
        TRUE ~ label.main)) %>%
    mutate(new.fine = factor(case_when(
        new.fine == "CD8+ T cells" ~ "CD8 T",
        new.fine == "CD4+ T cells" ~ "CD4 T",
        new.fine == "B cells" ~ "B",
        TRUE ~ new.fine), levels = names(colors.monaco))) %>%
    rename(labels = label.fine)

DefaultAssay(MULT) <- "RNA"
diet_MULT <- DietSeurat(MULT, assays = "RNA")
MULT.se <- as.SingleCellExperiment(diet_MULT)
library(SingleR)
pred.MULT.monaco <- SingleR(test = MULT.se, ref = monaco, labels = monaco$label.fine, BPPARAM=MulticoreParam(8))
pred.MULT.monaco <- data.frame(pred.MULT.monaco) %>% left_join(fine_map)
saveRDS(pred.MULT.monaco, "SingleR.monaco.fine_MULT.rds")

pred.MULT.monaco <- readRDS("SingleR.monaco.fine_MULT.rds")
MULT$monaco <- pred.MULT.monaco$new.fine

############################# Markers ##############################################

GEX_sub <- subset(GEX, subset = sample %in% sub_samples)
MULT_sub <- subset(MULT, subset = sample %in% sub_samples)

DefaultAssay(GEX_sub) <- "integrated"
GEX_sub %<>% ScaleData()
GEX_sub %<>% RunPCA(npcs = 50)
GEX_sub %<>% RunUMAP(reduction = "pca", dims = 1:50)
GEX_sub %<>% FindNeighbors(reduction = "pca", dims = 1:50)
GEX_sub %<>% FindClusters(resolution = 0.2)

DefaultAssay(MULT_sub) <- "integrated_RNA"
MULT_sub %<>% ScaleData()
MULT_sub %<>% RunPCA(npcs = 50)
MULT_sub %<>% RunUMAP(reduction = "pca", dims = 1:50, reduction.name = "umap.irna", reduction.key = "irnaUMAP_")
MULT_sub %<>% FindNeighbors(reduction = "pca", dims = 1:50)
MULT_sub %<>% FindClusters(resolution = 0.4)

saveRDS(GEX_sub, "GEX_sub.rds")
saveRDS(MULT_sub, "MULT_sub.rds")

GEX_markers <- FindAllMarkers(GEX_sub, only.pos = TRUE)
MULT_markers <- FindAllMarkers(MULT_sub, only.pos = TRUE)

saveRDS(GEX_markers, "GEX_markers_uns.rds")
saveRDS(MULT_markers, "MULT_markers_uns.rds")

Idents(GEX_sub) <- "monaco"
GEX_markers <- FindAllMarkers(GEX_sub, only.pos = TRUE)
saveRDS(GEX_markers, "GEX_markers_monaco.rds")

Idents(MULT_sub) <- "monaco"
MULT_markers <- FindAllMarkers(MULT_sub, only.pos = TRUE)
saveRDS(MULT_markers, "MULT_markers_monaco.rds")

Idents(MULT_sub) <- "predicted.scRNA.monaco"
MULT_markers <- FindAllMarkers(MULT_sub, only.pos = TRUE)
saveRDS(MULT_markers, "MULT_markers_scRNA.monaco.rds")

####################################################################################

tab <- table(monaco_ref = factor(MULT$monaco, levels = names(colors.mn)), scRNA_ref = factor(MULT$predicted.scRNA.monaco, levels = names(colors.mn))) 

# temp <- data.frame(monaco_ref = factor(MULT$monaco, levels = names(colors.mn)), scRNA_ref = factor(MULT$predicted.scRNA.monaco, levels = names(colors.mn))) 

pdf(file='CellTypeIdent_1.pdf', width=8, height=8)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap::pheatmap(log10(tab+1), cluster_rows = FALSE, cluster_cols = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Using Monaco reference", y=-0.07, gp=gpar(fontsize=16))
grid.text("Using scRNA reference", x=-0.07, rot=90, gp=gpar(fontsize=16))
dev.off()

pdf(file='CellTypeIdent_2.pdf', width=20, height=8)
plotdata <- data.frame(monaco_ref = factor(MULT$monaco, levels = names(colors.mn)), scRNA_ref = factor(MULT$predicted.scRNA.monaco, levels = names(colors.mn))) %>% table() %>% data.frame() %>%
    mutate(monaco_ref = factor(monaco_ref, levels = rev(names(colors.mn))))
p1 <- ggplot(plotdata, aes(fill=scRNA_ref, y=monaco_ref, x=Freq) ) +
    geom_bar(position="fill", stat="identity") + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = colors.mn) +
    theme_cowplot()
p1
dev.off()

####################################################################################

celltype.table1 <- data.frame(celltype = GEX$monaco, sample = GEX$sample) %>% table() %>% data.frame() %>% mutate(modality = "scRNA") %>%
    bind_rows(data.frame(celltype = MULT$monaco, sample = MULT$sample) %>% table() %>% data.frame() %>% mutate(modality = "scMULT1")) %>%
    bind_rows(data.frame(celltype = MULT$predicted.scRNA.monaco, sample = MULT$sample) %>% table() %>% data.frame() %>% mutate(modality = "scMULT2")) %>%
    rename(count = Freq) %>%
    mutate(celltype = factor(celltype, levels=names(colors.mn))) %>%
    filter(sample %in% sub_samples) %>%
    mutate(sample = factor(sample, levels = c("B23", "B30", "B31", "B34", "B36", "B39", "A5", "A6", "A12", "A16", "A25", "A40", "PM001", "C2", "C3", "C4", "M4666", "C8")), modality = factor(modality, levels = rev(c("scRNA", "scMULT1", "scMULT2"))))
p1 <- ggplot(celltype.table1) + 
    geom_bar(aes(fill=celltype, y=count, x=modality), position=position_fill(reverse = TRUE), stat="identity", color = "black") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = colors.mn) +
    facet_grid(rows = vars(sample)) +
    coord_flip() + 
    theme_cowplot() + 
    theme(legend.position = "right", strip.text = element_blank(), strip.background = element_blank(), axis.title = element_blank(), axis.text.y = element_blank()) + theme(panel.spacing = unit(1, "lines"))
pdf(file='proportions.pdf', width=15, height=12)
p1
dev.off()

celltype.table3 <- data.frame(celltype = GEX$monaco, sample = GEX$sample) %>% filter(sample %in% sub_samples) %>% select(celltype) %>% table() %>% data.frame() %>% rename(scRNA = Freq) %>%
    left_join(data.frame(celltype = MULT$monaco, sample = MULT$sample) %>% filter(sample %in% sub_samples) %>% select(celltype) %>% table() %>% data.frame() %>% rename(scMULT1 = Freq)) %>%
    left_join(data.frame(celltype = MULT$predicted.scRNA.monaco, sample = MULT$sample) %>% filter(sample %in% sub_samples) %>% select(celltype) %>% table() %>% data.frame() %>% rename(scMULT2 = Freq)) %>%
    mutate(celltype = factor(celltype, levels=rev(names(colors.mn)))) %>%
    replace_na(list(scMULT2 = 0)) %>% 
    mutate(scRNA = scRNA / sum(scRNA), scMULT1 = scMULT1 / sum(scMULT1), scMULT2 = scMULT2 / sum(scMULT2)) %>%
    pivot_longer(!c(celltype), names_to = "group", values_to = "freq") %>% 
    mutate(group = factor(group, levels = c("scRNA", "scMULT1", "scMULT2"))) %>%
    mutate(freq = rescale(freq, to = c(0, 1), from = c(0, 0.25))) %>%
    pivot_wider(names_from = "celltype", values_from = "freq")

pdf(file='proportions_spider.pdf', width=12, height=12)
ggradar(celltype.table3, values.radar = c(0, 0.125, 0.25))
dev.off()

# GEX_sub <- readRDS("GEX_sub.rds")
# MULT_sub <- readRDS("MULT_sub.rds")

rnaobj <- GEX_sub
multobj <- MULT_sub
p1 <- ggplot(data.frame(rnaobj[["umap"]][[]], label = factor(rnaobj$monaco)), aes(x=UMAP_1, y=UMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.3, stroke=0.1, shape=16), dpi = 400) +
    scale_color_manual(values = colors.mn) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + ggtitle("scRNA using Monaco ref")
p2 <- table(rnaobj$monaco) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn)))) %>%
    ggplot(aes(fill=Var1, y=Freq, x="Blah")) +
    geom_bar(position="fill", stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) +
    coord_flip() + theme_cowplot() + 
    theme(aspect.ratio=0.1, legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p3 <- ggplot(data.frame(multobj[["umap.irna"]][[]], label = factor(multobj$predicted.scRNA.monaco)), aes(x=-irnaUMAP_1, y=irnaUMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16), dpi = 400) +
    scale_color_manual(values = colors.mn) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + ggtitle("Multiome (RNA) using Monaco ref")
p4 <- table(multobj$predicted.scRNA.monaco) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn)))) %>%
    ggplot(aes(fill=Var1, y=Freq, x="Blah")) +
    geom_bar(position="fill", stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) +
    coord_flip() + theme_cowplot() + 
    theme(aspect.ratio=0.1, legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
p5 <- ggplot(data.frame(multobj[["umap.irna"]][[]], label = factor(multobj$monaco)), aes(x=-irnaUMAP_1, y=irnaUMAP_2, color=label)) + 
    ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16), dpi = 400) +
    scale_color_manual(values = colors.mn) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_cowplot() + theme(legend.position = "right", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()) + ggtitle("Multiome (RNA) using scRNA ref")
p6 <- table(multobj$monaco) %>% data.frame() %>% mutate(Var1 = factor(Var1, levels = rev(names(colors.mn)))) %>%
    ggplot(aes(fill=Var1, y=Freq, x="Blah")) +
    geom_bar(position="fill", stat="identity", width = 0.6) + 
    scale_fill_manual(values = colors.mn) +
    coord_flip() + theme_cowplot() + 
    theme(aspect.ratio=0.1, legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

pdf(file='umaps_sub.pdf', width=15, height=6)
(p1 / p2) | (p3 / p4) | (p5 / p6)
dev.off()

GEX_markers_monaco <- readRDS("GEX_markers_monaco.rds") %>% filter(p_val_adj < 0.05)
MULT_markers_monaco <- readRDS("MULT_markers_monaco.rds") %>% filter(p_val_adj < 0.05)
MULT_markers_scRNA_monaco <- readRDS("MULT_markers_scRNA.monaco.rds") %>% filter(p_val_adj < 0.05)

GEX_markers_monaco <- sapply(names(colors.mn), function(x){GEX_markers_monaco %>% filter(cluster == x) %>% select(gene) %>% unlist() %>% unname()})
MULT_markers_monaco <- sapply(names(colors.mn), function(x){MULT_markers_monaco %>% filter(cluster == x) %>% select(gene) %>% unlist() %>% unname()})
MULT_markers_scRNA_monaco <- sapply(names(colors.mn), function(x){MULT_markers_scRNA_monaco %>% filter(cluster == x) %>% select(gene) %>% unlist() %>% unname()})

temp1 <- sapply(GEX_markers_monaco, function(x){
    sapply(MULT_markers_monaco, function(y){
        length(intersect(x, y))
    })
})
temp2 <- sapply(GEX_markers_monaco, function(x){
    sapply(MULT_markers_scRNA_monaco, function(y){
        length(intersect(x, y))
    })
})
pdf(file='marker_analysis1.pdf', width=8, height=8)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap::pheatmap(log10(temp1+1), cluster_rows = FALSE, cluster_cols = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Using Monaco reference", y=-0.07, gp=gpar(fontsize=16))
grid.text("Using scRNA reference", x=-0.07, rot=90, gp=gpar(fontsize=16))
dev.off()
pdf(file='marker_analysis2.pdf', width=8, height=8)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap::pheatmap(log10(temp2+1), cluster_rows = FALSE, cluster_cols = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Using Monaco reference", y=-0.07, gp=gpar(fontsize=16))
grid.text("Using scRNA reference", x=-0.07, rot=90, gp=gpar(fontsize=16))
dev.off()

##############################################################################

DefaultAssay(GEX_sub) <- "RNA"
DefaultAssay(MULT_sub) <- "RNA"
GEX_sub <- DietSeurat(GEX_sub, assays = "RNA")
MULT_sub <- DietSeurat(MULT_sub, assays = "RNA")

GEX_sub <- subset(GEX_sub, features = intersect(rownames(GEX_sub@assays$RNA@counts), rownames(MULT_sub@assays$RNA@counts)))
MULT_sub <- subset(MULT_sub, features = intersect(rownames(GEX_sub@assays$RNA@counts), rownames(MULT_sub@assays$RNA@counts)))

library(muscat)
sce1 <- SummarizedExperiment(assays=list(counts=GEX_sub@assays$RNA@counts), colData=GEX_sub@meta.data)
sce1 <- as(sce1, "SingleCellExperiment")
(sce1 <- prepSCE(sce1, 
                kid = "sample", # subpopulation assignments
                gid = "status",  # group IDs
                sid = "monaco",   # sample IDs 
                drop = TRUE))  # drop all other colData columns

sce2 <- SummarizedExperiment(assays=list(counts=MULT_sub@assays$RNA@counts), colData=MULT_sub@meta.data)
sce2 <- as(sce2, "SingleCellExperiment")
(sce2 <- prepSCE(sce2, 
                kid = "sample", # subpopulation assignments
                gid = "status",  # group IDs
                sid = "monaco",   # sample IDs 
                drop = TRUE))  # drop all other colData columns
sce3 <- SummarizedExperiment(assays=list(counts=MULT_sub@assays$RNA@counts), colData=MULT_sub@meta.data)
sce3 <- as(sce3, "SingleCellExperiment")
(sce3 <- prepSCE(sce3, 
                kid = "sample", # subpopulation assignments
                gid = "status",  # group IDs
                sid = "predicted.scRNA.monaco",   # sample IDs 
                drop = TRUE))  # drop all other colData columns

t(table(sce1$sample_id))
t(table(sce2$sample_id))
t(table(sce3$sample_id))

pb1 <- aggregateData(sce1,
    assay = "counts", fun = "sum",
    by = "sample_id")
pb2 <- aggregateData(sce2,
    assay = "counts", fun = "sum",
    by = "sample_id")
pb3 <- aggregateData(sce3,
    assay = "counts", fun = "sum",
    by = "sample_id")

colnames(pb1) <- paste0("scRNA_", colnames(pb1))
colnames(pb2) <- paste0("snRNA_", colnames(pb2))
colnames(pb3) <- paste0("snRNA_", colnames(pb3))

pb12 <- cbind(pb1, pb2, deparse.level=1)
s1 <- CreateSeuratObject(counts = assay(pb12))
s1 <- NormalizeData(s1) %>% ScaleData()
s1 <- RunPCA(s1, npcs = 15, features = rownames(s1))
forplot1 <- s1[["pca"]][[]] %>% as.data.frame() %>% rownames_to_column() %>% separate(col = "rowname", sep = "_", into = c("platform", "celltype"), remove = FALSE) %>% mutate(celltype = factor(celltype, levels = names(colors.mn)))

pb13 <- cbind(pb1, pb3, deparse.level=1)
s1 <- CreateSeuratObject(counts = assay(pb13))
s1 <- NormalizeData(s1) %>% ScaleData()
s1 <- RunPCA(s1, npcs = 15, features = rownames(s1))
forplot2 <- s1[["pca"]][[]] %>% as.data.frame() %>% rownames_to_column() %>% separate(col = "rowname", sep = "_", into = c("platform", "celltype"), remove = FALSE) %>% mutate(celltype = factor(celltype, levels = names(colors.mn)))

pdf(file='pca.pdf', width=12, height=10)
p1 <- ggplot(forplot1, aes(x = PC_1, y = PC_2)) +
    geom_point(aes(color = celltype, shape = platform), size = 3) +
    scale_color_manual(values = colors.mn) +
    theme_cowplot() + theme(legend.position = "none")
p2 <- ggplot(forplot1, aes(x = PC_3, y = PC_4)) +
    geom_point(aes(color = celltype, shape = platform), size = 3) +
    scale_color_manual(values = colors.mn) +
    theme_cowplot()
p3 <- ggplot(forplot2, aes(x = PC_1, y = PC_2)) +
    geom_point(aes(color = celltype, shape = platform), size = 3) +
    scale_color_manual(values = colors.mn) +
    theme_cowplot() + theme(legend.position = "none")
p4 <- ggplot(forplot2, aes(x = PC_3, y = PC_4)) +
    geom_point(aes(color = celltype, shape = platform), size = 3) +
    scale_color_manual(values = colors.mn) +
    theme_cowplot()
(p1 | p2) / (p3 | p4)
dev.off()

####################################################################################

# sub_samples <- GEX$sample %>% table() %>% data.frame() %>% set_colnames(c("sample", "GEX")) %>% full_join(MULT$sample %>% table() %>% data.frame() %>% set_colnames(c("sample", "MULT"))) %>% drop_na() %>% select(sample) %>% unlist() %>% as.character()

# GEX_sub <- subset(GEX, subset = sample %in% sub_samples)
# MULT_sub <- subset(MULT, subset = sample %in% sub_samples)

GEX_sub$sample.celltype <- paste(GEX_sub$sample, GEX_sub$monaco, sep = ".")
MULT_sub$sample.celltype <- paste(MULT_sub$sample, MULT_sub$predicted.scRNA.monaco, sep = ".")

GEX_avg <- AverageExpression(GEX_sub, assays = "RNA", group.by = "sample.celltype", return.seurat = FALSE)$RNA %>% 
    data.frame(check.names = FALSE) %>% rownames_to_column("gene") %>%
    pivot_longer(!gene, names_to = "sample.celltype", values_to = "GEX")

MULT_avg <- AverageExpression(MULT_sub, assays = "RNA", group.by = "sample.celltype", return.seurat = FALSE)$RNA %>% 
    data.frame(check.names = FALSE) %>% rownames_to_column("gene") %>%
    pivot_longer(!gene, names_to = "sample.celltype", values_to = "MULT")

avg_expr <- GEX_avg %>% inner_join(MULT_avg) %>% filter(!(GEX == 0 & MULT == 0))
gene_list <- avg_expr %>% group_by(gene) %>% summarize(n = n()) %>% filter(n > 50) %>% select(gene) %>% unlist() %>% unname()

corr_p <- function(x, y){
    cor.test(x, y, method = "spearman")[["p.value"]]
}

corrs <- avg_expr %>% filter(gene %in% gene_list) %>% group_by(gene) %>% summarize(corr = cor(GEX, MULT, method = "spearman"), corrp = corr_p(GEX, MULT)) %>% drop_na()

pdf(file='corr_histogram_sig.pdf', width=4, height=4)
ggplot(corrs %>% filter(corrp < 0.05), aes(x = corr)) +
    geom_histogram(bins = 100, fill = "white", color = "black") +
    geom_vline(xintercept = 0) +
    theme_cowplot()
dev.off()

write.table(corrs, file = "corrs.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

all.markers.seurat <- readRDS("/mforge/research/labs/microbiome/Beyond_DNA/shared/analysis/COMPARISON/GEX_markers_monaco.rds") %>%
    filter(p_val_adj < 0.05)
avg_expr <- GEX_avg %>% inner_join(MULT_avg) %>% filter(gene %in% unique(all.markers.seurat$gene))
corrs <- avg_expr %>% group_by(gene) %>% summarize(corr = cor(GEX, MULT, method = "spearman"), corrp = corr_p(GEX, MULT)) %>% drop_na()
pdf(file='corr_histogram_sig_1.pdf', width=4, height=4)
ggplot(corrs %>% filter(corrp < 0.05), aes(x = corr)) +
    geom_histogram(bins = 100, fill = "white", color = "black") +
    geom_vline(xintercept = 0) +
    theme_cowplot()
dev.off()

write.table(corrs %>% arrange(-corr), file = "corrs_1.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


################################################

plotdata <- avg_expr %>% separate(col=sample.celltype, sep="\\.", into = c("sample", "celltype"))
expr_scatter <- function(gene_in){
    p1 <- ggplot(plotdata %>% filter(gene == gene_in), aes(x = GEX, y = MULT, color = celltype)) +
        geom_point(size=1.5, stroke=0.1, shape=16) +
        scale_color_manual(values = colors.mn) +
        guides(color = guide_legend(override.aes = list(size=5))) +
        theme_cowplot() + ggtitle(gene_in) + theme(legend.position = "none")
    return(p1)
}

genes <- (corrs %>% filter(corr > 0) %>% arrange(corrp) %>% head(n = 3))$gene
pdf(file='expr_scatter.pdf', width=12, height=4)
wrap_plots(lapply(genes, expr_scatter), ncol = 3)
dev.off()

posgenes <- corrs %>% filter(corrp < 0.05 & corr > 0) %>% select(gene) %>% unlist() %>% unname()
neggenes <- corrs %>% filter(corrp < 0.05 & corr < 0) %>% select(gene) %>% unlist() %>% unname()

Idents(GEX_sub) <- "All"
expr <- AverageExpression(GEX_sub, assays = "RNA", return.seurat = FALSE)$RNA %>% 
    data.frame() %>% rownames_to_column("gene") %>% 
    set_colnames(c("gene", "GEX"))
Idents(MULT_sub) <- "All"
expr_MULT <- AverageExpression(MULT_sub, assays = "RNA", return.seurat = FALSE)$RNA %>% 
    data.frame() %>% rownames_to_column("gene") %>% 
    set_colnames(c("gene", "MULT"))
expr <- expr %>% left_join(expr_MULT)

expr_1 <- expr %>%
    pivot_longer(!gene, names_to = "modality", values_to = "expr")

pdf(file='avg_expr_distribution.pdf', width=6, height=6)
ggplot(expr_1, aes(x = expr, fill = modality)) +
    geom_histogram(binwidth = 0.01, color = "black") +
    coord_cartesian(xlim = c(0, 1)) +
    scale_color_manual(values = c("red", "blue")) +
    theme_cowplot()
dev.off()

expr_2 <- expr_1 %>% filter(gene %in% neggenes) %>% mutate(correl = "Negative") %>%
    bind_rows(expr_1 %>% filter(gene %in% posgenes) %>% mutate(correl = "Positive"))

pdf(file='corr_genes_avg_expr.pdf', width=3, height=4)
ggplot(expr_2, aes(x = correl, y = expr)) +
    ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
    geom_boxplot(outlier.shape = NA, color = "blue") +
    stat_compare_means(label = "p.format") +
    facet_grid(cols = vars(modality)) +
    scale_y_continuous(breaks = c(0, 10, 30, 50, 100, 150, 200)) +
    coord_trans(y = scales::log1p_trans()) +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()










