setwd("/MultiomePaper/")

tab <- read.table("scMULT_targeting_efficiency.txt", skip = 1, col.names = c("sample", "tissue", "cellcount", "viability", "concentration", "targeted", "captured"))

tab1 <- tab %>% mutate(cap_eff = (captured*100)/targeted) %>%
    # filter(!(sample %in% c("B34", "B39", "A40", "C5"))) %>% 
    mutate(targeted = factor(targeted)) %>%
    mutate(viability_range = cut(viability, breaks=c(60, 70, 80, 90, 100), labels=c("60-70","70-80","80-90","90-100")))

cols <- OkabeIto40[-4][1:6]
names(cols) <- levels(tab1$targeted)
pdf(file='capture_efficiency_paper.pdf', width=10, height=4)
p1 <- ggplot(tab1 %>% filter(tissue == "PBMC"), aes(x = cap_eff, y = concentration, color = targeted)) +
    geom_point(size = 3) +
    # geom_text_repel(aes(label = sample), size = 2) +
    lims(x = c(0, 350), y = c(0, 4000)) +
    scale_color_manual(values=cols) +
    theme_cowplot()
p2 <- ggplot(tab1 %>% filter(tissue == "OVCAR"), aes(x = cap_eff, y = concentration, color = targeted)) +
    geom_point(size = 3) +
    # geom_text_repel(aes(label = sample), size = 2) +
    lims(x = c(0, 400), y = c(0, 8000)) +
    scale_color_manual(values=cols) +
    theme_cowplot()
p1 | p2
dev.off()

tab1 <- read.table("MultiomeCellRangerSummary.txt", header = TRUE) %>% 
    pivot_longer(!Statistic, names_to = "Sample", values_to = "value") %>%
    pivot_wider(names_from = "Statistic", values_from = "value")
tab2 <- read.table("MultiomeMeta.txt", header = TRUE, na.strings = "NA")
tab1 %<>% left_join(tab2) %>% mutate(across(Estimated.number.of.cells:GEX.Valid.barcodes, as.numeric)) %>%
    mutate(Biological = factor(Biological, levels = c("PBMC", "OC_TISSUE"))) %>%
    mutate(Trials = factor(Trials, c("Collagenase", "ComplexTissue", "IncuTime_5", "IncuTime_10", "IncuTime_30", "IncuTime_10_FACS", "IncuTime_30_NOFACS")))

plotdata <- tab1 %>% mutate(Estimated.number.of.cells = Estimated.number.of.cells / 1000, ATAC.Mean.raw.read.pairs.per.cell = ATAC.Mean.raw.read.pairs.per.cell / 1000, ATAC.Number.of.peaks = ATAC.Number.of.peaks / 1000, GEX.Mean.raw.reads.per.cell = GEX.Mean.raw.reads.per.cell / 1000, GEX.Median.UMI.counts.per.cell = GEX.Median.UMI.counts.per.cell / 1000, GEX.Median.genes.per.cell = GEX.Median.genes.per.cell / 1000, GEX.Total.genes.detected = GEX.Total.genes.detected / 1000)

pdf(file='qc_by_sample.pdf', width=14, height=8)
variables <- colnames(plotdata)[c(4, 10, 8, 24, 17, 27, 29, 30, 31, 46)]
thresholds <- c(0, 0.4, 0, 5, 0, 0.6, 0.1, 0, 0, 0)
ax_labels <- c("Estimated number of cells (x 10^3)", "ATAC Fraction of high quality fragments in cells", "ATAC Confidently mapped read pairs", "ATAC TSS enrichment score", "ATAC Number of peaks (x 10^3)", "GEX Fraction of transcriptomic reads in cells", "GEX Median UMI counts per cell (x 10^3)", "GEX Median genes per cell (x 10^3)", "GEX Percent duplicates", "GEX Total genes detected (x 10^3)")
wrap_plots(mapply(function(variable, threshold, axlabel, setbreaks){
    variable <- as.name(variable)
    p <- ggplot(plotdata, aes(x = Group, y = {{variable}})) +
        geom_boxplot(outlier.shape = NA) +
        geom_beeswarm(size=2, stroke = 0, shape = 16, color = "black") +
        facet_grid(cols = vars(Biological), scales = "free", space = "free") +
        ylab(str_wrap(axlabel, width = 30)) +
        scale_y_continuous(limits = c(0, NA)) +
        theme_cowplot() + theme(axis.title.x = element_blank(), strip.background = element_blank(), strip.text = element_blank())
    if (threshold != 0){p = p + geom_hline(yintercept = threshold, color = "red", linetype = "dashed")}
    return(p)
}, variables, thresholds, ax_labels, SIMPLIFY = FALSE), ncol = 5, guides = "collect")
dev.off()

plotdata <- tab1 %>% mutate(ATAC.Sequenced.read.pairs = ATAC.Sequenced.read.pairs / 1000000, GEX.Sequenced.read.pairs = GEX.Sequenced.read.pairs / 1000000, ATAC.Mean.raw.read.pairs.per.cell = ATAC.Mean.raw.read.pairs.per.cell / 1000, GEX.Mean.raw.reads.per.cell = GEX.Mean.raw.reads.per.cell / 1000)
pdf(file='sequencing_depth.pdf', width=6, height=8)
p1 <- ggplot(plotdata, aes(x = Group, y = ATAC.Sequenced.read.pairs)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(size=2, stroke = 0, shape = 16, color = "black") +
    facet_grid(cols = vars(Biological), scales = "free", space = "free") +
    ylab(str_wrap("ATAC Sequenced read pairs (x 10^6)", width = 30)) +
    scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 50)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), strip.background = element_blank(), strip.text = element_blank())
p2 <- ggplot(plotdata, aes(x = Group, y = GEX.Sequenced.read.pairs)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(size=2, stroke = 0, shape = 16, color = "black") +
    facet_grid(cols = vars(Biological), scales = "free", space = "free") +
    ylab(str_wrap("GEX Sequenced read pairs (x 10^6)", width = 30)) +
    scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 50)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), strip.background = element_blank(), strip.text = element_blank())
p3 <- ggplot(plotdata, aes(x = Group, y = ATAC.Mean.raw.read.pairs.per.cell)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(size=2, stroke = 0, shape = 16, color = "black") +
    facet_grid(cols = vars(Biological), scales = "free", space = "free") +
    ylab(str_wrap("ATAC Mean raw read pairs per cell (x 10^3)", width = 30)) +
    scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 50)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), strip.background = element_blank(), strip.text = element_blank())
p4 <- ggplot(plotdata, aes(x = Group, y = GEX.Mean.raw.reads.per.cell)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(size=2, stroke = 0, shape = 16, color = "black") +
    facet_grid(cols = vars(Biological), scales = "free", space = "free") +
    ylab(str_wrap("GEX Mean raw reads per cell (x 10^3)", width = 30)) +
    scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 50)) +
    theme_cowplot() + theme(axis.title.x = element_blank(), strip.background = element_blank(), strip.text = element_blank())
wrap_plots(list(p1, p2, p3, p4), ncol = 2)
dev.off()

plotdata <- tab1 %>% filter(Trials != "Other") %>% select(TrialSample, Trials, Group, ATAC.Fraction.of.high.quality.fragments.in.cells, GEX.Fraction.of.transcriptomic.reads.in.cells) %>% pivot_longer(!c(TrialSample, Trials, Group), names_to = "statistic", values_to = "value")
pdf(file='trial_samples.pdf', width=6, height=5)
p1 <- ggplot(plotdata %>% filter(TrialSample == "PH1199"), aes(x = Trials, y = value, color = statistic)) +
    geom_point(size=2, stroke = 0, shape = 16) +
    ylab("Fraction") +
    coord_cartesian(ylim = c(0, 1)) +
    ggtitle("PH1199") +
    theme_cowplot() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(plotdata %>% filter(TrialSample == "PH1214"), aes(x = Trials, y = value, color = statistic)) +
    geom_point(size=2, stroke = 0, shape = 16) +
    ylab("Fraction") +
    coord_cartesian(ylim = c(0, 1)) +
    ggtitle("PH1214") +
    theme_cowplot() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 <- ggplot(plotdata %>% filter(TrialSample == "PH1217"), aes(x = Trials, y = value, color = statistic)) +
    geom_point(size=2, stroke = 0, shape = 16) +
    ylab("Fraction") +
    facet_grid(cols = vars(Group), scales = "free", space = "free") +
    coord_cartesian(ylim = c(0, 1)) +
    ggtitle("PH1217") +
    theme_cowplot() + theme(legend.position = "none", axis.title.x = element_blank(), strip.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1 | p2 | p3
dev.off()

################################################################################

annotations <- readRDS("/mforge/research/labs/microbiome/Beyond_DNA/shared/analysis/FINAL/scMULT/annotations.rds")
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
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
library(BiocParallel)
register(SerialParam())
library(celldex)
library(SingleR)

source("SignacFunctions_mod.R")

############ Merge ##########

datadir1 <- "/COVIDCHIP/"
samplenames1 <- c("B23", "B30", "B31", "B34", "B36", "B39", "A5", "A6", "A12", "A16", "A25", "A40", "C1", "C2", "C3", "C4", "C5", "C8")
datadir2 <- "/OvCancerMULT/"
samplenames2 <- c("PH1177", "PH1178", "PH1183", "PH1184", "PH1194", "PH1203", "PH1199R", "PH1199N", "PH1214-5", "PH1214-30", "PH1217-10", "PH1217-30", "PH1217-S10", "PH1217-F30", "PH1230", "PH1238", "PH1239", "PH1243")

paths <- c(paste0(datadir1, samplenames1, "/outs/seurat_object.rds"), paste0(datadir2, samplenames2, "/outs/seurat_object.rds"))
objs <- lapply(paths, readRDS)

names(objs) <- c(samplenames1, samplenames2)
# Add sample names
add.samplenames <- function(obj, samplename){
    obj$sample <- samplename
    return(obj)
}
objs <- mapply(add.samplenames, objs, names(objs))
objs <- lapply(objs, function(x){RenameCells(x, add.cell.id = x$sample[[1]])})

objs_merged <- merge(objs[[1]], y = objs[2:length(objs)], project = "MULT")

saveRDS(objs_merged, "objs_merged.rds")

objs_merged <- NucleosomeSignal(object = objs_merged, assay = "ATAC")
objs_merged <- TSSEnrichment(object = objs_merged, fast = FALSE, assay = "ATAC")
objs_merged$project <- "MULT"

saveRDS(objs_merged, "objs_merged_1.rds")

# objs_merged <- readRDS("objs_merged_1.rds")

meta <- read.table("MultiomeMeta.txt", header = TRUE, na.strings = "NA") %>%
    mutate(Sample = str_replace(Sample, "\\.", "-")) %>%
    rename(sample = Sample)

temp <- objs_merged@meta.data %>% left_join(meta)
objs_merged$Group <- temp$Group

DefaultAssay(objs_merged) <- "ATAC"
Idents(objs_merged) <- "Group"

TSSPlottab <- TSSPlotdata(objs_merged)
FragmentHistogramtab <- FragmentHistogramdata(objs_merged)
TSSPlottab %<>% mutate(group = factor(group, levels = c("A", "B", "C", "D", "E")))
FragmentHistogramtab %<>% mutate(group = factor(group, levels = c("A", "B", "C", "D", "E")))

pdf(file='TSSPlot_FragmentHistogram.pdf', width=8, height=6)
p1 <- ggplot(TSSPlottab %>% filter(group == "A"), aes(x = position, y = norm.value)) +
    geom_line(stat = "identity", linewidth = 0.2) +
    xlab("Distance from TSS (bp)") +
    ylab(label = "Mean TSS enrichment score") +
    # coord_cartesian(ylim = c(0, 12.5)) +
    theme_cowplot()
p2 <- ggplot(TSSPlottab %>% filter(group != "A"), aes(x = position, y = norm.value, color = group)) +
    geom_line(stat = "identity", linewidth = 0.2) +
    xlab("Distance from TSS (bp)") +
    # coord_cartesian(ylim = c(0, 12.5)) +
    theme_cowplot() + theme(axis.title.y = element_blank())
p3 <- ggplot(data = FragmentHistogramtab %>% filter(group == "A"), aes(length)) +
    geom_freqpoly(bins = 1200) + 
    scale_y_continuous(labels = function(l) {trans = l / 1000}) +
    coord_cartesian(xlim = c(0, 800)) +
    theme_cowplot() +
    xlab("Fragment length (bp)") +
    ylab("Count (x 10^3)")
p4 <- ggplot(data = FragmentHistogramtab %>% filter(group != "A"), aes(length, color = group)) +
    geom_freqpoly(bins = 1200) + 
    scale_y_continuous(labels = function(l) {trans = l / 1000}) +
    coord_cartesian(xlim = c(0, 800)) +
    theme_cowplot() +
    xlab("Fragment length (bp)") +
    theme(axis.title.y = element_blank())
wrap_plots(p3, p4, p1, p2, ncol = 2)
dev.off()

plotdata <- data.frame(
    ident = "all", 
    nFeature_RNA = objs_merged$nFeature_RNA, 
    nCount_RNA = objs_merged$nCount_RNA,
    percent.mito = objs_merged$percent.mt,
    nFeature_ATAC = objs_merged$nFeature_ATAC, 
    nCount_ATAC = objs_merged$nCount_ATAC,
    TSS.enrichment = objs_merged$TSS.enrichment,
    nucleosome_signal = objs_merged$nucleosome_signal,
    sample = objs_merged$sample,
    group = objs_merged$Group
)

options(scipen = 999)
p1 <- ggplot(plotdata, aes(x=group, y=nFeature_RNA)) + 
        ggrastr::rasterise(geom_jitter(size = 0.05, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "blue") + 
        # scale_color_manual(values = OkabeIto40) +
        scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000)) +
        coord_trans(y = scales::log1p_trans()) +
        ylab("No. of unique genes detected (RNA)") +
        theme_cowplot()
p2 <- ggplot(plotdata, aes(x=group, y=nFeature_ATAC)) + 
        ggrastr::rasterise(geom_jitter(size = 0.05, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "blue") +
        # scale_color_manual(values = OkabeIto40) +
        scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000, 100000)) +        
        coord_trans(y = scales::log1p_trans()) +
        ylab("No. of unique peaks detected (ATAC)") +
        theme_cowplot()
p3 <- ggplot(plotdata, aes(x=group, y=percent.mito)) + 
        ggrastr::rasterise(geom_jitter(size = 0.05, stroke = 0, shape = 16, show.legend = TRUE), dpi = 400) +
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "blue") + 
        # scale_color_manual(values = OkabeIto40) +
        scale_y_continuous(breaks = seq(0, 100, 10)) +
        ylab("Percent reads mapping to mito (RNA)") +
        # guides(color = guide_legend(override.aes = list(size=3))) +
        theme_cowplot() + theme(legend.text=element_text(size=5))

pdf(file='qc_by_cell.pdf', width=12, height=4)
(p1 | p2 | p3) & theme(axis.title.x=element_blank()) 
dev.off()

plotdata1 <- plotdata %>% filter(sample %in% c("PH1199R", "PH1199N", "PH1214-5", "PH1214-30", "PH1217-10", "PH1217-30", "PH1217-S10", "PH1217-F30"))

options(scipen = 1)
temp <- plotdata1 %>% filter(sample %in% c("PH1199R", "PH1199N")) %>% mutate(sample = factor(sample, levels = c("PH1199R", "PH1199N")))
p1 <- ggplot(temp, aes(x=sample, y=nFeature_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "black") + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        # stat_compare_means(aes(label = ..p.format..)) +
        scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000)) +
        coord_trans(y = "log10", ylim = c(1, 100000)) +
        geom_hline(yintercept = 200, color = "red", linetype = "dashed") +
        ylab("No. of unique genes detected (RNA)") +
        theme_cowplot() + theme(axis.title.x=element_blank(), axis.text.x = element_blank())
p2 <- ggplot(temp, aes(x=sample, y=nFeature_ATAC)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "black") + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        # stat_compare_means(aes(label = ..p.format..)) +
        scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000)) +
        coord_trans(y = "log10", ylim = c(1, 100000)) +
        geom_hline(yintercept = 200, color = "red", linetype = "dashed") +
        ylab("No. of unique peaks detected (ATAC)") +
        theme_cowplot() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
temp <- plotdata1 %>% filter(sample %in% c("PH1214-5", "PH1214-30")) %>% mutate(sample = factor(sample, levels = c("PH1214-5", "PH1214-30")))
p3 <- ggplot(temp, aes(x=sample, y=nFeature_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "black") + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        # stat_compare_means(aes(label = ..p.format..)) +
        scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000)) +
        coord_trans(y = "log10", ylim = c(1, 100000)) +
        geom_hline(yintercept = 200, color = "red", linetype = "dashed") +
        ylab("No. of unique genes detected (RNA)") +
        theme_cowplot() + theme(axis.title=element_blank(), axis.text.x = element_blank())
p4 <- ggplot(temp, aes(x=sample, y=nFeature_ATAC)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "black") + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        # stat_compare_means(aes(label = ..p.format..)) +
        scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000)) +
        coord_trans(y = "log10", ylim = c(1, 100000)) +
        geom_hline(yintercept = 200, color = "red", linetype = "dashed") +
        ylab("No. of unique peaks detected (ATAC)") +
        theme_cowplot() + theme(axis.title=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
temp <- plotdata1 %>% filter(sample %in% c("PH1217-10", "PH1217-30", "PH1217-S10", "PH1217-F30")) %>% mutate(sample = factor(sample, levels = c("PH1217-10", "PH1217-30", "PH1217-S10", "PH1217-F30")))
p5 <- ggplot(temp, aes(x=sample, y=nFeature_RNA)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "black") + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000)) +
        coord_trans(y = "log10", ylim = c(1, 100000)) +
        geom_hline(yintercept = 200, color = "red", linetype = "dashed") +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
        ylab("No. of unique genes detected (RNA)") +
        theme_cowplot() + theme(axis.title=element_blank(), axis.text.x = element_blank(), strip.background = element_blank())
p6 <- ggplot(temp, aes(x=sample, y=nFeature_ATAC)) + 
        geom_violin(draw_quantiles = seq(0, 1, 0.5), color = "black") + 
        ggrastr::rasterise(geom_jitter(size = 0.2, stroke = 0, shape = 16, show.legend = FALSE), dpi = 400) +
        # stat_compare_means(aes(label = ..p.format..)) +
        scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000)) +
        coord_trans(y = "log10", ylim = c(1, 100000)) +
        geom_hline(yintercept = 200, color = "red", linetype = "dashed") +
        facet_grid(cols = vars(group), scales = "free", space = "free") +
        ylab("No. of unique peaks detected (ATAC)") +
        theme_cowplot() + theme(axis.title=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_blank())

pdf(file='qc_by_trial.pdf', width=10, height=10)
wrap_plots(p1, p3, p5, p2, p4, p6, ncol = 3, widths = c(1, 1, 2))
dev.off()

####################################################################

datadir1 <- "/cellranger_output/scMULT_v2/"
samplenames1 <- c("C3", "B23")
datadir2 <- "/OvCancerMULT/"
samplenames2 <- c("PH1199R", "PH1214-5", "PH1214-30", "PH1217-F30", "PH1230")

paths <- c(paste0(datadir1, samplenames1), paste0(datadir2, samplenames2))
projection <- "umap"

cellranger_umap <- function(samplepath){
    atac_projections <- read.csv(paste0(samplepath, "/outs/analysis/dimensionality_reduction/atac/", projection, "_projection.csv"))
    gex_projections <- read.csv(paste0(samplepath, "/outs/analysis/dimensionality_reduction/gex/", projection, "_projection.csv"))
    atac_clusters <- read.csv(paste0(samplepath, "/outs/analysis/clustering/atac/graphclust/clusters.csv"))
    gex_clusters <- read.csv(paste0(samplepath, "/outs/analysis/clustering/gex/graphclust/clusters.csv"))
    atac_projections <- atac_projections %>% left_join(atac_clusters) %>% mutate(Cluster = factor(Cluster))
    gex_projections <- gex_projections %>% left_join(gex_clusters) %>% mutate(Cluster = factor(Cluster))

    xaxis <- as.name(paste0(toupper(projection), ".1"))
    yaxis <- as.name(paste0(toupper(projection), ".2"))
    p1 <- ggplot(atac_projections, aes(x = {{xaxis}}, y = {{yaxis}}, color = Cluster)) +
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = OkabeIto40) +
        guides(color = guide_legend(override.aes = list(size=3))) +
        ggtitle("ATAC") +
        theme_cowplot() + theme(aspect.ratio=1)
    p2 <- ggplot(gex_projections, aes(x = {{xaxis}}, y = {{yaxis}}, color = Cluster)) +
        ggrastr::rasterise(geom_point(size=0.5, stroke=0.1, shape=16), dpi = 400) +
        scale_color_manual(values = OkabeIto40) +
        guides(color = guide_legend(override.aes = list(size=3))) +
        ggtitle("GEX") +
        theme_cowplot() + theme(aspect.ratio=1)
    return((p1 | p2) + theme(plot.background = element_rect(fill = NA, colour = 'black', size = 1)))
}

plots <- lapply(paths, cellranger_umap)

pdf(file='projections.pdf', width=20, height=20)
layout <- '
AB
C#
DE
FG
'
names(plots) <- LETTERS[1:7]
wrap_plots(plots, design = layout)
dev.off()


datadir2 <- "/mforge/research/labs/experpath/maia/m237371/OvCancerMULT/"
samplenames2 <- c("PH1199R", "PH1199N", "PH1214-5", "PH1214-30", "PH1217-10", "PH1217-30", "PH1217-S10", "PH1217-F30")

paths <- paste0(datadir2, samplenames2)
projection <- "umap"

plots <- lapply(paths, cellranger_umap)

pdf(file='projections_supp.pdf', width=20, height=20)
layout <- '
AB
CD
EF
GH
'
names(plots) <- LETTERS[1:7]
wrap_plots(plots, design = layout)
dev.off()

