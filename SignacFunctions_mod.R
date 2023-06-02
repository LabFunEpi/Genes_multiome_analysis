SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

GetGroups <- function(
  object,
  group.by,
  idents
) {
  if (is.null(x = group.by)) {
    obj.groups <- Idents(object = object)
  } else {
    obj.md <- object[[group.by]]
    obj.groups <- obj.md[, 1]
    names(obj.groups) <- rownames(x = obj.md)
  }
  if (!is.null(idents)) {
    obj.groups <- obj.groups[obj.groups %in% idents]
  }
  return(obj.groups)
}

ApplyMatrixByGroup <- function(
  mat,
  groups,
  fun,
  normalize = TRUE,
  group.scale.factors = NULL,
  scale.factor = NULL
) {
  if (normalize) {
    if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
      stop("If normalizing counts, supply group scale factors")
    }
  }
  all.groups <- as.character(x = unique(x = groups))
  if (any(is.na(x = groups))) {
    all.groups <- c(all.groups, NA)
  }
  ngroup <- length(x = all.groups)
  npos <- ncol(x = mat)

  group <- unlist(
    x = lapply(X = all.groups, FUN = function(x) rep(x, npos))
  )
  position <- rep(x = as.numeric(x = colnames(x = mat)), ngroup)
  count <- vector(mode = "numeric", length = npos * ngroup)

  for (i in seq_along(along.with = all.groups)) {
    grp <- all.groups[[i]]
    if (is.na(x = grp)) {
      pos.cells <- names(x = groups)[is.na(x = groups)]
    } else {
      pos.cells <- names(x = groups)[groups == all.groups[[i]]]
    }
    if (length(x = pos.cells) > 1) {
      totals <- fun(x = mat[pos.cells, ])
    } else {
      totals <- mat[pos.cells, ]
    }
    count[((i - 1) * npos + 1):((i * npos))] <- totals
  }

  # construct dataframe
  coverages <- data.frame(
    "group" = group, "position" = position, "count" = count,
    stringsAsFactors = FALSE
  )

  if (normalize) {
    scale.factor <- SetIfNull(
      x = scale.factor, y = median(x = group.scale.factors)
    )
    coverages$norm.value <- coverages$count /
      group.scale.factors[coverages$group] * scale.factor
  } else {
    coverages$norm.value <- coverages$count
  }
  return(coverages)
}

TSSPlotdata <- function(
  object,
  group.by = NULL,
  idents = NULL
) {
  assay <- DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
  # get the normalized TSS enrichment matrix
  positionEnrichment <- GetAssayData(
    object = object,
    assay = assay,
    slot = "positionEnrichment"
  )
  if (!("TSS" %in% names(x = positionEnrichment))) {
    stop("Position enrichment matrix not present in assay")
  }
  enrichment.matrix <- positionEnrichment[["TSS"]]

  # remove motif and expected
  if (nrow(x = enrichment.matrix) == (ncol(x = object) + 2)) {
    enrichment.matrix <- enrichment.matrix[1:(nrow(x = enrichment.matrix) - 2), ]
  }

  # average the signal per group per base
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  groupmeans <- ApplyMatrixByGroup(
    mat = enrichment.matrix,
    groups = obj.groups,
    fun = colMeans,
    normalize = FALSE
  )
  return(groupmeans)
}

TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  if (record.ident) {
    nrep <- elementNROWS(x = reads)
  }
  reads <- unlist(x = reads, use.names = FALSE)
  if (length(x = reads) == 0) {
    df <- data.frame(
      "chr" = "",
      "start" = "",
      "end" = "",
      "cell" = "",
      "count" = ""
    )
    df <- df[-1, ]
    return(df)
  }
  reads <- stringi::stri_split_fixed(str = reads, pattern = "\t")
  n <- length(x = reads[[1]])
  unlisted <- unlist(x = reads)
  e1 <- unlisted[n * (seq_along(along.with = reads)) - (n - 1)]
  e2 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 2)])
  e3 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 3)])
  e4 <- unlisted[n * (seq_along(along.with = reads)) - (n - 4)]
  e5 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 5)])
  df <- data.frame(
    "chr" = e1,
    "start" = e2,
    "end" = e3,
    "cell" = e4,
    "count" = e5,
    stringsAsFactors = FALSE,
    check.rows = FALSE,
    check.names = FALSE
  )
  if (record.ident) {
    df$ident <- rep(x = seq_along(along.with = nrep), nrep)
  }
  return(df)
}

GetReadsInRegion <- function(
  cellmap,
  region,
  tabix.file,
  cells = NULL,
  verbose = TRUE,
  ...
) {
  file.to.object <- names(x = cellmap)
  names(x = file.to.object) <- cellmap

  if (verbose) {
    message("Extracting reads in requested region")
  }
  if (!is(object = region, class2 = "GRanges")) {
    region <- StringToGRanges(regions = region, ...)
  }
  # remove regions that aren't in the fragment file
  common.seqlevels <- intersect(
    x = seqlevels(x = region),
    y = Rsamtools::seqnamesTabix(file = tabix.file)
  )
  region <- keepSeqlevels(
    x = region,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  reads <- Rsamtools::scanTabix(file = tabix.file, param = region)
  reads <- TabixOutputToDataFrame(reads = reads)
  reads <- reads[
    fmatch(x = reads$cell, table = cellmap, nomatch = 0L) > 0,
  ]
  # convert cell names to match names in object
  reads$cell <- file.to.object[reads$cell]
  if (!is.null(x = cells)) {
    reads <- reads[reads$cell %in% cells, ]
  }
  if (nrow(reads) == 0) {
    return(reads)
  }
  reads$length <- reads$end - reads$start
  return(reads)
}
MultiGetReadsInRegion <- function(
  object,
  region,
  fragment.list = NULL,
  assay = NULL,
  ...
) {
  if (inherits(x = object, what = "Seurat")) {
    # pull the assay
    assay <- DefaultAssay(object = object)
    object <- object[[assay]]
  }
  fragment.list <- Fragments(object = object)
  if (length(x = fragment.list) == 0) {
    # no fragments set
    stop("No fragment files found")
  }
  res <- data.frame()
  for (i in seq_along(along.with = fragment.list)) {
    tbx.path <- GetFragmentData(object = fragment.list[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragment.list[[i]], slot = "cells")
    tabix.file <- Rsamtools::TabixFile(file = tbx.path)
    open(con = tabix.file)
    reads <- GetReadsInRegion(
      cellmap = cellmap,
      region = region,
      tabix.file = tabix.file,
      ...
    )
    res <- rbind(res, reads)
    close(con = tabix.file)
  }
  return(res)
}

FragmentHistogramdata <- function(
  object,
  region = "chr1-1-2000000",
  group.by = NULL,
  log.scale = FALSE,
  ...
) {
  cells <- colnames(x = object)
  assay <- DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
  reads <- MultiGetReadsInRegion(
    object = object,
    assay = assay,
    region = region,
    cells = cells,
    verbose = FALSE,
    ...
  )
  # add group information
  if (is.null(x = group.by)) {
    groups <- Idents(object = object)
  } else {
    md <- object[[]]
    groups <- object[[group.by]]
    groups <- groups[, 1]
    names(x = groups) <- rownames(x = md)
  }
  reads$group <- groups[reads$cell]
  return(reads)
}
