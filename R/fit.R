
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Calculate PolyA Residuals
#'
#'
#' @param object Seurat object containing a polyAsiteAssay
#' @param assay Name of polyAsiteAssay to be used in calculating polyAresiduals
#' @param features Features to include in calculation of polyA residuals.
#' Default is to use all features.
#' @param background Identity of cells to use as background.
#' Default is to use all cells as a background.
#' @param gene.names Column containing the gene where each polyA site is annotated.
#' Default is symbol.
#' @param min.counts.background Features with at least this many counts in the background cells are included in calculation
#' @param min.variance Sets minimum variance. Default is 0.1.
#' @param sample.n Max number of observations to sample in each bin when performing regularization. Default is 1000.
#' @param mc.cores Number of cores used to fit the Dirichlet-multinomial across
#' genes in parallel via \code{\link[parallel]{mclapply}}. Default is 1 (serial).
#' Forking is not available on Windows, where values > 1 fall back to serial.
#' @param do.center Return the centered residuals. Default is TRUE.
#' @param do.scale Return the scaled residuals. Default is TRUE.
#' @param residuals.max Clip residuals above this value. Default is NULL (no clipping).
#' @param residuals.min Clip residuals below this value. Default is NULL (no clipping).
#' @param number.bins Number of bins used to perform regularization. Default is 30.
#' @param chunk.size Number of cells processed per block when computing residuals.
#' Residuals are streamed in column blocks so the large per-cell intermediates are
#' never all held at once; smaller values lower peak memory. Default is 2000.
#' @param verbose Print messages.
#'
#'
#' @return Returns a Seurat object with polyAresiduals assay
#'
#' @importFrom parallel mclapply
#' @importFrom Matrix sparseMatrix
#' @export
#' @concept residuals
#'
CalcPolyAResiduals <- function(object,
                               assay="polyA",
                               features=NULL,
                               background = NULL,
                               gene.names = "Gene_Symbol",
                               min.counts.background = 5,
                               min.variance = 0.1,
                               sample.n = 1000,
                               mc.cores = 1L,
                               do.scale = FALSE,
                               do.center = FALSE,
                               residuals.max = NULL,
                               residuals.min = NULL,
                               number.bins = 30,
                               chunk.size = 2000,
                               verbose=TRUE)
 {
  if(verbose) {
    message("Calculating background distribution")
  }

  #if features in NULL, then specify all features in polyA assay
  if (is.null(features)) {
    features <- rownames(LayerData(object, assay=assay, layer="counts"))
  }

  #if background is NULl, then make a dummy variable
  if (is.null(background)) {
    message("Using all cells in order to estimate background distribution")
    object$dummy <- "all"
    Idents(object) <- object$dummy
    background.use = "all"
  }  else {
    background.use = background
    if (!(background %in% unique(Idents(object)))) {
      stop("background must be one of the Idents of seurat object")
    }
    message(paste0("Using ", background, " as background distribution"))
  }


  #check if symbols are contained in meta features
  if (!(gene.names %in% colnames(object[[assay]]@meta.features))) {
    stop("Gene.names column not found in meta.features, please make sure
         you are specific gene.names correctly")
  }

  if (sum(is.na(object[[assay]]@meta.features[features,gene.names])) > 0) {
    features.no.anno <- features[is.na(object[[assay]]@meta.features[features,gene.names])]
    message(paste0("Removing ", length(features.no.anno), " sites without a gene annotation"))
    features <- setdiff(features, features.no.anno)
  }

  ##############################################################################
  #get pseudobulked fraction of reads from background
  background.dist <- GetBackgroundDist(object = object, features = features,
                                       background = background.use,
                                       gene.names = gene.names,
                                       assay = assay,
                                       min.counts.background = min.counts.background)
  #remove features without gene annotation
  background.dist <- subset(background.dist, gene!="_-")
  background.dist <- subset(background.dist, gene!="_+")
  features.use <- background.dist$peak

  ##############################################################################
  #calculate sum of counts within each gene
  m <- LayerData(object = object, assay=assay, layer="counts")
  m <- m[background.dist$peak,]
  #m <- m[order(match(rownames(m), background.dist$peak)), ]
  # sum counts within each gene using a sparse gene-by-peak indicator matrix.
  # This keeps m sparse instead of densifying it to a peaks x cells matrix as
  # rowsum() would.
  gene.factor <- factor(background.dist$gene)
  gene.indicator <- sparseMatrix(i = as.integer(gene.factor),
                                 j = seq_along(gene.factor),
                                 x = 1,
                                 dims = c(nlevels(gene.factor), nrow(m)))
  gene.sum <- as.matrix(gene.indicator %*% m)
  rownames(gene.sum) <- levels(gene.factor)
  genes <- rownames(gene.sum)

  ##############################################################################
  #fit dirichlet multinomial distribution
  if(verbose) {
    message("Running Dirichlet Multionmial Regression")
  }

  ncells = dim(object)[2]
  background.cells <- WhichCells(object, idents=background)
  # keep the background matrix sparse; DirichletMultionmial densifies only the
  # small per-gene submatrix it needs, rather than holding a full dense copy
  m.background <- m[, background.cells, drop = FALSE]

  # precompute the peaks belonging to each gene once, so the per-gene fit does
  # not rescan every rowname (avoids O(genes x peaks) lookups)
  peaks.by.gene <- split(background.dist$peak, background.dist$gene)

  # fit dirichlet multinomial for each gene. Genes are independent, so this is
  # parallelized across cores when mc.cores > 1. Each fit returns only its
  # parameters; expected counts and variances are generated on demand below.
  fits <- parallel::mclapply(genes, DirichletMultionmial,
                             peaks.by.gene = peaks.by.gene,
                             m.background = m.background, mc.cores = mc.cores)

  successful <- !vapply(fits, is.null, logical(1))
  if (!any(successful)) {
    stop("No genes successfully fitted with Dirichlet-multinomial model. ",
         "Try reducing min.counts.background or using fewer features.")
  }
  if (verbose && sum(successful) < length(genes)) {
    failed_count <- length(genes) - sum(successful)
    message(paste0("Warning: ", failed_count, " out of ", length(genes),
                   " genes failed Dirichlet-multinomial fitting"))
    message(paste0("Proceeding with ", sum(successful), " successful genes"))
  }
  fits <- fits[successful]

  # Per-peak parameter vectors (aligned to peaks.kept). Gene totals are looked up
  # by name via gene.idx, so the calculation does not depend on background.dist
  # being ordered the same way as the fitted genes.
  peaks.kept <- unlist(lapply(fits, `[[`, "peaks"), use.names = FALSE)
  param.peak <- unlist(lapply(fits, `[[`, "param"), use.names = FALSE)
  sump.peak  <- unlist(lapply(fits, function(z) rep(z$sum.p, length(z$peaks))), use.names = FALSE)
  pfrac      <- param.peak / sump.peak
  gene.of.peak <- background.dist$gene[match(peaks.kept, background.dist$peak)]
  gene.idx     <- match(gene.of.peak, rownames(gene.sum))

  ##############################################################################
  ### regularize dirichlet multinomial variance
  if(verbose) {
    message("Regularizing Dirichlet Multionmial Variance")
  }

  # Build the variance-regularization grid once, from the background cells only.
  grid <- BuildVarGrid(param.peak = param.peak, sump.peak = sump.peak, pfrac = pfrac,
                       gene.idx = gene.idx, gene.sum = gene.sum,
                       background.cells = background.cells,
                       number.bins = number.bins, sample.n = sample.n)

  ##############################################################################
  # Compute residuals in blocks of cells so the peaks x cells intermediates
  # (expected counts, variance, regularized variance) never all coexist.
  residual.matrix <- matrix(0, nrow = length(peaks.kept), ncol = ncells,
                            dimnames = list(peaks.kept, colnames(object)))
  for (start in seq(1, ncells, by = chunk.size)) {
    cols <- start:min(start + chunk.size - 1L, ncells)
    n.block  <- gene.sum[gene.idx, cols, drop = FALSE]
    ec.block <- (n.block * param.peak) / sump.peak
    v.block  <- ec.block * (1 - pfrac) * (n.block + sump.peak) / (1 + sump.peak)
    obs.key  <- findInterval(as.vector(ec.block), grid$ec_grid) +
                findInterval(as.vector(n.block), grid$n_grid) * grid$key.base
    reg.var  <- grid$grid.var[match(obs.key, grid$grid.key)]
    na.var   <- is.na(reg.var)
    reg.var[na.var] <- as.vector(v.block)[na.var]
    reg.var[reg.var < min.variance] <- min.variance
    counts.block <- as.matrix(m[peaks.kept, cols, drop = FALSE])
    residual.matrix[, cols] <- (counts.block - ec.block) / sqrt(matrix(reg.var, nrow = length(peaks.kept)))
  }

  #change to same order as counts slot
  features.order <- rownames(residual.matrix)[order(
    match(
      rownames(residual.matrix),
      rownames(LayerData(object = object, assay=assay, layer="counts"))
    )
  )]
  residual.matrix <- residual.matrix[features.order,]
  # scale() with center = FALSE and scale = FALSE just copies the matrix, so only
  # call it when centering or scaling is actually requested
  if (do.center || do.scale) {
    residual.matrix <- scale(residual.matrix, center=do.center, scale= do.scale )
  }
  if (!is.null(residuals.max)) {
    residual.matrix[residual.matrix > residuals.max] <- residuals.max
  }

  if (!is.null(residuals.min)) {
    residual.matrix[residual.matrix < residuals.min] <- residuals.min
  }
  #change default assay
  DefaultAssay(object = object) <- assay
  #need to met SetAssayData, GetAssayData for residuals
  LayerData(object, layer = "scale.data") <- residual.matrix
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' Get Background Distribution
#'
#' Calculated Pseudobulk Ratios of Each Isoform within a gene for background distribution.
#' @param object Seurat object containing a polyAsiteAssay
#' @param assay Name of polyAsiteAssay to be used in calculating polyAresiduals
#' @param features Features to include in calculation of polyA residuals.
#' If NULL, use all features.
#' @param background Identity of cells to use as background
#' If NULL, uses all cells combined as a background
#' @param gene.names Name of column containing gene annotations
#' @param min.counts.background Features with at least this many counts in the background cells are included in calculation
#'
#' @return Returns a data frame containing all peaks within genes that have multiple polyA sites that meet min.counts.background criteria
#'
#' @importFrom stats aggregate
#' @concept residuals
#'
GetBackgroundDist <- function(object, features, background, gene.names, assay,  min.counts.background) {
  # returns the pseudobulked background distribution for peaks specified
  # must contain gene information in meta data

  suppressMessages(nt.pseudo <- AverageExpression(object, features = features, assays = assay, slot="counts"))
  nt.pseudo <- data.frame(background = nt.pseudo[[1]][,background]) #subset just the background
  nt.pseudo$background <- nt.pseudo$background * sum(Idents(object)==background)
  nt.pseudo$gene <- paste0(object[[assay]]@meta.features[features, gene.names], "_", object[[assay]]@meta.features[features, "strand"])
  nt.pseudo$peak <- rownames(nt.pseudo)

  nt.pseudo <- nt.pseudo[nt.pseudo$background>min.counts.background,] #subset to peaks with min number of counts
  genes.use <- nt.pseudo$gene[duplicated(nt.pseudo$gene)] #only use genes with at least 2 peaks per gene
  nt.pseudo <- nt.pseudo[nt.pseudo$gene %in% genes.use,]

  if ( length(genes.use)  ==  0) {
    stop("Found no genes with more than 2 features within a gene. Please make sure you are including
         all peaks within a gene you would like to include.")
  }

  tmp <- aggregate(nt.pseudo$background, list(nt.pseudo$gene), FUN=sum)
  colnames(tmp) <- c("gene", "sum")
  nt.pseudo <- merge(nt.pseudo, tmp, by="gene")
  nt.pseudo$frac <- nt.pseudo$background/ nt.pseudo$sum
  return(nt.pseudo)
}


#' Fit Dirichlet Multinomial Distribution
#'
#' Fit a Dirichlet-multinomial distribution on the peaks within a gene using the
#' background cells, and return the fitted parameters. Expected counts and variances
#' are generated later, on demand, from these parameters.
#'
#' @param gene.test which gene to use
#' @param peaks.by.gene named list mapping each gene to its peaks
#' @param m.background sparse counts matrix restricted to the background cells
#'
#' @return A list with the gene's \code{peaks}, the fitted concentration parameters
#' \code{param}, and their sum \code{sum.p}. Returns NULL if the fit fails or the
#' optimizer drops a category.
#'
#' @importFrom MGLM MGLMfit
#' @concept residuals
#'
DirichletMultionmial <- function(
  gene.test,
  peaks.by.gene,
  m.background
) {
  peaks <- peaks.by.gene[[gene.test]]
  # densify only this gene's (few peaks) x (background cells) submatrix
  t <- as.matrix(m.background[peaks, , drop = FALSE])
  t <- t(t)
  # drop background cells with no counts for this gene: a zero-total observation
  # has Dirichlet-multinomial likelihood 1 and contributes nothing to the fit,
  # so removing it leaves the estimate unchanged while shrinking the optimizer input
  t <- t[rowSums(t) > 0, , drop = FALSE]
  fit <- try(compareFit <- suppressWarnings(MGLMfit(t, dist="DM")), silent=TRUE)

  if (!inherits(fit, "try-error")) {
    param <- compareFit@estimate
    # MGLMfit can drop an all-zero category, returning fewer parameters than
    # peaks; that gene cannot be mapped back, so skip it (treated as a failed fit)
    if (length(param) != length(peaks)) {
      return(NULL)
    }
    return(list(peaks = peaks, param = param, sum.p = sum(param)))
  }
}


#' Build the Dirichlet-multinomial variance-regularization grid
#'
#' Computes expected counts and Dirichlet-multinomial variances for the background
#' cells, then fits a smooth variance surface over a regular
#' (gene-total, expected-count) grid via kernel regression. The returned grid is
#' used to regularize the variance of every cell when residuals are computed.
#'
#' @param param.peak fitted concentration parameter for each (kept) peak
#' @param sump.peak sum of concentration parameters for each peak's gene
#' @param pfrac \code{param.peak / sump.peak}, the background fraction for each peak
#' @param gene.idx row of \code{gene.sum} giving each peak's gene total
#' @param gene.sum genes x cells matrix of per-gene counts
#' @param background.cells cells used to estimate the background distribution
#' @param number.bins number of bins along each grid axis
#' @param sample.n max observations sampled per bin before kernel regression
#'
#' @return A list with the bin edges (\code{ec_grid}, \code{n_grid}), the regularized
#' variance at each grid cell (\code{grid.var}), the integer keys of those grid cells
#' (\code{grid.key}), and the key base (\code{key.base}) used to encode (ec, n) bins.
#'
#' @importFrom stats quantile
#' @importFrom gplm kreg
#' @concept residuals
#'
BuildVarGrid <- function(param.peak, sump.peak, pfrac, gene.idx, gene.sum,
                         background.cells, number.bins, sample.n) {
  # expected counts and variances for the background cells (peaks x background cells)
  n.bg  <- gene.sum[gene.idx, background.cells, drop = FALSE]
  ec.bg <- (n.bg * param.peak) / sump.peak
  v.bg  <- ec.bg * (1 - pfrac) * (n.bg + sump.peak) / (1 + sump.peak)

  ec.v <- as.vector(ec.bg); n.v <- as.vector(n.bg); md.var <- as.vector(v.bg)
  keep <- n.v > 0
  ec.v <- ec.v[keep]; n.v <- n.v[keep]; md.var <- md.var[keep]

  cutoff.ec <- quantile(ec.v, 0.99)
  cutoff.n  <- quantile(n.v, 0.99)
  max.n  <- max(n.v[n.v < cutoff.n]);   max.ec <- max(ec.v[ec.v < cutoff.ec])
  min.n  <- min(n.v[n.v < cutoff.n]);   min.ec <- min(ec.v[ec.v < cutoff.ec])

  lx <- ly <- number.bins
  n_grid  <- min.n  + (max.n  - min.n)  / lx * 0:lx
  ec_grid <- min.ec + (max.ec - min.ec) / ly * 0:ly

  # subsample within each (ec.bin, n.bin) cell before the kernel regression
  ec_n <- paste0(findInterval(ec.v, ec_grid), "_", findInterval(n.v, n_grid))
  sub  <- unlist(lapply(split(seq_along(ec_n), ec_n), sample_within_groups, sample.n = sample.n))

  n_grid_midpoints  <- calculate_midpoints(min.n,  max.n,  lx)
  ec_grid_midpoints <- calculate_midpoints(min.ec, max.ec, ly)
  grid <- matrix(cbind(rep(n_grid_midpoints, length(ec_grid_midpoints)),
                       rep(ec_grid_midpoints, each = length(n_grid_midpoints))), ncol = 2)
  mh <- kreg(x = matrix(cbind(n.v[sub], ec.v[sub]), ncol = 2), y = md.var[sub], grid = grid)

  # integer key per grid cell; matches the (ec.bin, n.bin) encoding used at lookup.
  # key.base exceeds any ec bin index so the encoding is a bijection.
  key.base <- length(ec_grid) + 1L
  ec_bin <- rep(1:length(ec_grid_midpoints), length(n_grid_midpoints))
  n_bin  <- rep(1:length(n_grid_midpoints), each = length(ec_grid_midpoints))
  list(ec_grid = ec_grid, n_grid = n_grid, grid.var = mh$y,
       grid.key = ec_bin + n_bin * key.base, key.base = key.base)
}


#generate midpoints
calculate_midpoints <- function(a, b, n_intervals) {
  width <- (b - a) / n_intervals
  first_point <- a + width/2
  midpoints <- seq(first_point, by = width, length.out = n_intervals)
  return(midpoints)
}


sample_within_groups <- function(x, sample.n) {
  if (length(x) <=  sample.n) return(x)
  x[x %in% sample(x, sample.n)]
}

