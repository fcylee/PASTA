
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Find polyA sites that are differentially used across cells, utilizing
#' polyA residuals.
#' @param object An object
#' @param assay name of polyAsite assay to test. Default is polyA.
#' @param ident.1 Identity class to find polyA sites for
#' @param ident.2 A second identity class for comparison.
#' @param features Features to test. Default is all features. If a subset of
#' features is provided, pct.1 and pct.2 will be calculated using all features
#' with the same gene annotation where residuals have been calculated.
#' @param covariates Vector of covariates to adjust for (e.g. "sex"). Do NOT put
#' the group.by.sample column here.
#' @param group.by.sample Optional metadata column identifying which sample /
#' animal each cell came from (e.g. "sample_id"). Cells from the same sample are
#' not independent observations, so treating every cell as independent
#' overstates confidence (pseudo-replication). Setting this fits a mixed model
#' (via lmerTest) that accounts for that grouping, so the test reflects the
#' number of samples rather than the number of cells. When NULL (default) an
#' ordinary lm treating each cell as independent is used. For a given comparison
#' with too few samples or too few cells per sample to estimate the grouping
#' (see min.groups / min.cells.per.group), it falls back to lm automatically.
#' @param min.groups Minimum number of samples needed to fit the mixed model;
#' below this it falls back to lm. Default 3.
#' @param min.cells.per.group Minimum mean cells per sample needed to fit the
#' mixed model; below this it falls back to lm. Default 3.
#' @param gene.names Column name providing gene annotation of each polyA site.
#' Default is "Gene_Symbol"
#'
#' @importFrom stats lm relevel as.formula
#'
#' @rdname FindDifferentialPolyA
#' @concept differential_polyA
#' @export
#'
FindDifferentialPolyA <- function(
    object,
    assay = "polyA",
    ident.1,
    ident.2,
    features = NULL,
    covariates = NULL,
    group.by.sample = NULL,
    min.groups = 3,
    min.cells.per.group = 3,
    gene.names = "Gene_Symbol") {

  if( !inherits(object[[assay]], "polyAsiteAssay")){
    stop(paste0(assay," assay is not a polyAsiteAssay"))
  }

  if (dim(LayerData(object, layer="scale.data", assay = assay))[1] == 0)  {
    stop ("No features found in scale.data layer for specified assay. Run CalcPolyAResiduals prior to FindPolyASites")
  }

  if (!(gene.names %in% colnames(object[[assay]][[]]))) {
    stop("Gene.names column not found in meta.features, please make sure
         you are specific gene.names correctly")
  }

  if (!(ident.1 %in% unique(Idents(object)))) {
    stop("ident.1 not found in object")
  }
  df <- data.frame(ident = Idents(object))

  if (!is.null(covariates)) {
    for (i in 1:length(covariates)) {
      if (is.na(match(covariates[[i]], colnames(object[[]])))) {
        stop("Covariates are not found in meta data, please make sure you are specifying correctly.")
      }
      df[,i+1] <- object[[]][,match(covariates[[i]], colnames(object[[]]))]
    }
    colnames(df) <- c("ident", covariates)

  }

  if (!is.null(group.by.sample)) {
    if (is.na(match(group.by.sample, colnames(object[[]])))) {
      stop("group.by.sample '", group.by.sample, "' not found in meta data.")
    }
    df[[group.by.sample]] <- as.factor(
      object[[]][, match(group.by.sample, colnames(object[[]]))])
  }

  features <- features %||% rownames(x = object[[assay]]@scale.data)
  r.matrix <- object[[assay]]@scale.data

  # only test features that actually have residuals in scale.data
  missing.features <- setdiff(features, rownames(r.matrix))
  if (length(missing.features) > 0) {
    warning(length(missing.features), " of ", length(features),
            " requested features have no residuals in scale.data and will be skipped ",
            "(run CalcPolyAResiduals on them first). e.g. ",
            paste(head(missing.features, 3), collapse = ", "))
    features <- intersect(features, rownames(r.matrix))
  }
  if (length(features) == 0) {
    stop("None of the requested features have residuals in scale.data.")
  }

  df$ident <- relevel(df$ident, ref = ident.2)

  sub <- subset(df, ident %in% c(ident.1, ident.2))

  missing.cells <- setdiff(rownames(sub), colnames(r.matrix))
  if (length(missing.cells) > 0) {
    warning(length(missing.cells), " of ", nrow(sub),
            " cells in ident.1/ident.2 have no residuals in scale.data and will be ",
            "dropped from the test (they were likely filtered during CalcPolyAResiduals).")
    sub <- sub[setdiff(rownames(sub), missing.cells), , drop = FALSE]
  }
  if (nrow(sub) == 0) {
    warning("No cells with residuals remain for ", ident.1, " vs ", ident.2,
            "; returning NULL and skipping this comparison.")
    return(NULL)
  }

  # drop factor levels no longer present after subsetting/cell filtering
  sub <- droplevels(sub)

  # both groups must still be represented, or the ident contrast is undefined
  if (nlevels(sub$ident) < 2) {
    warning("Only one of ", ident.1, " / ", ident.2,
            " has cells with residuals; returning NULL and skipping this comparison.")
    return(NULL)
  }

  # covariates that are constant in this subset break lm() contrasts and add
  # nothing to the model; drop them (with a note) so the comparison still runs
  if (!is.null(covariates)) {
    const.cov <- covariates[vapply(covariates,
                                   function(cv) length(unique(sub[[cv]])) < 2,
                                   logical(1))]
    if (length(const.cov) > 0) {
      warning("Dropping covariate(s) constant in ", ident.1, " vs ", ident.2,
              ": ", paste(const.cov, collapse = ", "))
      sub <- sub[, setdiff(colnames(sub), const.cov), drop = FALSE]
    }
  }

  r.matrix.sub <- r.matrix[features, rownames(sub), drop = FALSE]

  # decide once per comparison whether to fit a mixed model: the sample grouping
  # must have enough samples, and enough cells per sample, to be estimable.
  # Otherwise (or if lmerTest is unavailable) fall back to a plain per-cell lm.
  fixed.terms <- setdiff(colnames(sub), c("residuals", group.by.sample))
  use.lmer <- FALSE
  if (!is.null(group.by.sample)) {
    n.groups <- length(unique(sub[[group.by.sample]]))
    cells.per.group <- nrow(sub) / n.groups
    if (!requireNamespace("lmerTest", quietly = TRUE)) {
      warning("group.by.sample set but lmerTest is not installed; ",
              "falling back to fixed-effects lm.")
    } else if (n.groups < min.groups || cells.per.group < min.cells.per.group) {
      warning("group.by.sample '", group.by.sample, "' has ", n.groups, " group(s), ",
              round(cells.per.group, 1), " cells/group for ", ident.1, " vs ",
              ident.2, "; too few to estimate - falling back to fixed-effects lm.")
    } else {
      use.lmer <- TRUE
    }
  }

  rhs.fixed <- paste(fixed.terms, collapse = " + ")
  form.fixed <- as.formula(paste("residuals ~", rhs.fixed))
  form.mixed <- if (use.lmer) {
    as.formula(paste("residuals ~", rhs.fixed, "+ (1 | ", group.by.sample, ")"))
  } else NULL

  model.used <- if (use.lmer) "lmer" else "lm"

  # Trial fit: some comparisons are so thin that lmer cannot be estimated at all
  # (non-positive-definite vcov, which makes lmerTest error on every peak). Probe
  # once here with the highest-variance peak - the strongest-signal peak, i.e.
  # the best case for estimating a positive sample variance. If even that peak
  # fails, weaker peaks certainly will, so run the whole comparison as lm and
  # flag it "lm_fallback". If it succeeds, keep lmer and let the per-peak net
  # skip any individual stragglers.
  if (use.lmer) {
    row.var <- apply(r.matrix.sub, 1, stats::var)
    trial.i <- if (all(is.na(row.var) | row.var == 0)) 1L else which.max(row.var)
    trial.sub <- sub
    trial.sub$residuals <- as.numeric(r.matrix.sub[trial.i, ])
    trial.ok <- tryCatch(
      suppressWarnings(suppressMessages({
        co <- summary(lmerTest::lmer(form.mixed, data = trial.sub))$coefficients
        all(c("Estimate", "Std. Error") %in% colnames(co))
      })),
      error = function(e) FALSE
    )
    if (!isTRUE(trial.ok)) {
      warning("lmer not estimable for ", ident.1, " vs ", ident.2,
              " (e.g. non-positive-definite vcov); running lm for all peaks ",
              "(model = 'lm_fallback').")
      use.lmer <- FALSE
      model.used <- "lm_fallback"
    }
  }

  form <- if (use.lmer) form.mixed else form.fixed

  # provenance carried on every result row (constant within a comparison)
  n.samples.used <- if (is.null(group.by.sample)) NA_integer_
                    else length(unique(sub[[group.by.sample]]))

  all.models <- lapply(
    X = 1:nrow(x = r.matrix.sub),
    FUN = function(i) {
      sub$residuals <- as.numeric(r.matrix.sub[i, ])
      peak.name <- rownames(r.matrix.sub)[i]

      # a single peak can fail (e.g. lmer singular / non-positive-definite
      # vcov); catch it so one bad peak skips itself rather than aborting the
      # whole comparison. warnings are logged and muffled as before. Also record
      # whether an lmer fit was singular (variance ~0) - if so the sample
      # grouping wasn't really estimated and the p-value reverts toward the
      # over-confident per-cell result, so flag it rather than trust it blindly.
      fit <- tryCatch(
        withCallingHandlers(
          {
            model <- if (use.lmer) {
              # suppressMessages() silences lme4's un-labelled "boundary
              # (singular) fit" note; we capture singularity in a column below
              suppressMessages(lmerTest::lmer(form, data = sub))
            } else {
              lm(form, data = sub)
            }
            list(co        = summary(model)$coefficients,
                 singular  = if (use.lmer) lme4::isSingular(model) else NA,
                 converged = if (use.lmer)
                   is.null(model@optinfo$conv$lme4$messages) else NA)
          },
          warning = function(w) {
            message(sprintf("[%s vs %s | %s] %s",
                            ident.1, ident.2, peak.name, conditionMessage(w)))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          message(sprintf("[%s vs %s | %s] fit failed, skipping peak: %s",
                          ident.1, ident.2, peak.name, conditionMessage(e)))
          NULL
        }
      )

      # skip the peak if the fit failed or the coefficient table is unusable
      co <- fit$co
      if (is.null(co) || !all(c("Estimate", "Std. Error") %in% colnames(co))) {
        return(NULL)
      }

      # index columns by name so this works for both lm (Estimate/Std. Error/
      # t value/Pr(>|t|)) and lmerTest (which inserts an extra df column);
      # tolerate a missing p-value column (can happen on degenerate vcov)
      tcol <- if ("t value"  %in% colnames(co)) co[, "t value"]  else NA_real_
      pcol <- if ("Pr(>|t|)" %in% colnames(co)) co[, "Pr(>|t|)"] else NA_real_
      data.frame(
        Estimate    = co[, "Estimate"],
        std_error   = co[, "Std. Error"],
        t           = tcol,
        p.value     = pcol,
        coefficient = rownames(co),
        peak        = peak.name,
        model       = model.used,
        n_samples   = n.samples.used,
        is_singular = fit$singular,
        converged   = fit$converged,
        stringsAsFactors = FALSE
      )
    }
  )
  results <- do.call(rbind, all.models)
  if (is.null(results) || nrow(results) == 0) {
    warning("No peaks could be fit for ", ident.1, " vs ", ident.2,
            "; returning NULL and skipping this comparison.")
    return(NULL)
  }
  main.effects <- results[grep("ident", results$coefficient),]

  gene.idx = match(gene.names, colnames(object[[assay]][[]]))
  main.effects$symbol <- object[[assay]][[]][main.effects$peak,gene.idx]

  #get all features within genes to calculate percentage.usage
  #also require that they have residuals calculated
  features.all.genes <- rownames(object[[assay]][[]][object[[assay]][[]][[gene.names]] %in% unique(main.effects$symbol), ])
  features.all.genes <- intersect(features.all.genes, rownames(LayerData(object[[assay]], layer="scale.data")))
  #calculate percentage usage in group 1 and group 2
  percent.1 <- percentage.usage(object,
                                assay = assay,
                                cells = WhichCells(object, idents = ident.1),
                                features = features.all.genes,
                                gene.names = gene.names)
  main.effects$percent.1 <- percent.1[main.effects$peak]
  percent.2 <- percentage.usage(object,
                                assay = assay,
                                cells = WhichCells(object, idents = ident.2),
                                features = features.all.genes,
                                gene.names = gene.names)
  main.effects$percent.2 <- percent.2[main.effects$peak]
  main.effects$p_val_adj <- main.effects$p.value * nrow(object[[assay]]@scale.data)
  main.effects$p_val_adj[main.effects$p_val_adj  > 1] <- 1

  rownames(main.effects) <- main.effects$peak
  main.effects.return <- main.effects[,c("Estimate", "std_error", "p.value", "p_val_adj", "percent.1", "percent.2", "symbol", "model", "n_samples", "is_singular", "converged")]

  #order by p-value
  main.effects.return <- main.effects.return[ with(main.effects.return, order(p_val_adj, -Estimate)),]
  return(main.effects.return)
}


percentage.usage <- function( object,
                              assay = "polyA",
                              cells,
                              features,
                              gene.names = "Gene_Symbol") {
  df <- data.frame(peak=features)
  meta <- object[[assay]]@meta.features
  df$symbol <- meta[features,gene.names]
  df$counts1 <- rowSums(object[[assay]]@counts[df$peak, cells, drop = FALSE])
  sum1 <- aggregate(df$counts1, by=list(gene=df$symbol), FUN=sum)
  colnames(sum1) <- c("symbol", "sum")
  df <- merge(df, sum1, by="symbol")
  df$frac <- df$counts1/df$sum
  rownames(df) <- df$peak
  df <- df[features,]
  df$frac[df$sum==0] <- 0
  v <- as.vector(df$frac)
  names(v) <- features
  return(v)
}



