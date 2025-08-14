

#' Overlap PAS with curated polyA database
#' @param object Seurat object containing a polyAsiteAssay
#' @param assay Name of polyAsiteAssay
#' @param polyAdb.file Location to find curated polyA atlas file (tsv expected).
#' @param max.dist Keep sites within this distance to a curated polyA site in the atlas. Default is 50 nucleotides.
#' @param atlas Atlas versions. Supported for "polyAdbv3-hg38", "PolyASitev3.0", "PolyASitev2.0".
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges distanceToNearest trim
#' @importFrom plyranges anchor_3p mutate
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @importFrom methods slot
#' @importFrom utils read.table
#'
#' @return Seurat object with annotated meta.features for polyAsiteAssay. Keeps all features
#' in original matrix, but annotations columns will be populated with NAs if a feature does not
#' overlap with a curated polyA site within the distance specified my max.dist.
#' @concept annotation
#' @export
#'
GetPolyADbAnnotation <- function(
  object,
  assay = "polyA",
  polyAdb.file = NULL, #
  max.dist = 50,
  atlas = "polyAdbv3-hg38")
{
  #readin polyAdb file and make GRanges file
  if (!( file.exists(polyAdb.file))) {
    stop("Please check that you have specified the location of polyAdb.file correctly")
  }

  anno <- read.table(file = polyAdb.file, header = TRUE, sep = "\t", quote = "")

  if( !atlas %in% c("polyAdbv3-hg38", "PolyASitev3.0", "PolyASitev2.0") ){
    stop(paste0("Atlas formats supported are polyAdbv3-hg38, PolyASitev3.0 or PolyASitev2.0."))
  }
  if( atlas %in% c("polyAdbv3-hg38") ){
    GR.polyA.db = makeGRangesFromDataFrame( anno,
                                            keep.extra.columns = TRUE,
                                            seqnames.field = "hg38_Chromosome_format",
                                            start.field = "hg38_Position",
                                            end.field = "hg38_Position",
                                            strand.field = "Strand")  
  }
  if( atlas %in% c("PolyASitev3.0") ){
    if(min(anno$stringency_level <= 62) ){ 
      message(paste0("Important! Consider pre-filtering the atlas by stringency level. See polyAsite v3.0 paper Figure 1B."))  
    }
  }
  if( atlas %in% c("PolyASitev3.0", "PolyASitev2.0") ){
    GR.polyA.db = makeGRangesFromDataFrame( anno,
                                            keep.extra.columns = TRUE,
                                            seqnames.field = "chrom",
                                            start.field = "chromStart",
                                            end.field = "chromEnd",
                                            strand.field = "strand")  
  }

  if( !assay %in% Assays(object) ){
    stop(paste0(assay," assay is not present in object"))
  }
  if( !inherits(object[[assay]], "polyAsiteAssay")){
    stop(paste0(assay," assay is not a polyAsiteAssay"))
  }
  if( !assay == DefaultAssay(object)){
    DefaultAssay(object) <- assay
    message(paste0("Setting default assay to ",assay))
  }
  ranges = slot(object = object[[assay]], name = "ranges")
  if ( ! "strand" %in% colnames(object[[assay]]@meta.features)){
    #object[[assay]] <- AddMetaData(object = object[[assay]] , metadata = as.character(strand(ranges)) , col.name = "strand")
    object[[assay]]@meta.features$strand <- as.character(strand(ranges))
  }

  features <- rownames(object[[assay]]@counts)
  if (  "*" %in% unique(strand(ranges)) ){
      warning("\n Cannot annotate unstranded PAS.\n")  }
  if( !all(seqlevelsStyle(ranges) == seqlevelsStyle(GR.polyA.db)) ) {
    if( seqlevelsStyle(ranges)[1] == "UCSC"){
      warning("\n Annotation does not match between ranges and polyAdb.\n Annotation set to USCS.\n")
      seqlevelsStyle(GR.polyA.db) <- "UCSC"
    }else{
      warning("\n Annotation does not match between ranges and polyAdb.\n Annotation set to Ensembl.\n")
      seqlevelsStyle(GR.polyA.db) <- "Ensembl"
    }
  }
  GR = anchor_3p(ranges)
  GR_cleavage = mutate(GR, width = 1)
  #GR = trim(GR)
  #gr = ranges %>% anchor_3p() %>% mutate(width=1) %>% trim()
  OL = suppressWarnings(distanceToNearest(x = GR_cleavage, subject = GR.polyA.db,
                                          ignore.strand=FALSE))
  keep <- subset(as.data.frame(OL), distance<= max.dist)

  #write this to handle if there are 2 sites equidistant
  #if (sum(duplicated(queryHits(OL))) >0 ) {
  #  dups <- OL[duplicated(queryHits(OL))]
  #}

  #re-write this
  peak.df <- data.frame(GR)
  peak.df <- peak.df[,c("seqnames", "start", "end", "width", "strand")]
  peak.df$peak <- rownames(object[[assay]]@counts)

  #add information about what feature in object matches curated polyA peaks
  tmp <- peak.df[keep$queryHits,]
  tmp2 <- data.frame(GR.polyA.db)[keep$subjectHits,]
  if( atlas %in% c("polyAdbv3-hg38") ){
    tmp2$hg38_Position <- ifelse(tmp2$strand=="+", tmp2$end, tmp2$start)  
  }
  if( atlas %in% c("PolyASitev3.0", "PolyASitev2.0") ){
    tmp2$chromStart <- tmp2$start
    tmp2$chromEnd <- tmp2$end  
  }  
  tmp2 <- tmp2[,!(names(tmp2) %in%  c("strand", "seqnames", "start", "end", "width"))]
  tmp3 <- cbind(tmp, tmp2)

  meta.new <- left_join(peak.df, tmp3, by = c("seqnames", "strand", "peak"))
  rownames(meta.new) <- meta.new$peak
  for( i in 1:ncol(meta.new)){
    object[[assay]] <- AddMetaData(object = object[[assay]] , metadata = meta.new[,colnames(meta.new)[i]] , col.name = colnames(meta.new)[i] )
  }
  return(object)
}
