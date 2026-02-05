#' Normalize Nanostring gene expression data
#'
#' @description
#' `normalizeruvg` is used to normalize the data using RUVg method.
#'
#' @details
#' This function is useful when nanostring data has been loaded, quality-checked and is ready to be normalized.
#' Only RUVg method of normalization can be applied.
#'
#' @param rawdata A matrix of raw expression values of genes.
#' @param fData Feature data.
#' @param pData Phenotype data.
#' @param hkgenes A vector of HouseKeeping gene names.
#' @param exclude A vector of HouseKeeping gene to exclude.
#' @param k The number of dimensions of unwanted variation estimated.
#' @return A matrix of log2 transformed normalized genes.
#' @importFrom BiocGenerics counts counts<-
#' @importFrom EDASeq newSeqExpressionSet betweenLaneNormalization
#' @importMethodsFrom EDASeq counts
#' @importFrom RUVSeq RUVg
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersionsGeneEst estimateDispersionsFit varianceStabilizingTransformation
#' @importMethodsFrom DESeq2 counts
#' @importFrom Biobase assayData assayData<- phenoData phenoData<- fData fData<- pData pData<-
#' @importFrom magrittr %>%
#' @importFrom MatrixGenerics rowVars rowMeans
#' @importFrom SummarizedExperiment assay assay<- rowData rowData<- colData
#' @importFrom limma removeBatchEffect
#' @importFrom stats lm sd
#' @importFrom utils read.csv
#' @importFrom S4Vectors mcols mcols<-
#' @export
#' @examples
#'
#' path_to_rcc <- system.file("extdata", "GSE146530", package = "rhythmCI")
#' meta <- system.file("extdata", "meta_GSE146530.csv", package = "rhythmCI")
#'
#' # Load and No normalization -----------------------------------------------------------------
#'
#' exp_df1 <- loadnormalizercc(dir = path_to_rcc, meta_dir = meta,
#'   hkgenes = c("LOC106606270", "LOC106600007"), file_number = FALSE)
#' df1_raw <- exp_df1$rawdata
#' df1_pdata <- exp_df1$pData
#' df1_fdata <- exp_df1$fData
#' df1_endblod <- exp_df1$num_endogenous_blod
#'
#' # Perform nSolver normalization ------------------------------------------------------------
#' suppressMessages(
#'
#' exp_df3 <- normalizeruvg(rawdata = df1_raw,
#'   fData = df1_fdata,
#'   pData = df1_pdata,
#'   hkgenes = c("LOC106606270", "LOC106600007"),
#'   exclude = NULL,
#'   k=3)
#' )
#'
#' df3_set <- exp_df3$set
#' df3_vsd <- exp_df3$vsd
#' df3_ruvgnormlog2 <- exp_df3$ruvgnormlog2
#'
#'
#' @author Jayanth Sreekanth

normalizeruvg  <- function(rawdata,
                              fData,
                              pData,
                              hkgenes="",
                              exclude = NULL,
                              k)
{
  if (all(nchar(hkgenes)) == 0) {
    stop("Housekeeping genes are not provided. Please provide them.")
  } else if (missing(rawdata)) {
    stop("'rawdata' for normalization is missing. Please provide it.")
  } else if (missing(fData)) {
    stop("'fData' is missing. Please provide it.")
  } else if (missing(pData)) {
    stop("'pData' is missing. Please provide it.")
  } else if (missing(k)) {
    stop("'k' should be defined for RUVg normalization")
  }

  else {

      fData = fData[rownames(rawdata), ]
      int = intersect(rownames(rawdata), rownames(fData))
      fData = fData[int, ]
      rawdata = rawdata[int, ]

      set <- newSeqExpressionSet(counts = as.matrix(round(rawdata)),
                                 phenoData = pData,
                                 featureData = fData)

      if (!is.null(hkgenes)) {
        fData(set)$Class[rownames(set) %in% hkgenes] = "Housekeeping"
      }

      cIdx <- rownames(set)[fData(set)$Class == "Housekeeping"]
      if (!is.null(exclude)) {
        cat("Excluding:", exclude, "\n")
      }
      cIdx = cIdx[!(cIdx %in% exclude)] # any housekeeping genes to exclude

      set <- betweenLaneNormalization(set, which = "upper")
      counts <- BiocGenerics::counts
      set <- RUVSeq::RUVg(set, cIdx, k = k)

      dds <- DESeqDataSetFromMatrix(EDASeq::counts(set), colData = pData(set), #counts(set)
                                    design = ~1)

      rowData(dds) <- fData
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersionsGeneEst(dds)

      cts <- DESeq2::counts(dds, normalized = TRUE)

      disp <- pmax((rowVars(cts) - rowMeans(cts)), 0)/rowMeans(cts)^2
      mcols(dds)$dispGeneEst <- disp
      dds <- estimateDispersionsFit(dds, fitType = "mean")
      vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
      mat <- assay(vsd)
      covars <- as.matrix(colData(dds)[, grep("W", colnames(colData(dds))),
                                       drop = FALSE])
      mat <- removeBatchEffect(mat, covariates = covars)
      assay(vsd) <- mat
      ruvgnormlog2 <- log2(set@assayData[["normalizedCounts"]])

      return(list(set = set, vsd = vsd,
                  ruvgnormlog2 = ruvgnormlog2))
  }
}

