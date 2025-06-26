#' Load and Normalize Nanostring gene expression data
#'
#' @description
#' `loadnormalizercc` is used to read the raw nSolver data saved in RCC file and normalize the data using nSolver or RUVg method.
#'
#' @details
#' This function considers the directory of RCC files and meta data. In the meta data file, the first column
#' is assumed to contain the RCC file names and the subsequent columns are the sample characteristics.
#' Any rows with duplicated entries will be dropped and the first one will be retained.
#'
#' @param dir Path to the NanoString files (make sure the folder contains only the files with .RCC extension).
#' @param meta_dir Path to the meta data file (Consider meta_GSE146530.csv as a template).
#' @param file_number TRUE or FALSE. When TRUE, numbers will be printed when the RCC files are accessed.
#' @param normalization Anyone of the following should be provided: "none" or "nSolver" or "RUVg". Default is "none".
#' @param bd_low The lower limit for binding density (QC Flag). Default is 0.05.
#' @param bd_high The upper limit for binding density (QC Flag). Default is 2.25.
#' @param numSD The number of standard deviations to calibrate the LOD (QC Flag). Default is 0.
#' @param hkgenes A vector of HouseKeeping gene names.
#' @param exclude A vector of HouseKeeping gene to exclude (applicable only to RUVg normalization).
#' @param k The number of dimensions of unwanted variation estimated (applicable only to RUVg normalization).
#' @param rlf Path to the respective RLF file (applicable only to nsolver normalization).
#' @param phenoDataRccColName Column name in meta data that specifies the RCC files (applicable only to nSolver normalization).
#' @return A list of variables containing (not limited to) either a matrix of raw expression values of genes or a matrix
#'   of log2 transformed normalized genes, phenotype data (pdata), feature data (fdata), expression sets and number of
#'   endogenous genes in the samples that are below the limit of detection.
#' @importFrom stringr str_sort str_split
#' @importFrom NanoStringNCTools readNanoStringRccSet setQCFlags normalize
#' @importFrom BiocGenerics counts counts<-
#' @importFrom EDASeq newSeqExpressionSet betweenLaneNormalization
#' @importMethodsFrom EDASeq counts
#' @importFrom RUVSeq RUVg
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersionsGeneEst estimateDispersionsFit varianceStabilizingTransformation
#' @importMethodsFrom DESeq2 counts
#' @importFrom Biobase assayData assayData<- phenoData phenoData<- fData fData<- pData pData<-
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
#' @importFrom rlang .data
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
#' # No normalization -----------------------------------------------------------------
#'
#' exp_df1 <- loadnormalizercc(dir = path_to_rcc, meta_dir = meta,
#'   hkgenes = c("LOC106606270", "LOC106600007"))
#' df1_raw <- exp_df1$rawdata
#' df1_pdata <- exp_df1$pdata
#' df1_fdata <- exp_df1$fdata
#' df1_endblod <- exp_df1$num_endogenous_blod
#'
#' # nSolver normalization ------------------------------------------------------------
#'
#' exp_df2 <- loadnormalizercc(dir = path_to_rcc, meta_dir = meta,
#'   normalization = "nSolver",
#'   hkgenes = c("LOC106606270", "LOC106600007"), phenoDataRccColName = 'files.RCC')
#' df2_nsolverexpset <- exp_df2$nsolverexpset
#' df2_nsolvernormlog2 <- exp_df2$nsolvernormlog2
#' df2_pdata <- exp_df2$pdata
#' df2_fdata <- exp_df2$fdata
#' df2_endblod <- exp_df2$num_endogenous_blod
#'
#' # RUVg normalization ---------------------------------------------------------------
#'
#' exp_df3 <- loadnormalizercc(dir = path_to_rcc, meta_dir = meta,
#'   normalization = "RUVg",
#'   hkgenes = c("LOC106606270", "LOC106600007"), k=3)
#' df3_set <- exp_df3$set
#' df3_vsd <- exp_df3$vsd
#' df3_ruvgnormlog2 <- exp_df3$ruvgnormlog2
#' df3_pdata <- exp_df3$pdata
#' df3_fdata <- exp_df3$fdata
#' df3_endblod <- exp_df3$num_endogenous_blod
#'
#' @author Jayanth Sreekanth

loadnormalizercc  <- function(dir, meta_dir,
                              file_number=TRUE,
                              normalization="none",
                              bd_low=0.05,
                              bd_high=2.25,
                              numSD=0,
                              hkgenes="",
                              exclude = NULL,
                              k,
                              rlf = NULL,
                              phenoDataRccColName="")
{
  if (is.null(dir) | dir == "") {
    stop("The path is empty. Please provide a valid path.")
  } else if (missing(meta_dir)) {
    stop("Please provide the meta data path.")
  } else if (!normalization %in% c("none", "nSolver", "RUVg")) {
    stop("Either define 'none' or 'nSolver' or 'RUVg' to normalization.")
  } else if (all(nchar(hkgenes)) == 0) {
    stop("Housekeeping genes are not provided. Please provide them.")
  } else if (normalization=="RUVg" & missing(k)) {
    stop("'k' should be defined for RUVg normalization")
  } else if (normalization=="nSolver" & (nchar(phenoDataRccColName)== 0)) {
    stop("Please specify the column name in meta data which contains the RCC filenames")
  }

  else {
    # Some scripts adapted from https://github.com/bhattacharya-a-bt/NanoNormIter/tree/master/R
    # Citing https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138885/

    # Determining the RCC files
    files.RCC = list.files(path = dir, full.names = TRUE)
    files.RCC = files.RCC[grepl('RCC',files.RCC)]
    # Retrieving the information of RCC files
    ng = nrow(readRcc(files.RCC[1])$Code_Summary)
    ncol = length(files.RCC)

    raw_expression = as.data.frame(matrix(nrow = ng, ncol = length(files.RCC)+2))
    colnames(raw_expression)[1:2] = c('Gene','Class')
    pData = as.data.frame(matrix(nrow = length(files.RCC),ncol = 11))
    colnames(pData) = c('BCAC_ID','RCC_file_name','Owner','Comments','Date','GeneRLF','SystemAPF','imagingQC',
                        'bindingDensityQC','limitOfDetectionQC','positiveLinearityQC')
    raw_expression[,1:2] = readRcc(files.RCC[1])$Code_Summary[,c(2,1)]

    meta = read.csv(meta_dir)

    # To also include columns from meta data to pdata
    for (ad in 2:ncol(meta)) {
      met_col <- colnames(meta)[ad]
      pData[, met_col] <- NA
    }
    # Reading the RCC files
    for (i in 1:length(files.RCC)){

      if (file_number) {
        print(i)
      }

      rcc = readRcc(files.RCC[i])
      raw = rcc$Code_Summary

      raw_expression[,i+2] = as.numeric(raw$Count)

      colnames(raw_expression)[i+2] = paste0(meta[which(meta[,1] == basename(rcc$File_path)), 2],
                                             '_',
                                             meta[which(meta[,1] == basename(rcc$File_path)), 3])

      pData[i, 2] = basename(rcc$File_path)
      pData[i,3:7] = as.vector(rcc$Sample_Attributes)[2:6]
      # Quality flags
      pData$imagingQC[i] = imagingQC(rcc)
      pData$bindingDensityQC[i] = bindingDensityQC(rcc, bd_low, bd_high)
      pData$limitOfDetectionQC[i] = limitOfDetectionQC(rcc, numSD)
      pData$positiveLinearityQC[i] = positiveLinQC(rcc)
      for (adp in 2:ncol(meta)) {
        pData[i, 11+adp-1] <- meta[which(meta[,1] == basename(rcc$File_path)), adp]
      }
    }
    # Ordering the rows
    pData <- pData %>% arrange(.data$id,.data$time)
    raw_expression <- raw_expression[,c('Gene','Class', str_sort(colnames(raw_expression)[-2:0], numeric = TRUE))]
    # Working on the meta data
    raw = raw_expression[,-c(1:2)]
    fData = raw_expression[,c(1:2)]

    rownames(fData) = fData$Gene
    rownames(raw) = fData$Gene
    rownames(pData) = colnames(raw)
    # Checking the quality of measurements
    neg_raw = raw[grep('NEG_', rownames(raw)), ]
    lod = colMeans(neg_raw) - apply(neg_raw, 2, sd)
    num_endogenous_blod = colSums(raw[fData$Class == 'Endogenous',] < lod)

    rawdata <- as.matrix(raw)

    if (normalization == "none") {
      return(list(rawdata = rawdata, pData = pData, fData = fData,
                  num_endogenous_blod = num_endogenous_blod))
    } else if (normalization == 'nSolver') {
      # Reading the data and saving as nsolver expression set
      nsolver_raw_exp <- readNanoStringRccSet(files.RCC,
                                              rlfFile = rlf,
                                              phenoDataFile = meta_dir,
                                              phenoDataRccColName = phenoDataRccColName)

      # Finding the position of housekeeping gene in codeclass
      pos <- sapply(hkgenes, function(fchar) grep(fchar, rownames(nsolver_raw_exp)))

      # Retrieving the string of characters before the first '_'
      rnam <- str_split(rownames(nsolver_raw_exp), pattern = '\\_', simplify = TRUE)[,1]

      rnam[pos] <- 'Housekeeping'
      nsolver_raw_exp@featureData@data[['CodeClass']] <- rnam

      # Setting Sequencing QC Flags to the dataset as part of nSolver procedure. The default values are given.
      nsolver_raw_exp <- setQCFlags(nsolver_raw_exp,
                                    qcCutoffs = list(Housekeeper = c(failingCutoff = 32, passingCutoff = 100),
                                                     Imaging = c(fovCutoff = 0.75),
                                                     BindingDensity = c(minimumBD = 0.1, maximumBD = 2.25, maximumBDSprint = 1.8),
                                                     ERCCLinearity = c(correlationValue = 0.95),
                                                     ERCCLoD = c(standardDeviations = 2)))

      # Performing nSolver normalization
      nsolver_norm <- normalize(nsolver_raw_exp,
                                type='nSolver',
                                fromELT = "exprs" ,
                                toELT = "exprs_norm")

      # Assigning the normalized count data to a variable
      norm_nsolver_df <- transform(assayData(nsolver_norm)[["exprs_norm"]])
      # Defining proper column names
      colnames(norm_nsolver_df) <- paste0(phenoData(nsolver_norm)[['id']],'_',phenoData(nsolver_norm)[['time']])
      # Sorting the columns for easy visualization
      norm_nsolver_df <- norm_nsolver_df[,c(str_sort(colnames(norm_nsolver_df), numeric = TRUE))]

      # Small function to retrieve the gene names which could be defined as row_names
      extract_rn <- function(rw_name) {
        # if endogenous or housekeeping genes, extract the string between the first and second underscores
        if (grepl("^Endogenous|Housekeeping", rw_name)) {
          return(sub("^[^_]+_([^_]+)_.*", "\\1", rw_name))
          # if positive or negative, extract the string between the first and third underscores
        } else if (grepl("^(Positive|Negative)", rw_name)) {
          w_para <- sub("^[^_]+_([^_]+_[^_]+)_.*", "\\1", rw_name)
          # Removes the parentheses and its contents
          wo_para <- gsub("\\(.*?\\)", "", w_para)
          return(wo_para)
        } else {
          return(NA)
        }
      }

      row.names(norm_nsolver_df) <- sapply(row.names(norm_nsolver_df), extract_rn)

      return(list(nsolverexpset = nsolver_norm,
                  nsolvernormlog2 = as.matrix(norm_nsolver_df),
                  pData = pData, fData = fData,
                  num_endogenous_blod = num_endogenous_blod))

    } else if (normalization == 'RUVg') {

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
                  ruvgnormlog2 = ruvgnormlog2,
                  pData = pData, fData = fData,
                  num_endogenous_blod = num_endogenous_blod))
    }
  }
}

##' @title          Read an .RCC file
##' @description    Parse an .RCC file into a list with each part of the file (Header,
##'                 Sample_Attributes, Lane_Attributes, Code_Summary, etc) stored as a vector
##'                 or data frame.
##'
##' @param  rcc                 Path to the .RCC file.
##' @param  removeSpikeInLabels Logical. If TRUE (the default), RNA \dQuote{spike-in} input labels (if any)
##'                             in the GeneName for positive and negative control probes will be removed.
##'
##' @return
##' A list where each element holds the contents of one part of the .RCC file
##' (Header, Sample_Attributes, Lane_Attributes, Code_Summary, etc) as a
##' vector or data frame.
##'
##' @author NanoStringQCPro

readRcc <- function(rcc, removeSpikeInLabels=TRUE)
{
  # Source: https://github.com/warnes/NanoStringQCPro/blob/master/R/NanoStringConvenience.R
  # Copied from the 'NanoStringQCPro' R package
  # Original license: Artistic-2.0
  # Function modified: No

  suppressWarnings(
    lines <- readLines(rcc)
  )
  lines.df <-
    data.frame(stringsAsFactors=FALSE,
               line   = lines
               ,tag    = rep("", length(lines))
    )

  tags <- c("Header", "Sample_Attributes", "Lane_Attributes", "Code_Summary", "Messages")
  for(i in 1:length(tags))
  {
    tag_start <- which(lines.df$line == paste0("<", tags[i], ">"))
    tag_end <- which(lines.df$line == paste0("</", tags[i], ">"))
    lines.df$tag[ tag_start:tag_end ] <- tags[i]
  }

  unrecognized_lines <- lines.df$line[ lines.df$tag == "" ]
  unrecognized_tags <- grep(value=TRUE, "^<", unrecognized_lines)
  unrecognized_tags <- grep(invert=TRUE, value=TRUE, "^</", unrecognized_tags)
  if (length(unrecognized_tags) > 0) {
    warning(paste0("Unrecognized tags: ", paste(collapse=" ", unrecognized_tags)))
  }

  lines.df <- lines.df[ (lines.df$tag != ""), ]
  lines.df <- lines.df[ !grepl("<", lines.df$line), ]
  lines.df <- lines.df[ !grepl("</", lines.df$line), ]

  ID_lnum <- grep("^ID,", lines.df$line)
  lines.df$line[ID_lnum] <- paste0( sub("_Attributes", "", lines.df$tag[ID_lnum]), lines.df$line[ID_lnum] )

  rcc.ls <- list()

  rcc.ls$File_path <- rcc

  for (attr_tag in c("Header", "Sample_Attributes", "Lane_Attributes")) {
    attr_lnum <- which(lines.df$tag == attr_tag)
    attr.strsplit <- strsplit(lines.df$line[attr_lnum], ",")
    attr.v <- vapply(attr.strsplit, function(x) {if (!is.na(x[2])) {x[2]} else {""}}, character(1))
    names(attr.v) <- vapply(attr.strsplit, function(x) {x[1]}, character(1))
    rcc.ls[[attr_tag]] <- attr.v
  }

  Code_Summary_lnum.all <- which(lines.df$tag == "Code_Summary")
  Code_Summary_lnum.header <- Code_Summary_lnum.all[1]
  Code_Summary_lnum.body <- Code_Summary_lnum.all[2:length(Code_Summary_lnum.all)]

  if (lines.df$line[Code_Summary_lnum.header] != "CodeClass,Name,Accession,Count") {
    stop(
      paste0("Unrecognized header line for Code_Summary section of \"", rcc.ls, "\"",
             " (expected \"CodeClass,Name,Accession,Count\" but got \"", lines.df$line[Code_Summary_lnum.header], "\")")
    )
  }

  Code_Summary_fields <- strsplit(lines.df$line[Code_Summary_lnum.body], ",")
  rcc.ls$Code_Summary <- data.frame(stringsAsFactors=FALSE
                                    ,CodeClass   = vapply(Code_Summary_fields, function(x) {x[1]}, character(1))
                                    ,Name        = vapply(Code_Summary_fields, function(x) {x[2]}, character(1))
                                    ,Accession   = vapply(Code_Summary_fields, function(x) {x[3]}, character(1))
                                    ,Count       = vapply(Code_Summary_fields, function(x) {x[4]}, character(1))
  )

  if(removeSpikeInLabels)
  {
    spikein <- getSpikeInInput(CodeClass=rcc.ls$Code_Summary$CodeClass, GeneName=rcc.ls$Code_Summary$Name)
    rcc.ls$Code_Summary$Name <- spikein$GeneName
  }

  return(rcc.ls)
}

##' @title
##' getSpikeInInput
##'
##' @description
##' Gets the RNA \dQuote{spike-in} input levels for positive and negative control probes from the label in their GeneName.
##' Note that this is a helper function for readRlf() and elsewhere and is not intended for external use.
##'
##' @param CodeClass        Character vector with code classes for each probe.
##' @param GeneName         Character vector with gene names/symbols for each probe.
##' @param nonCtrlProbeVal  Value to assign as the spike-in input for the non-control probes.
##'
##' @return
##' A data frame with the input CodeClass and GeneName but where the latter has been split into
##' two columns: one showing the GeneName for each probe with spike-in input labels removed
##' -- and another with the spike-in input levels.
##'
##' @author NanoStringQCPro

getSpikeInInput <- function(CodeClass, GeneName, nonCtrlProbeVal=NA)
{
  # Source: https://github.com/warnes/NanoStringQCPro/blob/master/R/NanoStringConvenience.R
  # Copied from the 'NanoStringQCPro' R package
  # Original license: Artistic-2.0
  # Function modified: No

  stopifnot(!is.null(CodeClass))
  stopifnot(!is.null(GeneName))
  stopifnot(length(CodeClass) == length(GeneName))

  GeneName.split <- strsplit(GeneName, '\\(')
  GeneName.without_parentheses <- vapply(GeneName.split, function(x) {x[1]}, character(1))
  SpikeInInput.char <- sub('\\)$', '', vapply(GeneName.split, function(x) {x[2]}, character(1)))

  ctrls <- (CodeClass %in% c("Positive", "Negative"))

  GeneName.adjusted <- GeneName
  GeneName.adjusted[ ctrls ] <- GeneName.without_parentheses[ ctrls ]     # Only strip the labels from the control probes (in some cases there will be GeneNames with parentheses!)

  SpikeInInput <- rep(nonCtrlProbeVal, length(GeneName))
  SpikeInInput[ ctrls ] <- as.numeric(SpikeInInput.char[ ctrls ])

  return(
    data.frame(stringsAsFactors=FALSE,
               CodeClass   = CodeClass,
               GeneName    = GeneName.adjusted,
               SpikeInInput = SpikeInInput
    )
  )
}

##' @title          Quality control
##' @description    nSolver imaging quality flag identified in the nSolver manual
##' Scripts copied from https://github.com/bhattacharya-a-bt/NanoNormIter/tree/master/R
##' Citing https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138885/
##'
##' @param  rcc     loaded and valid rcc file
##'
##' @return
##' A character of 'Flag' or 'No flag' for imaging quality

imagingQC <- function(rcc){

  fovRatio = as.numeric(rcc$Lane_Attributes[3]) / as.numeric(rcc$Lane_Attributes[2])
  if (!(fovRatio > .75)) {return('Flag')}
  if (fovRatio > .75) {return('No flag')}

}

##' @title          Quality control - Binding density
##' @description    nSolver binding density quality flag identified in the nSolver manual
##' Scripts copied from https://github.com/bhattacharya-a-bt/NanoNormIter/tree/master/R
##' Citing https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138885/
##'
##' @param  rcc     loaded and valid rcc file
##' @param  low     the lower limit for binding density
##' @param  high    the upper limit for binding density
##'
##' @return
##' A character of 'Flag' or 'No flag' for binding density quality

bindingDensityQC <- function(rcc, low, high){

  bd = as.numeric(rcc$Lane_Attributes[6])
  if(!(bd < high & bd > low)) {return('Flag')}
  if (bd < high & bd > low) {return('No flag')}


}

##' @title          Quality control - Positive line
##' @description    nSolver linearity of positive controls flag identified in the nSolver manual
##' Scripts copied from https://github.com/bhattacharya-a-bt/NanoNormIter/tree/master/R
##' Citing https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138885/
##'
##' @param  rcc     loaded and valid rcc file
##'
##' @return
##' A character of 'Flag' or 'No flag' for linearity of positive controls

positiveLinQC <- function(rcc){

  cnts = rcc$Code_Summary
  posControls = as.numeric(cnts$Count[grepl('POS_',cnts$Name)])
  known = c(128,128/4,128/16,128/64,128/256,128/(256*4))
  r2 = summary(lm(sort(posControls)~sort(known)))$r.squared
  if(!(r2 > .95) | is.na(r2)) {return('Flag')}
  if(r2 > .95) {return('No flag')}

}

##' @title          Quality control - limit of detection
##' @description    nSolver limit of detection flag identified in the nSolver manual
##' Scripts copied from https://github.com/bhattacharya-a-bt/NanoNormIter/tree/master/R
##' Citing https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8138885/
##'
##' @param  rcc     loaded and valid rcc file
##' @param  numSD   number of standard deviations to calibrate the LOD
##'
##' @return
##' A character of 'Flag' or 'No flag' for limit of detection

limitOfDetectionQC <- function(rcc, numSD){
  cnts = rcc$Code_Summary
  posE = as.numeric(cnts$Count[cnts$Name == 'POS_E'])
  negControls = as.numeric(cnts$Count[grepl('NEG',cnts$Name)])
  if(!(posE > mean(negControls) + numSD*sd(negControls))) {return('Flag')}
  if (posE > mean(negControls) + numSD*sd(negControls)) {return('No flag')}
}
