#' Create a RLE plot customized to time-series data
#'
#' @description
#' `create_RLEplot` generates a RLE plot with samples colored according to subjects.
#'
#' @details
#' Only a raw expression matrix of genes or a log2 transformed matrix of normalized genes
#' with properly formatted column names is needed to create a plot.
#'
#' @param exp_matrix Matrix of the raw expression of genes or log2 transformed normalized genes.
#' @param ylim Numeric vector of length 2 for the y-coordinate range.
#' @return A RLE plot
#' @importFrom stringr str_split
#' @importFrom scales hue_pal
#' @importFrom MatrixGenerics rowMedians
#' @importFrom graphics boxplot title text axis abline box par
#' @export
#' @examples
#'
#' # Path of the sample RCC files
#' path_to_rcc <- system.file("extdata", "GSE146530", package = "rhythmCI")
#' meta <- system.file("extdata", "meta_GSE146530.csv", package = "rhythmCI")
#'
#' # Load and No normalization -----------------------------------------------------------------
#' suppressWarnings(
#'
#' exp_df1 <- loadnormalizercc(dir = path_to_rcc, meta_dir = meta,
#'                             hkgenes = c("LOC106606270", "LOC106600007"),
#'                             file_number = FALSE
#'                             )
#' )
#'
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
#' create_RLEplot(df3_ruvgnormlog2, ylim = c(-2,2))
#'
#'@author Jayanth Sreekanth

create_RLEplot <- function(exp_matrix, ylim){

  grp <- str_split(colnames(exp_matrix), pattern = '\\_', simplify = TRUE)

  mat_centered <- sweep(exp_matrix,
                        1,
                        rowMedians(exp_matrix),
                        FUN = '-')

  grp_fac <- factor(grp)

  colors <- hue_pal(l=45, c = 100)(nlevels(grp_fac))

  box_cols <- colors[grp_fac]

  boxplot(mat_centered, axes = F,
          lwd=2, ylim = ylim, col = "white", border = box_cols, outline = F)
  title('Relative Log Expression plot', ylab = substitute(paste(bold('Median deviation of log expression'))),
        line=2.5, cex.main = 2, cex.lab = 1.3)
  axis(2, mgp = c(1, 1, 0), lwd.ticks=2, las=2)
  abline(h = 0, col = 'black', lwd = 2, lty=2)
  box(lwd = 2, bty = 'l')
  par(lwd = 2)

  # Dynamically position labels slightly below the lower y-axis limit
  y_pos <- ylim[1] - (diff(ylim) * 0.05)

  text(x = seq_along(colnames(exp_matrix)), y = y_pos,
       xpd = TRUE,
       labels = colnames(exp_matrix),
       srt = 90,
       adj = 1)

}
