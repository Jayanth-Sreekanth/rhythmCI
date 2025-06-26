#' Create a RLE plot customized to time-series data
#'
#' @description
#' `create_RLEplot` generates a RLE plot with samples colored according to subjects.
#'
#' @details
#' Only two variables are needed to create a plot, a raw expression matrix of genes or a log2 transformed
#' matrix of normalized genes and the number of unique subjects.
#'
#' @param exp_matrix Matrix of the raw expression of genes or log2 transformed normalized genes.
#' @param num The number of unique subjects (used to generate hue).
#' @return A RLE plot
#' @author Jayanth Sreekanth
#' @importFrom stringr str_split
#' @importFrom scales hue_pal
#' @importFrom MatrixGenerics rowMedians
#' @importFrom graphics boxplot title text axis abline box par
#' @export
#' @examples
#'
#' # Create a sample gene expression matrix
#' genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")
#' samples <- c("Pat1_0", "Pat1_4", "Pat1_8", "Pat2_0", "Pat2_4", "Pat2_8")
#'
#' df3_ruvgnormlog2 <- matrix(c(6.2, 6.5, 6.0, 7.0, 5.8, 6.1,
#'                             7.5, 7.8, 7.2, 8.0, 7.0, 7.6,
#'                             5.9, 6.1, 5.6, 6.3, 5.4, 5.8,
#'                             8.1, 8.3, 8.0, 8.5, 7.9, 8.2,
#'                             4.2, 4.5, 4.0, 4.6, 4.1, 5.0),
#'                             nrow = 5,
#'                             byrow = TRUE)
#'
#' rownames(df3_ruvgnormlog2) <- genes
#' colnames(df3_ruvgnormlog2) <- samples
#'
#' create_RLEplot(df3_ruvgnormlog2, num=2)
#'
create_RLEplot <- function(exp_matrix, num){
  grp <- str_split(colnames(exp_matrix), pattern = '\\_', simplify = TRUE)

  colors <- hue_pal(l=25)(num)
  boxplot(exp_matrix - rowMedians(exp_matrix), axes = F,
          lwd=2, ylim = c(-5,5), boxcol = colors[factor(grp)], border = colors[factor(grp)], outline = F,
          col = 'white')
  title('Relative Log Expression plot', ylab = substitute(paste(bold('Median deviation of log expression'))),
        line=2.5)
  axis(2, mgp = c(1, 1, 0), lwd.ticks=2, las=2)
  abline(h = 0, col = 'black', lwd = 2, lty=2)
  box(lwd = 2, bty = 'l')
  par(lwd = 2)

  text(x = seq_along(colnames(exp_matrix)), y = -5.5,
       xpd = TRUE,
       labels = colnames(exp_matrix),
       srt = 90,
       adj = 1)
}
