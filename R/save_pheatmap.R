#' Save pheatmap
#'
#' @description Save the pheatmap object as a png fileßß
#'
#' @param pheatmap pheatmap oject
#' @param filename file name or file path name to write the pheatmap object to
#' @param width Width of plot
#' @param height Height of plot
#' @param res Resolution of plot
#'
#' @examples
#' \dontrun{
#'
#' ph <- pheatmap(log2_mean_differential_genes
#'
#' save_pheatmap_png(ph, "heatmap.png")
#' }
#'
#' @return No return object
#'
#' @export
#'

save_pheatmap_png <- function(pheatmap = NULL,
                              filename = NULL,
                              width = 1200,
                              height = 1000,
                              res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(pheatmap$gtable)
  dev.off()
}
