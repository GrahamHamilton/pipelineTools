#' Save pheatmap
#'
#' @description Save the pheatmap object as a png file
#'
#' @import grDevices
#'
#' @param plot pheatmap oject
#' @param filename file name or file path name to write the pheatmap object
#' @param path Path of the directory to save plot
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

save_pheatmap_png <- function(plot = NULL,
                              filename = NULL,
                              path = NULL,
                              width = 1200,
                              height = 1000,
                              res = 150) {
  # Path to save file
  if (!is.null(path)){
    file_path <- file.path(path, filename, fsep = .Platform$file.sep)
  }else{
    file_path <- filename
  }

  png(file_path, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(plot$gtable)
  dev.off()
}
