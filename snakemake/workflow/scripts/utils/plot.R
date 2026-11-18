# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

#' Save a plot to PDF and PNG
#'
#' Saves a ggplot2 (or base R) plot object to both a \code{.pdf} and a
#' \code{.png} file at 600 dpi with the specified dimensions.
#'
#' @param plot     A plot object (e.g. ggplot2 or base R).
#' @param filename character. Base filename without extension.
#' @param width    numeric. Output width in inches.
#' @param height   numeric. Output height in inches.
#'
#' @return NULL (files written to disk as a side effect).
#'
#' @examples
#' \dontrun{
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
#' plot_pdf_and_png(p, "my_plot", width = 6, height = 4)
#' }
plot_pdf_and_png <- function(plot, filename, width, height, res = 600, png_only = FALSE) {

    if (!png_only) {
        pdf(file = paste0(filename, ".pdf"), width = width, height = height)
        print(plot)
        dev.off()
    }

    png(
        file  = paste0(filename, ".png"),
        width = width, height = height,
        res   = res, units = "in"
    )
    print(plot)
    dev.off()
}
