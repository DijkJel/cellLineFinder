#' Load summarizedExperiment
#'
#' This function combines multiple summarizedExperiment objects into a single object.
#'
#' @details This function is used to merge multiple summarizedExperiment objects into a single object by row-binding them together.
#'
#' @return A summarizedExperiment object containing the merged data from multiple summarizedExperiment objects. It contains cosmic and depmap data on cell lines.
#'
#' @examples
#' # Example usage
#' se <- load_se()
#'
#' @export
load_se = function(){

  se = rbind(se1, se2)
  return(se)
}
