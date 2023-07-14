#' Plot Heatmap
#'
#' This function generates a heatmap plot based on mutation data and gene expression data from cell lines.
#'
#' @param se A summarizedExperiment object containing mutation data.
#' @param mutations A character vector specifying the mutations of interest.
#' @param genes A character vector specifying the genes of interest. Default is NULL.
#' @param group_var A variable (column) name included in colData(se). Can be used to subset the dataset.
#' @param group_lvl A level present in group_var. Only cell lines that are included in this subgroup are taken into account. Required if group_var is provided.
#' @param cell_lines A character vector specifying the cell lines to include in the plot. Default is NULL, which means that all cell lines will be included (except if group_Var/group_lvl is specified)
#' @param remove_na A logical value indicating whether to remove cell lines that have no data on mutations. Default is set to TRUE
#'
#' @return Plots oncoplots showing which mutations are present in cell lines and optionally also gene expression data. In addition,  list containing the mutations and expression data if genes is provided. Otherwise, returns the mutations.
#'
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @importFrom grid gpar grid.rect grid.points unit
#'
#' @export
#'
#' @examples
#' plotHeatmap(se, mutations = c("APC", "KRAS"), genes = c("TP53", "CDK9"),
#' group_var = "oncoTreeCode", group_lvl = "COAD")

plotHeatmap = function(se, mutations, genes = NULL, group_var = NULL, group_lvl = NULL, cell_lines = NULL, remove_na = TRUE){

  if (!is.null(cell_lines)){

    cell_lines = colnames(se)[colnames(se) %in% cell_lines]
    se = se[,cell_lines]
  }

  else if (!is.null(group_var)){

    if (!(group_var %in% colnames(colData(se)))){
      message('group_var not present in colData.')
      col_vars = paste(colnames(colData(se)), collapse = '\n')
      message(paste('The following options are available:\n', col_vars))
      return()
    }

    if (is.null(group_lvl)) {
      message('Please specify group_lvl.')
      return()
    }

    else if (!(group_lvl %in% colData(se)[[group_var]])){
      message(paste('group_lvl not present in', group_var))
      group_lvls = paste(unique(colData(se)[[group_var]]), collapse = '\n')
      message(paste('The following options are available:\n', group_lvls))
      return()
    }

    se = se[,se[[group_var]] %in% group_lvl]
  }

  mutations = SummarizedExperiment::assay(se, 'mutations')[mutations,]

  if (remove_na){
    mutations = mutations[,colSums(is.na(mutations)) != nrow(mutations)]
  }

  mutations[is.na(mutations)|mutations == 'wt'] = ' '
  mutations[mutations == 'both'] = 'mut;both'
  mutations[mutations == 'cosmic'] = 'mut;cosmic'
  mutations[mutations == 'depmap'] = 'mut;depmap'

  alter_fun = list(
    # red rectangles
    mut = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red4", col = NA)),
    # dots
    cosmic = function(x, y, w, h)
      grid.points(x, y, pch = 15, size = unit(0.8, 'char'), gp = gpar(col = 'dodgerblue4')),
    # crossed lines
    depmap = function(x, y, w, h)
      grid.points(x, y, pch = 15, size = unit(0.8, 'char'), gp = gpar(col = 'seagreen')),
    both = function(x, y, w, h)
      grid.points(x, y, pch = 15, size = unit(0.8, 'char'), gp = gpar(col = 'wheat3'))

  )

  col = c('mut' = 'darkred')

  op = ComplexHeatmap::oncoPrint(mutations, alter_fun = alter_fun, col = col)

  if (!is.null(genes)){

    expr = SummarizedExperiment::assay(se, 'reads_vst')
    expr = t(scale(t(expr)))
    expr = expr[genes, colnames(expr) %in% colnames(mutations)]

    hm = ComplexHeatmap::Heatmap(expr, name = 'Z-scores', show_row_dend = F)

    hm_list = op %v% hm
    draw(hm_list)
    return(list(mutations = mutations, expression = expr))
  }
  else {

    draw(op)
    return(mutations)
  }



}
