
returnCellLines = function(se, mutations, genes = NULL, group_var = NULL, group_lvl = NULL){


  if (!is.null(group_var)){

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

  se_mutations = se[rownames(se) %in% mutations,]

  mutations_m = assay(se_mutations, 'mutations')
  mutations_TF = mutations_m
  mutations_TF[] = F
  mutations_TF[mutations_m != 'wt'&!is.na(mutations_m)] = T

  mutations_TF = matrix(as.logical(mutations_TF), nrow = length(mutations), dimnames = dimnames(mutations_TF))
  mutations_TF = mutations_TF[,colSums(mutations_TF) == nrow(mutations_TF), drop = F]

  se_expression = se[rownames(se) %in% genes,colnames(se) %in% colnames(mutations_TF)]
  reads = assay(se_expression, 'reads_vst')
  reads = t(scale(t(reads)))

  df = data.frame(cell_line = colnames(mutations_TF))
  df = cbind(df, t(reads))

  return(df)


}
