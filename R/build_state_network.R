#============================================================
# Build state-specific weighted PPI network
#
# Args:
#   expr        : gene × sample expression matrix
#   phenotype   : phenotype data.frame
#   ppi         : STRING interaction table
#   state_col   : phenotype state column
#   sample_col  : phenotype sample column
#   method      : correlation method
#
# Returns:
#   weighted PPI network
#============================================================

build_state_network <- function(
    expr,
    phenotype,
    ppi,
    state_col = "phase",
    sample_col = "sampleID",
    method = "pearson"
){
  
  ## sample order check
  stopifnot(
    all(phenotype[[sample_col]] ==
          colnames(expr))
  )
  
  ## retain genes in expression matrix
  edge <- ppi[
    ppi$protein1 %in% rownames(expr) &
      ppi$protein2 %in% rownames(expr),
  ]
  
  cat(
    "Network nodes:",
    length(unique(c(edge$protein1, edge$protein2))),
    "\n"
  )
  
  cat(
    "Expression genes:",
    nrow(expr),
    "\n"
  )
  
  ## disease stages
  states <- sort(unique(phenotype[[state_col]]))
  
  i=1
  for(st in states){
    
    cat("Processing state",i,":", st, "\n")
    
    samples <- phenotype[
      phenotype[[state_col]] == st,
      sample_col
    ]
    
    cor_mat <- cor(
      t(expr[, samples]),
      method = method
    )
    
    idx1 <- match(edge$protein1,
                  rownames(cor_mat))
    
    idx2 <- match(edge$protein2,
                  rownames(cor_mat))
    
    cor_value <- cor_mat[
      cbind(idx1, idx2)
    ]
    
    edge[[st]] <- cor_value
    
    edge[[paste0(st,".w")]] <-
      abs(
        atanh(cor_value)
      )
    
    i=i+1
  }
  
  return(edge)
  
}