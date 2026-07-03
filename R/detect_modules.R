#============================================================
# Detect network modules and identify module evolution
#
# Args:
#   network      : state-specific weighted PPI network
#   resolution   : Louvain resolution
#   baseline     : baseline state
#   min_module   : minimum module size
#   method       : "Ochiai" or "Jaccard"
#   seed         : random seed
#
# Returns:
#   list(
#      modules,
#      gene_module,
#      module_change
#   )
#============================================================


detect_modules <- function(network,
                           resolution = 2,
                           seed = 1,
                           baseline = 1,
                           min_module = 20,
                           method = "Ochiai"){
  
  ## edge list
  edge <- network[, 1:2]
  
  ## graph
  g <- graph_from_data_frame(
    d = edge,
    directed = FALSE
  )
  
  ## 所有权重列（以 .w 结尾）
  weight_cols <- grep("\\.w$", colnames(network), value = TRUE)
  
  modules <- vector("list", length(weight_cols))
  names(modules) <- sub("\\.w$", "", weight_cols)
  
  for(i in seq_along(weight_cols)){
    
    cat("Detecting modules:", names(modules)[i], "\n")
    
    set.seed(seed)
    
    clu <- cluster_louvain(
      g,
      weights = network[[weight_cols[i]]],
      resolution = resolution
    )
    
    modules[[i]] <- data.frame(
      gene = clu$names,
      module = clu$membership,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    
  }
  
  resTable <- data.frame(
    gene = modules[[1]]$gene,
    do.call(
      cbind,
      lapply(modules, function(x) x$module)
    ),
    check.names = FALSE
  )
  colnames(resTable) <-
    c("gene", names(modules))
  
  cat("Calculate the evolutionary path and conservative score\n")
  resScore <- .search_best_match(
    modules,
    baselineIndex = baseline,
    minmodule = min_module,
    method = method
  )
  
  return(list(
    modules = modules,
    gene_module = resTable,
    module_change = resScore
  ))
}


.search_best_match <- function(modules,
                               baselineIndex = 1,
                               minmodule = 20,
                               method = "Ochiai") {
  
  ## baseline module
  baseline <- modules[[baselineIndex]]
  
  baselinetable <- data.frame(table(baseline$module))
  baselinetable <- baselinetable[baselinetable$Freq >= minmodule, ]
  baselinetable <- baselinetable[order(baselinetable$Freq, decreasing = TRUE), ]
  
  baseline <- baseline[
    baseline$module %in% as.character(baselinetable$Var1),
    ,
    drop = FALSE
  ]
  
  baselinematrix <- matrix(
    0,
    nrow = nrow(baselinetable),
    ncol = 4 * length(modules)
  )
  
  baselinematrix[, 1] <- as.character(baselinetable$Var1)
  
  seqline <- seq_along(modules)
  seqline <- seqline[-baselineIndex]
  
  for (i in seqline) {
    
    indexi <- modules[[i]]
    
    indexitable <- data.frame(table(indexi$module))
    indexitable <- indexitable[indexitable$Freq >= minmodule, ]
    indexitable <- indexitable[order(indexitable$Freq, decreasing = TRUE), ]
    
    indexi <- indexi[
      indexi$module %in% indexitable$Var1,
      ,
      drop = FALSE
    ]
    
    ##--------------------------------------------------
    ## similarity matrix
    ##--------------------------------------------------
    
    matrixi <- matrix(
      0,
      nrow = nrow(baselinetable),
      ncol = nrow(indexitable)
    )
    
    matrixCount <- matrix(
      0,
      nrow = nrow(baselinetable),
      ncol = nrow(indexitable)
    )
    
    for (modulesi in seq_len(nrow(baselinetable))) {
      
      modulesiCode <- as.character(baselinetable$Var1[modulesi])
      
      modulesiGene <- rownames(
        baseline[
          baseline$module == modulesiCode,
          ,
          drop = FALSE
        ]
      )
      
      for (modulesj in seq_len(nrow(indexitable))) {
        
        modulesjCode <- as.character(indexitable$Var1[modulesj])
        
        modulesjGene <- rownames(
          indexi[
            indexi$module == modulesjCode,
            ,
            drop = FALSE
          ]
        )
        
        inters <- length(intersect(modulesiGene, modulesjGene))
        
        if (method == "Jaccard") {
          
          tempScore <- inters /
            length(union(modulesiGene, modulesjGene))
          
        } else if (method == "Ochiai") {
          
          tempScore <- inters /
            sqrt(length(modulesiGene) * length(modulesjGene))
          
        } else {
          
          stop("method should be 'Jaccard' or 'Ochiai'.")
          
        }
        
        matrixCount[modulesi, modulesj] <- inters
        matrixi[modulesi, modulesj] <- tempScore
      }
    }
    
    rownames(matrixi) <- baselinetable$Var1
    colnames(matrixi) <- indexitable$Var1
    
    rownames(matrixCount) <- baselinetable$Var1
    colnames(matrixCount) <- indexitable$Var1
    
    ##--------------------------------------------------
    ## best match
    ##--------------------------------------------------
    
    nodev <- NULL
    edgev <- NULL
    nodeC <- NULL
    intersC <- NULL
    
    for (modulesi in seq_len(nrow(baselinetable))) {
      
      maxIndex <- which.max(matrixi[modulesi, ])
      
      nodev <- c(nodev, colnames(matrixi)[maxIndex])
      
      edgev <- c(
        edgev,
        round(matrixi[modulesi, maxIndex], 4)
      )
      
      nodeC <- c(
        nodeC,
        indexitable$Freq[
          indexitable$Var1 == colnames(matrixi)[maxIndex]
        ]
      )
      
      intersC <- c(
        intersC,
        matrixCount[modulesi, maxIndex]
      )
    }
    
    baselinematrix[, i] <- nodev
    baselinematrix[, i + length(modules)] <- edgev
    baselinematrix[, i + length(modules) * 2] <- nodeC
    baselinematrix[, 1 + length(modules) * 2] <- baselinetable$Freq
    baselinematrix[, i + length(modules) * 3] <- intersC

  }
  
  baselinematrix <- data.frame(baselinematrix)
  
  colnames(baselinematrix) <- c(
    names(modules),    #each state
    paste0("consScore", seq_along(modules)-1),
    paste0("nodeCount", seq_along(modules)-1),
    paste0("intersectCount", seq_along(modules)-1)
  )
  
  score_cols <- grep("^consScore", colnames(baselinematrix))
  baselinematrix[, score_cols] <-
    lapply(
      baselinematrix[, score_cols, drop = FALSE],
      as.numeric
    )
  
  ## 去除 baseline（第一个）
  score_cols <- score_cols[-baselineIndex]
  baselinematrix$consScoreAve <-
    rowMeans(
      baselinematrix[, score_cols, drop = FALSE],
      na.rm = TRUE
    )
  baselinematrix <- baselinematrix[order(baselinematrix$consScoreAve,
                                         decreasing = TRUE),]
  
  return(baselinematrix)
}
