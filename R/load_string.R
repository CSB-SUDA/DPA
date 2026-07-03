
#============================================================
# Load STRING PPI database
#
# Args:
#   ppi_dir       : directory containing STRING database
#   score_cutoff  : STRING confidence score cutoff
#
# Returns:
#   data.frame(protein1, protein2)
#============================================================

load_string <- function(ppi_dir = "data",
                        score_cutoff = 400){
  
  ## read STRING interactions
  ppi <- read.table(
    file.path(ppi_dir, "9606.protein.links.v11.5.txt.gz"),
    header = TRUE,
    sep = " "
  )
  
  ppi <- ppi[ppi$combined_score >= score_cutoff, ]
  
  ## read protein annotation
  map <- read.delim(
    file.path(ppi_dir, "9606.protein.info.v11.5.txt"),
    header = TRUE,
    quote = ""
  )
  
  map <- map[, 1:2]
  
  colnames(map) <- c(
    "ensembl_peptide_id",
    "geneSYMBOL"
  )
  
  ## convert ENSP → Gene Symbol
  ppi$protein1 <- map[
    match(ppi$protein1, map$ensembl_peptide_id),
    "geneSYMBOL"
  ]
  
  ppi$protein2 <- map[
    match(ppi$protein2, map$ensembl_peptide_id),
    "geneSYMBOL"
  ]
  
  ## remove NA
  ppi <- ppi[
    !is.na(ppi$protein1) &
      !is.na(ppi$protein2),
  ]
  
  ## remove duplicated edges
  g <- graph_from_data_frame(
    ppi[, c("protein1","protein2")],
    directed = FALSE
  )
  
  g <- simplify(g)
  
  ppi <- data.frame(
    as_edgelist(g),
    stringsAsFactors = FALSE
  )
  
  colnames(ppi) <- c(
    "protein1",
    "protein2"
  )
  
  return(ppi)
  
}