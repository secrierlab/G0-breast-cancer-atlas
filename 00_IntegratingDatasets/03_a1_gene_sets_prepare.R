gene_sets_prepare <- function(path_to_db_file, cell_type) {
  cell_markers = path_to_db_file
  cell_markers = cell_markers[cell_markers$tissueType == cell_type, ]
  cell_markers$geneSymbolmore1 = gsub(" ", "", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub(" ", "", cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i) {
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i], ",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if (length(markers_all) > 0) {
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse = ",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i) {
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i], ",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if (length(markers_all) > 0) {
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse = ",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///", ",", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore1 = gsub(" ", "", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///", ",", cell_markers$geneSymbolmore2)
  cell_markers$geneSymbolmore2 = gsub(" ", "", cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]), "," ))))
  names(gs) = cell_markers$cellName
  
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]), ","))))
  names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}
