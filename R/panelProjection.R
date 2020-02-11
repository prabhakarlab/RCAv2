#' Compute Reference Component features for clustering analysis
#'
#' @param rca.obj RCA object.
#' @param panel Reference panel as a matrix
#' @param corMeth Any of the correlation measures supported by R, defaults to pearson
#' @param power power to raise up to for the RCA features before clustering, default is 4
#' @param scale True if the data should be scaled, False otherwise
#' @param apply_threshold True if a threshold to raise the lowest expression values needs to be applied on panel data, False otherwise.
#' @param threshold Minimal wished expression value for panel data. All lower expression values in panel will be raised to this value if apply_threshold is TRUE. 
#' @return RCA object.
#' @export
#'

panelProjection = function(sc_data, panel, corMeth="pearson", power=4, scale=TRUE, apply_threshold=FALSE, threshold=NULL) {
    
    # Select genes that are shared by the input data and the panel
    shared_genes <- intersect(rownames(sc_data), rownames(panel))
    print(paste0("Projection on ", length(shared_genes), " genes."))
    
    # Reduce the panel and input data to the shared genes
    subset_panel = panel[shared_genes, ]
    subset_data = sc_data[shared_genes, , drop = FALSE]
    
    # For values in the panel below the minimum threshold, set those values to threshold
    if (apply_threshold) {
        subset_panel[subset_panel <= threshold] = threshold
    }
    
    # Compute projection of input data with the panel
    if(corMeth == "pearson") {
        subset_panel = as.matrix(subset_panel)
        projection <- qlcMatrix::corSparse(X = subset_panel, Y = subset_data)
    } else {
        projection <- cor(subset_panel, subset_data, method = corMeth)
    }
    rownames(projection) <- colnames(subset_panel)
    colnames(projection) <- colnames(subset_data)
    
    # Raise the projection to power
    projection = abs(projection) ^ (power) * sign(projection)
    
    # If scaling is required
    if (scale) {
        # Scale
        projection = scale(projection,
                           center = TRUE,
                           scale = TRUE)
    }
    
    return(projection)
}