#' Function to combine base reference panels for one projection.
#' 
#' @param rca.obj RCA object. 
#' @param panel.name Vector of methods from dataProject to select panels ("GlobalPanel", "ColonEpitheliumPanel", "ENCODEPanel" and "Custom").
#' @param customPath Character vector with corresponding customPath to "Custom" panel name for dataProject.
#' @param corMeth Character vector with corresponding customPath to "Custom" panel name for dataProject.
#' @param power Numeric vector for corresponding power parameter to pass to dataProject for each method in panel.name.
#' @param scale Boolean vector for corresponding scale parameter to pass to dataProject for each method in panel.name.
#' @param ctSelection List of character vector for corresponding ctSelect parameter to pass to dataProject for each method in panel.name.
#'
#' @export
#' 
multiRCAPanel = function(rca.obj, 
                         panel.name=c("GlobalPanel", "ENCODEPanel", "ColonEpitheliumPanel"),
                         customPath = rep(NULL, length(panel.name)),
                         corMeth = rep("pearson", length(panel.name)),
                         power = rep(4, length(panel.name)),
                         scale = rep(TRUE, length(panel.name)),
                         ctSelection = rep(NULL, length(panel.name))) {
    
    full.projection = list()
    for (i in 1:length(panel.name)) {
        panel = panel.name[i]
        print(panel)
        if (panel %in% c("GlobalPanel", "ENCODEPanel", "ColonEpitheliumPanel", "Custom")) {
            dataProject(rca.obj = rca.obj, method = panel,
                               corMeth = corMeth[i], power = power[i],
                               scale = scale[i], ctSelection = ctSelection[[i]])
            full.projection[[panel]] = rca.obj$projection.data
        } else {
            print("panel.name must be methods from dataProject().")
        }
    }
    
    rca.obj$projection.data = do.call("rbind", full.projection)
    return(rca.obj)
}

