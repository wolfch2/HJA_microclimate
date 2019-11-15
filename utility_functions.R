# Buffer large clumps of NA cells in raster to remove artifacts
# around edges.  Clumps only buffered if their size exceeds threshold.
# "drop" is the number of times to repeat buffering
buffer_inward = function(raster,drop,threshold=1e4){
        require(raster)

        rast_NA = rast
        rast_NA[] = as.numeric(is.na(rast[]))
        rast_clump = clump(rast_NA, directions=8) # NA and 0 are bg values for clumping

        outside = rast * 0
        outside[rast_clump[] %in% names(which(table(rast_clump[]) > 1e4))] = NA

        for(x in 1:drop){ # drop outer cells this many times
                print(x)
                edge = boundaries(outside, type='inner', classes=FALSE, directions=8)
                outside[edge[] %in% 1] = NA
        } # expands outside inward by Drop_cells # of cells

        rast[is.na(outside)] = NA

        # we can also have errors along borders at the extents of the raster -- just drop these parts directly:
        for(drop_row in c(1:drop, nrow(rast):(nrow(rast)-drop+1))) rast[drop_row,] = NA
        for(drop_col in c(1:drop, ncol(rast):(ncol(rast)-drop+1))) rast[,drop_col] = NA 

        return(rast)
}

# pretty formatting for response and predictor variable names
format_names = function(names){
	plyr::revalue(names, replace = c(
		"Apr_Jun_mean_min"="Spring min",
		"Apr_Jun_mean_max"="Spring max",
		"Apr_Jun_mean_mean"="Spring mean",
		"Jul_Sep_mean_max"="Summer max",
		"Jul_Sep_mean_mean"="Summer mean",
		"GDD_winter_5"=paste0("Winter GDD (5 ","\U00B0", "C)"),
		"Closure_2"="Closure (> 2 m)",
		"Closure_10"="Closure (> 10 m)",
		"Closure_40"="Closure (> 40 m)",
		"Density_0_2"="Density (0-2 m)",
		"Density_2_10"="Density (2-10 m)",
		"Height_mean"="Height (mean)",
		"Height_75"="Height (75%)",
		"Height_95"="Height (95%)",
		"Height_SD"="Height (s.d.)",
		"Height_SD_7x7"="Height (7x7 s.d.)",
		"Height_80_1m"="Height (80%, 1 m)",
		"Height_95_1m"="Height (95%, 1 m)",
		"Height_SD_1m"="Height (s.d., 1 m)",
		"Biomass_Mg_Ha"="Biomass (Mg/Ha)",
		"Veg_height_DEM"="Veg. height (DEM)",
		"Elevation"="Elevation",
		"Slope"="Slope",
		"Aspect_eastness"="Aspect (eastness)",
		"Aspect_northness"="Aspect (northness)",
		"Convergence_index"="Convergence index",
		"Position_index"="Position index",
                "PC1"="PC 1",
                "PC2"="PC 2"
	), warn_missing = TRUE)
}

