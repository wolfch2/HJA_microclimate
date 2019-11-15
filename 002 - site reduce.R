# see setup description in Sophie folder
LD_LIBRARY_PATH=/home/stats/wolfch/lib:/cm/local/apps/gcc/5.1.0/lib64:/cm/shared/apps/gdal/2.1.1/lib:/cm/shared/apps/gdal/2.1.1/bin:/home/stats/wolfch/my_libs:/home/stats/wolfch/R_libs/3.5.1:/cm/shared/apps/R/3.5.1/lib64/R:/home/stats/wolfch/wood/curl/wolf_curl/lib:/cm/shared/apps/R/3.2.2/lib64/R:/cm/shared/apps/slurm/14.11.11/lib64/slurm:/cm/shared/apps/slurm/14.11.11/lib64:/cm/shared/apps/sge/2011.11p1/lib/linux-x64:/cm/local/apps/gcc/5.1.0/lib:/cm/local/apps/gcc/5.1.0/lib64
export LD_RUN_PATH=/cm/shared/apps/gdal/2.1.1/bin
export PKG_CONFIG_PATH=/home/stats/wolfch/lib/pkgconfig
export R_LIBS=/home/stats/wolfch/R_libs/3.5.1

/cm/shared/apps/R/3.5.1/bin/R

require(ggplot2)
require(RColorBrewer)
require(scales)
require(raster)
require(rgdal)
require(sf)
require(sp)
require(foreach)
require(doSNOW)
require(pbapply)
require(reshape2)
require(grid)
require(gridExtra)
require(hues)
require(plyr)
require(readxl)
require(velox)
require(unixtools)
require(pdist)
require(doMC)
require(SpaDES)

setwd("/home/stats/wolfch/Andrews")
source("scripts - Andrews/utility_functions.R")

set.tempdir("/home/stats/wolfch/temp2")

######################################## clean up vegetation rasters

rast_df = data.frame(read_excel("data/HJA variables_v3.xlsx", na="NA")) # mean(file.exists(rast_df$Path))
rast_df = rast_df[rast_df$Status == "In use" & rast_df$Group != "PC",]

# the 5-m rasters need to be aligned - cleanest to use aggregated extent of a 1-m raster, I think
# (so we can disaggregate to 1-m later, if needed)
elev = raster("data/LIDAR/GI010/gi01001.e00")
elev_5 = aggregate(elev, fact=5)

res_1 = stack(foreach(i=which(rast_df$Resolution_m == 1 & rast_df$Group == "Vegetation"), .packages=c("raster")) %do% {
        out = raster(rast_df$Path[i])
        names(out) = rast_df$Variable[i]
        return(out)
})

cl = makeCluster(8, type = "SOCK", outfile="")
registerDoSNOW(cl)
clusterEvalQ(cl, {unixtools::set.tempdir("/home/stats/wolfch/temp2")})
res_5 = stack(foreach(i=which(rast_df$Resolution_m == 5 & rast_df$Group == "Vegetation"), .packages=c("raster")) %dopar% {
        print(i)
        rast = raster(rast_df$Path[i]) # some (e.g. mean_height) have bad pixels near edges.. buffer NA inward to be safe.
        drop = rast_df$Drop_cells[i]
        if(drop > 0){
                rast_NA = rast
                rast_NA[] = as.numeric(is.na(rast[]))
                rast_clump = clump(rast_NA, directions=8) # NA and 0 are bg values for clumping
                sort(table(rast_clump[]))
                outside = setValues(rast, 0) # had a bug - was "rast * 0"
                outside[rast_clump[] %in% names(which(table(rast_clump[]) > 1e4))] = NA
                for(x in 1:drop){ # drop outer cells this many times
                        edge = boundaries(outside, type='inner', classes=FALSE, directions=8)
                        outside[edge[] %in% 1] = NA
                } # expands outside inward by Drop_cells # of cells
                rast[is.na(outside)] = NA
                # we can also have errors along borders at the extents of the raster -- just drop these parts directly:
                for(drop_row in c(1:drop, nrow(rast):(nrow(rast)-drop+1))) rast[drop_row,] = NA
                for(drop_col in c(1:drop, ncol(rast):(ncol(rast)-drop+1))) rast[,drop_col] = NA                
        }
        out = projectRaster(rast, elev_5, method="ngb")
        names(out) = rast_df$Variable[i]
        return(out)
})
stopCluster(cl)

cl = makeCluster(3, type = "SOCK", outfile="") # doing all at once takes too much memory (could request more..)
registerDoSNOW(cl)
clusterEvalQ(cl, {unixtools::set.tempdir("/home/stats/wolfch/temp2")})
res_1_agg = stack(foreach(i=which(rast_df$Resolution_m == 1 & rast_df$Group == "Vegetation"), .packages=c("raster")) %dopar% {
        print(i)
        out = aggregate(raster(rast_df$Path[i]), fact=5) # could switch to velox for speed
        names(out) = rast_df$Variable[i]
        return(out)
})
stopCluster(cl)

res_5_joined = stack(res_1_agg, res_5) # can use for gridded prediction, mapping...
res_5_joined[["PC1"]] = raster("data_processed/veg_index_rasts/all_PC/PC_1.tif")
res_5_joined[["PC2"]] = raster("data_processed/veg_index_rasts/all_PC/PC_2.tif")

saveRDS(res_5_joined, "data_processed/res_5_joined_veg.RDS")
res_5_joined = readRDS("data_processed/res_5_joined_veg.RDS")

######################################## get PC1 from Sarah data (to add to reduce df in next step)

sites = read.csv("data/temperature/locMS045.csv", as.is=TRUE)[-(1:4),] # first few are for regions
sites_Sarah = read_excel("data/Sarah/PC_point.xlsx")
PC_data = read_excel("data/Sarah/PCA_results_Frey_etal_SciAdv.xlsx")

sites_Sarah = st_as_sf(sites_Sarah, coords=c("X","Y"),crs=as.character(crs(raster("data/LIDAR/GI010/gi01001.e00"))))
sites_Sarah_WGS84 = st_transform(sites_Sarah, crs=4326)
boundaries = read_sf("data/boundaries/boundaries_WGS84.shp")

pdf("temp.pdf")
ggplot(boundaries) +
        geom_sf() +
        geom_sf(data=sites_Sarah_WGS84)
dev.off() # looks like we have the projection correct for Sarah's data.

sites_WGS84 = st_as_sf(sites, coords=c("WEST_BOUND_COORD_decdeg","NORTH_BOUND_COORD_decdeg"),crs=4326)
sites_UTM_spatial = st_transform(sites_WGS84,crs=st_crs(projection(raster("data/LIDAR/GI010/gi01001.e00"))), asText=TRUE)

a = st_coordinates(sites_Sarah)
b = st_coordinates(sites_UTM_spatial)

pdf("temp.pdf")
ggplot() +
        geom_sf(data=sites_WGS84, col="red", size=3) +
        geom_sf(data=sites_Sarah_WGS84, col="blue")
dev.off() # appears that one Sarah site is not in the main temp. dataset... probably just won't use

sites_Sarah$LOCATION_CODE = sapply(1:nrow(sites_Sarah), function(i){
        dists = pdist(st_coordinates(sites_Sarah[i,,drop=FALSE]), st_coordinates(sites_UTM_spatial))
        print(min(dists@dist)) # looks good
        if(min(dists@dist) > 1) return(NA)
        matched = trimws(sites_UTM_spatial$LOCATION_CODE[which(dists@dist < 1)])
})

PC_data$LOCATION_CODE = sites_Sarah$LOCATION_CODE[match(PC_data$point, sites_Sarah$POINT)]
sum(is.na(PC_data$LOCATION_CODE)) # perhaps not a real Sarah site??

sites_Sarah_WGS84$PC1 = PC_data$PC1[match(sites_Sarah_WGS84$POINT, PC_data$point)]

source("https://raw.githubusercontent.com/paleolimbot/ggspatial/9d8e671379f3373235d321ddd8abe8b5f6d14ca4/R/annotation-scale.R") # couldn't install package... w/e

p = ggplot() +
        geom_sf(data=boundaries) +
        geom_sf(data=sites_Sarah_WGS84[! is.na(sites_Sarah_WGS84$PC1),], aes(fill=PC1), size=2.75, shape=21, color="black") +
        scale_fill_gradientn(colors=brewer.pal(11,"RdYlGn"),
                              breaks=pretty_breaks(n=5),
                              guide=guide_colorbar(title="Forest structure gradient (PC1)",
                                                   title.position = "top")) +
        theme_bw() +
        theme(axis.ticks=element_blank(),
                      axis.title=element_blank(),
                      axis.text=element_blank(),
                      panel.grid.major = element_line(colour="transparent"),
                      panel.grid.minor = element_line(colour="transparent"),
                      legend.position=c(0.05,0.95),
                      legend.justification=0:1,
                      legend.direction="horizontal",
                      legend.key.width=unit(2,"lines")) +
#        annotation_scale(location = "tl", pad_y = unit(2.5, "cm"), pad_x = unit(1, "cm"))
        annotation_scale(location = "br", width_hint=0.3)

pdf("output/veg_PC/PC_1_map.pdf", width=7, height=5)
p
dev.off()

######################################## reduce by sites

sites = read.csv("data/temperature/locMS045.csv", as.is=TRUE)
weather_sites = read.csv("data/weather_station/locMS001.csv", as.is=TRUE)
sites = rbind(sites, weather_sites)
sites$is_station = sites$LOCATION_CODE %in% weather_sites$LOCATION_CODE
sites = sites[which(sites$WEST_BOUND_COORD_decdeg == sites$EAST_BOUND_COORD_decdeg),] # others are regional
sites$LOCATION_CODE = trimws(sites$LOCATION_CODE)
sites_WGS84 = st_as_sf(sites, coords=c("WEST_BOUND_COORD_decdeg","NORTH_BOUND_COORD_decdeg"),crs=4326)
sites_UTM_spatial = st_transform(sites_WGS84,crs=st_crs(projection(raster("data/LIDAR/GI010/gi01001.e00"))), asText=TRUE)
sites_UTM_sp = as(sites_UTM_spatial, 'Spatial')
sites_UTM = data.frame(sites_UTM_sp, st_coordinates(sites_UTM_spatial))
sites_buffer = st_buffer(sites_UTM_spatial, 500) # for cropping

micro_files = setdiff(list.files("data/Corrected_LIDAR_Dave/microtopography"),"elev_5")
microtopo = stack(pblapply(micro_files, function(scale){
        stack = readRDS(paste0("data/Corrected_LIDAR_Dave/microtopography/", scale))
        names(stack) = paste0(names(stack), "XX", scale)
        return(stack)
}))
elev_5 = readRDS("data/Corrected_LIDAR_Dave/microtopography/elev_5")
names(elev_5) = "ElevationXX5"
topo = stack(microtopo, elev_5)

cl = makeCluster(8, type = "SOCK", outfile="")
registerDoSNOW(cl)
clusterEvalQ(cl, {
        unixtools::set.tempdir("/home/stats/wolfch/temp2")
        require(foreach); require(raster); require(reshape2); require(sf); # faster to load here
        # note: sf required for cropping w/ sf object!
})
rast_reduce = foreach(i = 1:nrow(sites_buffer), .packages=c("raster","tidyr"), .combine="rbind") %dopar% {
        print(i)
        # veg reduce        
        rast_crop_5 = crop(res_5_joined, sites_buffer[i,])
        dist_rast_5 = distanceFromPoints(rast_crop_5, sites_UTM_sp[i,])
        reduce_5 = t(sapply(seq(10,100,by=10), function(x){
                apply(rast_crop_5[dist_rast_5 <= x], 2, function(x) mean(x,na.rm=TRUE))
        }))
        rownames(reduce_5) = seq(10,100,by=10)
        reduce_5 = melt(reduce_5)
        out_veg = data.frame(site=sites$LOCATION_CODE[i], scale=reduce_5$Var1, variable=reduce_5$Var2,
                             value=reduce_5$value)
        # topo reduce
        topo_crop = crop(topo, sites_buffer[i,])
        dist_rast_topo = distanceFromPoints(topo_crop, sites_UTM_sp[i,])
        reduce_topo = melt(raster::extract(topo_crop, sites_UTM_sp[i,]))
        reduce_topo = separate(reduce_topo,"Var2",c("variable","scale"), sep="XX")
        out_topo = data.frame(site=sites$LOCATION_CODE[i], reduce_topo[,c("variable","scale","value")],
                              stringsAsFactors=FALSE)
        #
        out = na.omit(rbind(out_veg, out_topo))
        saveRDS(out, paste0("/home/stats/wolfch/temp3/",i))
        return(out)
}
stopCluster(cl)

rast_reduce = do.call("rbind", lapply(1:nrow(sites_buffer),function(i) readRDS(paste0("/home/stats/wolfch/temp3/",i))))
# rast_reduce$variable = factor(rast_reduce$variable, levels=rast_df$Variable)
rast_reduce$scale = as.numeric(rast_reduce$scale)

# add PC1...
#temp = data.frame(site=PC_data$LOCATION_CODE, scale=NA, variable="PC1", value=PC_data$PC1)
#rast_reduce = rbind(rast_reduce, temp)
saveRDS(rast_reduce, "data_processed/rast_reduce_v2.RDS")

rast_reduce = readRDS("data_processed/rast_reduce_v2.RDS")
site_colors = as.character(iwanthue(length(unique(rast_reduce$site)), 0, 360, 54, 180, 27, 67))

p1 = ggplot(rast_reduce[rast_reduce$variable != "Elevation",], aes(x=scale,y=value,group=site,color=site)) +
        facet_wrap(~ variable, scales="free_y") +
        geom_line(alpha=0.25, show.legend=FALSE) +
        scale_color_manual(values=site_colors) +
        theme_bw() +
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              panel.border=element_rect(color="black")) +
        xlab("Buffer radius (m)") +
        ylab("Mean value") +
        scale_x_log10() +
        geom_vline(xintercept=c(25,250),color="black",linetype="dashed") +
        annotation_logticks(sides="b")
png("output/scale_plot_1.png", width=7.5, height=7.5, units="in", res=300)
p1
dev.off()

######################################## map of predictors (5 m res.)

sites = read.csv("data/temperature/locMS045.csv", as.is=TRUE)
sites = sites[which(sites$WEST_BOUND_COORD_decdeg == sites$EAST_BOUND_COORD_decdeg),] # others are regional
sites$LOCATION_CODE = trimws(sites$LOCATION_CODE)
sites_WGS84 = st_as_sf(sites, coords=c("WEST_BOUND_COORD_decdeg","NORTH_BOUND_COORD_decdeg"),crs=4326)
sites_UTM_spatial = st_transform(sites_WGS84,crs=st_crs(projection(raster("data/LIDAR/GI010/gi01001.e00"))), asText=TRUE)
sites = data.frame(sites_UTM_spatial, st_coordinates(sites_UTM_spatial))
sites$geometry = NULL

res_5_joined = readRDS("data_processed/res_5_joined_veg.RDS")
stack = readRDS("data/Corrected_LIDAR_Dave/microtopography/10") # maybe 5m is better but w/e
pred_mask = readRDS("data_processed/pred_mask.RDS")
stack[pred_mask == 0] = NA
res_5_joined = stack(res_5_joined, stack)

rast_df$Group_color = muted(c("blue","brown","green"))[factor(rast_df$Group)]
rast_df$Group_linetype = c("dotted","dashed","solid")[factor(rast_df$Group)]

# note sure why this point shape has stroke http://www.sthda.com/english/wiki/ggplot2-point-shapes
# maybe try 19?...

elev_5 = readRDS("data/Corrected_LIDAR_Dave/microtopography/elev_5")
names(elev_5) = "Elevation"
res_5_plot = addLayer(res_5_joined, Elevation = elev_5)

plot_list = lapply(setdiff(rast_df$Variable,c("PC1","PC2")), function(var){
        print(var)                           
        rast = res_5_plot[[var]]
        rast_pts <- data.frame(rasterToPoints(rast))
        # rast_pts = rast_pts[sample(1:nrow(rast_pts),1000),]
        colnames(rast_pts) <- c('x','y','value')

        if(var == "Position_index"){
                rast_pts$value[rast_pts$value > 1] = 1
                rast_pts$value[rast_pts$value < -1] = -1
        }
        if(var == "Density_0_2"){
                rast_pts$value[rast_pts$value < 95] = 95
        }

        p = ggplot(data=rast_pts,aes(x=x,y=y,fill=value)) +
           geom_raster() +
           coord_equal() +
           geom_point(data=sites, aes(x=X,y=Y,fill=NULL), shape=20, size=0.75, stroke=0) +
           scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),
                                breaks=pretty_breaks(n=3),                                
                                guide=guide_colorbar(title=format_names(var))) +
           theme_bw() +
           theme(axis.ticks=element_blank(),
                      axis.title=element_blank(),
                      axis.text=element_blank(),
                      panel.grid.major = element_line(colour="transparent"),
                      panel.grid.minor = element_line(colour="transparent"),
                      panel.border = element_rect(color=rast_df$Group_color[rast_df$Variable == var],
                                                  # linetype=rast_df$Group_linetype[rast_df$Variable == var],
                                                  size=2),
                      plot.margin=margin(0,1,0,0),
                      legend.position=c(0.02,0.98),
                      legend.justification=0:1,
                      legend.title=element_text(size=8),
                      legend.text=element_text(size=7),
                      legend.key.height=unit(0.4,"lines"),
                      legend.key.width=unit(0.6,"lines"),
                      legend.background=element_rect(fill=NA))
})

# for speed could write individual maps to pngs and then tile (or just resample to 10-m res. or so -- zooming
# in on the final image and comparing to rasters in Arc suggests we aren't getting 5-m detail)
png("output/covar_map.png", width=7.5, height=7.75, units="in", res=200)
do.call("grid.arrange", c(plot_list, nrow=6))
dev.off()

# NOT RUN
pdf("output/covar_map.pdf", width=7.5, height=7.75)
do.call("grid.arrange", c(plot_list, nrow=6))
dev.off()

######################################## smooth veg. rasters for spatial prediction (and join with elev and microtopo)

rast_stack = readRDS("data_processed/res_5_joined_veg.RDS")
# rast_stack = dropLayer(rast_stack, setdiff(names(rast_stack), rast_df$Variable[rast_df$Group != "Microtopography"]))

scales = seq(10,100,by=10)
wt_list = lapply(scales, function(scale){
        wt  = raster:::.circular.weight(raster(vals=1), scale/5)
        wt[wt > 0] = 1 # velox multiples by weight matrix (elementwise) and sums
        return(wt)
})                
names(wt_list) = as.character(scales)

registerDoMC(8) # faster than usual makeCluster SOCK for single node (fork shares memory, I think)
# yep: https://www.glennklockwood.com/data-intensive/r/foreach-parallelism.html

smoothed_stack = stack(lapply(names(rast_stack), function(rast_name){                                      
        print(rast_name)
#        if(file.exists(paste0("data_processed/rasters/", rast_name)))
#                return(readRDS(paste0("data_processed/rasters/", rast_name)))
        # best when num. tiles is a small multiple of num. cores                                      
        rast_split = splitRaster(rast_stack[[rast_name]], 4, 6, buffer=c(max(scales)/5,max(scales)/5))
        smoothed_rast_split = foreach(rast_small=rast_split, .packages=c("velox","raster"),
                                      .export=c("wt_list","scales")) %dopar% {
                rast_small_not_NA = rast_small
                rast_small_not_NA[] = 1 - is.na(rast_small[])
                rast_small[is.na(rast_small[])] = 0 # set NA to 0 so velox sum will work (tracking NAs in rast_small_not_NA)
                rast_small_vel = velox(rast_small); rast_small_not_NA_vel = velox(rast_small_not_NA)                
                #
                out = lapply(scales, function(scale){
                        print(scale)                                     
                        r_sum = rast_small_vel$copy()
                        NA_sum = rast_small_not_NA_vel$copy()
                        r_sum$sumFocal(weights=wt_list[[as.character(scale)]])
                        NA_sum$sumFocal(weights=wt_list[[as.character(scale)]])
                        out = r_sum$as.RasterStack()[[1]] / NA_sum$as.RasterStack()[[1]]
                        names(out) = paste0(rast_name, "XX", scale)
                        return(out)
                })
                return(out)
        }
        smoothed_rast_split = purrr::transpose(smoothed_rast_split)
        out = stack(foreach(rast_list=smoothed_rast_split, .packages=c("raster", "SpaDES")) %dopar% {
                # setPaths(cachePath="temp", inputPath="temp", modulePath="temp", outputPath="temp", silent = FALSE)                    
                out = mergeRaster(rast_list)
                return(out)
        })
        saveRDS(out, paste0("data_processed/rasters/", rast_name))
        return(out)
}))
# great: https://stackoverflow.com/questions/16384140/how-to-avoid-duplicating-objects-with-foreach
# https://stackoverflow.com/questions/18028452/reading-global-variables-using-foreach-in-r
# https://stackoverflow.com/questions/45767416/how-to-export-many-variables-and-functions-from-global-environment-to-foreach-lo
# I think we can do better... nest foreach inside an outer foreach, set .noexport=TRUE, use iter over rast list or w/e, ...

micro_files = setdiff(list.files("data/Corrected_LIDAR_Dave/microtopography"),"elev_5")
microtopo = stack(pblapply(micro_files, function(scale){
        stack = readRDS(paste0("data/Corrected_LIDAR_Dave/microtopography/", scale))
        names(stack) = paste0(names(stack), "XX", scale)
        return(stack)
}))
elev_5 = readRDS("data/Corrected_LIDAR_Dave/microtopography/elev_5")
names(elev_5) = "ElevationXX5"
predictors = stack(dropLayer(smoothed_stack, grep("Elevation", names(smoothed_stack), value=TRUE)),
                   microtopo,
                   elev_5)

saveRDS(predictors, "data_processed/predictors.RDS")

