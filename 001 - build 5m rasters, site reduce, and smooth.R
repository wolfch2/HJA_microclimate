rast_df = data.frame(read_excel("data_input/HJA variables_final.xlsx", na="NA"))

######################################## set up elevation raster (buffer inward 30/5=6 pixels) and prediction mask

# the 5-m rasters need to be aligned - cleanest to use aggregated extent of a 1-m raster, I think
# (so we can disaggregate to 1-m later, if needed)
elev = raster(rast_df$Path[rast_df$Variable == "Elevation"])

elev_5 = aggregate(elev, fact=5)
elev_5 = buffer_inward(elev_5, rast_df$Drop_cells[rast_df$Variable == "Elevation"] / 5)

pred_mask = elev_5
pred_mask[] = ! is.na(pred_mask[])

saveRDS(pred_mask, "data_processed/pred_mask.RDS")

######################################## clean up and resample vegetation rasters

res_1_agg = stack(foreach(i=which(rast_df$Resolution_m == 1 & rast_df$Group == "Vegetation"), .packages=c("raster")) %do% {
                          out = raster(rast_df$Path[i])
                          out = projectRaster(out, elev_5, method="ngb")        
                          names(out) = rast_df$Variable[i]
                          return(out)
})

# 51 sec on desktop, 149 sec on cluster
res_5 = stack(foreach(i=which(rast_df$Resolution_m == 5 & rast_df$Group == "Vegetation"), .packages=c("raster")) %dopar% {
                      print(i)
                      rast = raster(rast_df$Path[i]) # some (e.g. mean_height) have bad pixels near edges.. buffer NA inward to be safe.
                      drop = rast_df$Drop_cells[i]
                      if(drop > 0){
                              rast = buffer_inward(rast, drop)
                      }
                      out = projectRaster(rast, elev_5, method="ngb")
                      names(out) = rast_df$Variable[i]
                      return(out)
})

######################################## PC analysis

veg_rasts = stack(res_1_agg, res_5) # input veg rasters
predictor_mat = as.matrix(veg_rasts)

df = data.frame(predictor_mat, row_num = 1:nrow(predictor_mat))
df = na.omit(df)
row_num = df$row_num
df$row_num = NULL
df = scale(df) # since more variable variables (hehe) are not necessarily more important

PCA = prcomp(df)
apply(PCA$rotation,2,function(x) sum(x^2)) # each LC is a unit vector
PCA_mat = PCA$x

PC_rasts = pblapply(1:ncol(PCA_mat), function(i){
                            out = NA * pred_mask
                            out[row_num] = PCA_mat[,i]
                            names(out) = paste0("PC", i)
                            return(out)
})

sites = read.csv("data_input/locMS045.csv", as.is=TRUE)[-(1:4),] # first few are for regions
sites$code = trimws(sites$LOCATION_CODE)
sites$code_group = substr(sites$code,1,2)
sites$code_short = substr(sites$code,3,nchar(sites$code))
sites_WGS84 = st_as_sf(sites, coords=c("WEST_BOUND_COORD_decdeg","NORTH_BOUND_COORD_decdeg"),crs=4326)
sites_UTM_spatial = st_transform(sites_WGS84,crs=st_crs(projection(elev_5)), asText=TRUE)

PC1 = PC_rasts[[1]]
sites$PC1 = raster::extract(PC1, sites_UTM_spatial, buffer=25, fun=mean)

PC2 = PC_rasts[[2]]
sites$PC2 = raster::extract(PC2, sites_UTM_spatial, buffer=25, fun=mean)

harvest = read_sf("data_input/harvest_layer/harvest.shp")
harvest_proj = st_transform(harvest, st_crs(sites_UTM_spatial)) # original proj. info may be slightly different
harvest_union = st_union(harvest_proj)

sites$plantation = factor(st_intersects(sites_UTM_spatial, harvest_union, sparse=FALSE)[,1],
                          levels=c(TRUE,FALSE),
                          labels=c("Plantation","Mature forest/old growth"))
sites = data.frame(sites, st_coordinates(sites_UTM_spatial))

p1 = ggplot(sites, aes(x=PC1, y=PC2, color=plantation)) +
        geom_point() +
        geom_density_2d(alpha=0.15, show.legend=FALSE) +
        scale_color_manual(values=brewer.pal(3,"Set1")[3:2],
                           guide=guide_legend(title=NULL, nrow=1)) +
theme_bw() +
expand_limits(y=min(sites$PC2,na.rm=T)-0.6) +
theme(axis.text=element_text(color="black"),
      axis.ticks=element_line(color="black"),
      panel.border=element_rect(color="black"),
      legend.position=c(0,0),
      legend.key.size=unit(0.5,"lines"),
      legend.justification=c(0,0),
      aspect.ratio=1,
      legend.background=element_rect(color="black")) +
xlab("Principal component 1") +
ylab("Principal component 2")

loading_df = data.frame(PCA$rotation)
loading_df$var_name = rownames(loading_df)

p2 = ggplot(loading_df, aes(x=PC1, y=PC2)) +
        geom_hline(yintercept=0, linetype="dashed", color="gray") +
        geom_vline(xintercept=0, linetype="dashed", color="gray") +
        geom_point(color="red") +
        theme_bw() +
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              panel.border=element_rect(color="black"),
              legend.position=c(1,1),
              legend.justification=c(1,1),aspect.ratio=1,
              legend.background=element_rect(color="black")) +
geom_text_repel(aes(label=var_name), size=3, force=2) +
xlab("Principal component 1") +
ylab("Principal component 2")

prop_var = data.frame(t(summary(PCA)$importance["Proportion of Variance",,drop=FALSE]))
prop_var$PC = as.numeric(gsub("PC","",rownames(prop_var)))

p3 = ggplot(prop_var, aes(x=PC, y=Proportion.of.Variance)) +
        geom_point(color="black") +
        geom_line() +
        theme_bw() +
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              panel.border=element_rect(color="black"),
              panel.grid.minor.x=element_blank(),
              legend.position=c(1,1),
              legend.justification=c(1,1),aspect.ratio=1,
              legend.background=element_rect(color="black")) +
xlab("Principal component") +
ylab("Proportion of variance explained") +
scale_x_continuous(breaks=0:30)

rast_pts_PC1 = data.frame(rasterToPoints(PC1))
colnames(rast_pts_PC1) = c('x','y','value')

p_PC1 = ggplot(data=rast_pts_PC1,aes(x=x,y=y,fill=value)) +
        geom_raster() +
        coord_equal() +
        scale_fill_gradientn(colors=brewer.pal(9,"YlGn"),
                             guide=guide_colorbar(title="PC1")) +
theme_bw() +
theme(axis.ticks=element_blank(),
      axis.title=element_blank(),
      axis.text=element_blank(),
      panel.grid.major = element_line(colour="transparent"),
      panel.grid.minor = element_line(colour="transparent"),
      plot.margin=margin(20,1,0,15),
      legend.position=c(0.02,0.98),
      legend.justification=0:1,
      legend.background=element_rect(fill=NA)) +
geom_point(data=sites, aes(x=X,y=Y,fill=NULL,), color="black", shape=20, size=1.5, stroke=0)

rast_pts_PC2 = data.frame(rasterToPoints(PC2))
colnames(rast_pts_PC2) = c('x','y','value')

p_PC2 = ggplot(data=rast_pts_PC2,aes(x=x,y=y,fill=value)) +
        geom_raster() + #alpha=0.5) +
        coord_equal() +
        scale_fill_gradientn(colors=brewer.pal(9,"YlGn"),
                             guide=guide_colorbar(title="PC2")) +
theme_bw() +
theme(axis.ticks=element_blank(),
      axis.title=element_blank(),
      axis.text=element_blank(),
      panel.grid.major = element_line(colour="transparent"),
      panel.grid.minor = element_line(colour="transparent"),
      plot.margin=margin(20,1,0,15),
      legend.position=c(0.02,0.98),
      legend.justification=0:1,
      legend.background=element_rect(fill=NA)) +
geom_point(data=sites_XY, aes(x=X,y=Y,fill=NULL,), color="black", shape=20, size=1.5, stroke=0)

p_top = plot_grid(p2, p3, p1, labels = c('A', 'B', 'C'), label_size = 12, nrow=1, align='v')
p_bot = plot_grid(p_PC1, p_PC2, labels = c('D', 'E'), label_size = 12, nrow=1)

png("output/veg_PC/biplot.png", width=12, height=8.5, units="in", res=300)
print(grid.arrange(p_top, p_bot))
dev.off()

######################################## build microtopo vars. at multiple scales

dir.create("data_processed/microtopography/")

scales = c(10, 25, 50, 100, 250, 500)
registerDoMC(24)

# https://stackoverflow.com/questions/30904740/issue-using-saverds-with-raster-objects
foreach(scale=scales, .packages=c("raster")) %dopar% {
        print(scale)        
        elev_agg = aggregate(elev_5, fact=scale/5)

        aspect = terrain(elev_agg, opt='aspect')
        slope = terrain(elev_agg, opt='slope')
        TPI  = terrain(elev_agg, opt='TPI')
        upslope = upslope.area(elev_agg, atb=TRUE)

        rast_stack = stack(list(Aspect_eastness=sin(aspect),
                                Aspect_northness=cos(aspect),
                                Slope=slope,
                                Position_index=TPI,
                                Convergence_index=upslope[["atb"]]))

        out = resample(rast_stack, elev_5)
        saveRDS(readAll(out), paste0("data_processed/microtopography/", scale))
}

######################################## smooth vegetation rasters

rast_stack = stack(veg_rasts, PC1, PC2)

scales = seq(10,100,by=10)
wt_list = lapply(scales, function(scale){
                         wt  = raster:::.circular.weight(raster(vals=1), scale/5)
                         wt[wt > 0] = 1 # velox multiplies by weight matrix (elementwise) and sums
                         return(wt)
      })                
names(wt_list) = as.character(scales)

dir.create("data_processed/rasters/")
smoothed_stack = stack(lapply(names(rast_stack), function(rast_name){                                      
                                      print(rast_name)
                                      if(file.exists(paste0("data_processed/rasters/", rast_name)))
                                              return(readRDS(paste0("data_processed/rasters/", rast_name)))
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
                                                          out = mergeRaster(rast_list)
                                                          return(out)
                                                                    })
                                      saveRDS(out, paste0("data_processed/rasters/", rast_name))
                                      return(out)
                                }))

######################################## map of predictors (5 m res.)

stack = readRDS("data_processed/microtopography/10")
pred_mask = readRDS("data_processed/pred_mask.RDS")
stack[pred_mask == 0] = NA
names(elev_5) = "Elevation"
res_5_joined = stack(veg_rasts, elev_5, stack)

rast_df = rast_df[rast_df$Group != "PC",] # plotted separately
rast_df$Group_color = muted(c("blue","brown","green"))[factor(rast_df$Group)]

plot_list = lapply(rast_df$Variable, function(var){
                           print(var)                           
                           rast = res_5_joined[[var]]
                           rast_pts <- data.frame(rasterToPoints(rast))
                           colnames(rast_pts) <- c('x','y','value')

                           if(var == "Position_index"){
                                   rast_pts$value[rast_pts$value > 1] = 1
                                   rast_pts$value[rast_pts$value < -1] = -1
                           }
                           if(var == "Density_0_2"){
                                   rast_pts$value[rast_pts$value < 97.5] = 97.5
                           }
                           p = ggplot(data=rast_pts,aes(x=x,y=y,fill=value)) +
                                   geom_raster() +
                                   coord_equal(xlim=extent(res_5_joined)[1:2], ylim=extent(res_5_joined)[3:4]) +
                                   geom_point(data=sites, aes(x=X,y=Y,fill=NULL), shape=20, size=0.75, stroke=0) +
                                   scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),
                                                        breaks=pretty_breaks(n=3),                                
                                                        guide=guide_colorbar(title=var)) +
                           theme_bw() +
                           theme(axis.ticks=element_blank(),
                                 axis.title=element_blank(),
                                 axis.text=element_blank(),
                                 panel.grid.major = element_line(colour="transparent"),
                                 panel.grid.minor = element_line(colour="transparent"),
                                 panel.border = element_rect(color=rast_df$Group_color[rast_df$Variable == var],
                                                             size=2),
                                 plot.margin=margin(0,1,0,0),
                                 legend.position=c(0,1),
                                 legend.justification=0:1,
                                 legend.title=element_text(size=8),
                                 legend.text=element_text(size=7),
                                 legend.key.height=unit(0.4,"lines"),
                                 legend.key.width=unit(0.6,"lines"),
                                 legend.background=element_rect(fill=NA))
      })

png("output/covar_map.png", width=7.5, height=7.75, units="in", res=400)
print(do.call("grid.arrange", c(plot_list, nrow=6)))
dev.off()

########################################




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

predictor_mat = as.matrix(predictors)
pred_mask = readRDS("data_processed/pred_mask.RDS")































res_5_joined = stack(res_1_agg, res_5) # can use for gridded prediction, mapping...
res_5_joined[["PC1"]] = raster("data_processed/veg_index_rasts/all_PC/PC_1.tif")
res_5_joined[["PC2"]] = raster("data_processed/veg_index_rasts/all_PC/PC_2.tif")












pred_mask = elev_5
pred_mask[] = ! is.na(pred_mask[])
saveRDS(pred_mask, "data_processed/pred_mask.RDS") # mask for showing predicted values


























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

