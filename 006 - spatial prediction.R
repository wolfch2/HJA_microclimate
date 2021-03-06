############################## set up data

rast_reduce = readRDS("data_processed/rast_reduce.RDS")
rast_reduce = rast_reduce[! rast_reduce$variable %in% c("PC1","PC2"),] # not using PC vars for main models
rast_reduce$var_scale = paste0(rast_reduce$variable, "XX", rast_reduce$scale)
rast_spread = spread(rast_reduce[,c("site","var_scale","value"),],key="var_scale",value="value")

temp_merged = readRDS("data_processed/temperature_metrics.RDS")
temp_merged$site_year = paste(temp_merged$LOCATION_CODE, temp_merged$Year)

data = merge(temp_merged, rast_spread, by="site")

temp_sd = aggregate(cbind(value,delta_metrics) ~ site+variable,FUN=sd,temp_merged)
data_sd = merge(temp_sd, rast_spread, by="site")

############################## build predictor matrix for newdata

predictors = readRDS("data_processed/predictors.RDS")
predictor_mat = as.matrix(predictors)
pred_mask = readRDS("data_processed/pred_mask.RDS")

############################## fit models!

pred_mat = expand.grid(var=unique(data$variable),
                       response=c("value","delta_metrics"),
                       type=c("quantile","sd"),
                       stringsAsFactors=FALSE)

pred_rast_list = foreach(i=1:nrow(pred_mat)) %do%{  # startup is slow since predictor_mat big
        print(i)
        gc() # seems to help w/ memory allocation errors
        var = pred_mat$var[i]; response = pred_mat$response[i]; type = pred_mat$type[i]
        if (type == "quantile"){
                data_var = data[data$variable==var,]
        } else {
                data_var = data_sd[data_sd$variable==var,]
        }
        data_var$geometry = NULL
        covar_data = data_var

        # https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html
        # can likely predict faster w/ treelite https://github.com/Microsoft/LightGBM/issues/2094
        mod <- lightgbm(data = as.matrix(covar_data[,unique(rast_reduce$var_scale)]),
                label = covar_data[,response],
                num_leaves = round(0.75*2^5),
                learning_rate = 1,
                objective = ifelse(type=="quantile","quantile","regression"),
                alpha = 0.9,
                nrounds = 25,
                nthread = 8) # nthread > 1 also enables parallel prediction!

        pred = predict(mod, data = predictor_mat[,unique(rast_reduce$var_scale)])

        out = pred_mask
        out[] = pred
        out[pred_mask[] == 0] = NA
        return(out)
}

pbsapply(unique(pred_mat$var), function(var){ # check correlations
  cor(na.omit(as.matrix(stack(pred_rast_list[which(pred_mat$type == "quantile" & pred_mat$var == var)]))))[1,2]
})

############################## plot unadjusted

temps = readRDS("data_processed/temperature_metrics.RDS")
sites = readRDS("data_processed/sites.RDS")
sites= merge(temps, sites, by.x="SITECODE", by.y="LOCATION_CODE")

map_theme = theme(axis.ticks=element_blank(),
                      axis.title=element_blank(),
                      axis.text=element_blank(),
                      panel.grid.major = element_line(colour="transparent"),
                      panel.grid.minor = element_line(colour="transparent"),
                      plot.margin=margin(0,1,0,0),
                      legend.position=c(0.02,0.98),
                      legend.justification=0:1,
                      legend.title=element_text(size=6.5),
                      legend.text=element_text(size=6.5),
                      legend.key.height=unit(0.4,"lines"),
                      legend.key.width=unit(0.6,"lines"),
                      legend.background=element_rect(fill=NA)) 

plot_list = lapply(sort(unique(pred_mat$var))[c(1:3,5:6,4)], function(var){
        print(var)
        rast = pred_rast_list[[which(pred_mat$var == var & pred_mat$response == "value" & pred_mat$type == "quantile")]]
        rast = aggregate(rast, 2) # optional - aggregate to reduce file size

        rast_pts <- data.frame(rasterToPoints(rast))
        colnames(rast_pts) <- c('x','y','value')

        sites_var = sites[sites$variable == var,]
        sites_var = sites_var[! duplicated(sites_var$LOCATION_CODE),]
        print(nrow(sites_var))

        p = ggplot(data=rast_pts,aes(x=x,y=y,fill=value)) +
           geom_raster() +
           coord_equal() +
           geom_point(data=sites_var, aes(x=X,y=Y,fill=NULL), shape=20, size=0.5, stroke=0) +
           scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),
                                breaks=pretty_breaks(n=3),                                
                                guide=guide_colorbar(title=format_names(var))) +
           theme_bw() +
           map_theme +
           expand_limits(y=min(sites$Y) - 300)

        if(var == "Apr_Jun_mean_max")
                p = p + annotation_scale(location = "bl", pad_x=unit(0.5,"in"), pad_y=unit(0.05, "in"), width_hint=0.25, height=unit(0.2,"lines"), text_cex=0.5)

        return(p)
})

elevation = aggregate(predictors[["ElevationXX5"]], 2)
elevation_pts <- data.frame(rasterToPoints(elevation))
colnames(elevation_pts) <- c('x','y','value')
sites_var = sites[! duplicated(sites$LOCATION_CODE),]

p_elev = ggplot(data=elevation_pts,aes(x=x,y=y,fill=value)) +
           geom_raster() +
           coord_equal() +
           geom_point(data=sites_var, aes(x=X,y=Y,fill=NULL), shape=20, size=0.5, stroke=0) +
           scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),
                                breaks=pretty_breaks(n=3),                                
                                guide=guide_colorbar(title="Elevation (m)")) +
           theme_bw() +
           map_theme +
           expand_limits(y=min(sites$Y) - 300)

PC1 = aggregate(predictors[["PC1XX10"]], 2)
PC1_pts <- data.frame(rasterToPoints(PC1))
colnames(PC1_pts) <- c('x','y','value')
sites_var = sites[! duplicated(sites$LOCATION_CODE),]

p_PC1 = ggplot(data=PC1_pts,aes(x=x,y=y,fill=value)) +
           geom_raster() +
           coord_equal() +
           geom_point(data=sites_var, aes(x=X,y=Y,fill=NULL), shape=20, size=0.5, stroke=0) +
           scale_fill_gradientn(colors=brewer.pal(9,"YlGn"),
                                breaks=pretty_breaks(n=3),                                
                                guide=guide_colorbar(title="Forest structure (PC1)")) +
           theme_bw() +
           map_theme +
           expand_limits(y=min(sites$Y) - 300)

pred_plot = plot_grid(plot_grid(plotlist=plot_list, nrow=2), plot_grid(p_elev, p_PC1, nrow=1), nrow=2, rel_heights=c(1,0.8))

png("output/quantile/pred_map_90.png", width=6-0.2, height=7/4*7.75/6*2, units="in", res=300)
plot(pred_plot)
dev.off()

pdf("output/quantile/pred_map_90.pdf", width=6-0.2, height=7/4*7.75/6*2)
plot(pred_plot)
dev.off()

############################## plot adjusted

plot_list = lapply(sort(unique(pred_mat$var))[c(1:3,5:6,4)], function(var){
        print(var)
        rast = pred_rast_list[[which(pred_mat$var == var & pred_mat$response == "delta_metrics" & pred_mat$type == "quantile")]]
        rast = aggregate(rast, 2) # optional - aggregate to reduce file size

        rast_pts <- data.frame(rasterToPoints(rast))
        colnames(rast_pts) <- c('x','y','value')

        sites_var = sites[sites$variable == var,]
        sites_var = sites_var[! duplicated(sites_var$LOCATION_CODE),]
        print(nrow(sites_var))

        p = ggplot(data=rast_pts,aes(x=x,y=y,fill=value)) +
           geom_raster() +
           coord_equal() +
           geom_point(data=sites_var, aes(x=X,y=Y,fill=NULL), shape=20, size=0.5, stroke=0) +
           scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),
                                breaks=pretty_breaks(n=3),                                
                                guide=guide_colorbar(title=format_names(var))) +
           theme_bw() +
           map_theme +
           expand_limits(y=min(sites$Y) - 300)

        if(var == "Apr_Jun_mean_max")
                p = p + annotation_scale(location = "bl", pad_x=unit(0.5,"in"), width_hint=0.25, height=unit(0.2,"lines"), text_cex=0.5)

        return(p)
})

png("output/quantile/pred_map_90_delta_metrics.png", width=6-0.2, height=7.75/6*2, units="in", res=400)
print(do.call("grid.arrange", c(plot_list, nrow=2)))
dev.off()

############################## s.d. maps

sd_raw = lapply(sort(unique(pred_mat$var))[c(1:3,5:6,4)], function(var){
        print(var)
        rast = pred_rast_list[[which(pred_mat$var == var & pred_mat$response == "value" & pred_mat$type == "sd")]]
        rast = aggregate(rast, 2) # optional - aggregate to reduce file size

        rast_pts <- data.frame(rasterToPoints(rast))
        colnames(rast_pts) <- c('x','y','value')

        sites_var = sites[sites$variable == var,]
        sites_var = sites_var[! duplicated(sites_var$LOCATION_CODE),]
        print(nrow(sites_var))

        p = ggplot(data=rast_pts,aes(x=x,y=y,fill=value)) +
           geom_raster() +
           coord_equal() +
           geom_point(data=sites_var, aes(x=X,y=Y,fill=NULL), shape=20, size=0.5, stroke=0) +
           scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),
                                breaks=pretty_breaks(n=3),                                
                                guide=guide_colorbar(title=format_names(var))) +
           theme_bw() +
           map_theme +
           expand_limits(y=min(sites$Y) - 300)

        if(var == "Apr_Jun_mean_max")
                p = p + annotation_scale(location = "bl", pad_x=unit(0.5,"in"), width_hint=0.25, height=unit(0.2,"lines"), text_cex=0.5)

        return(p)
})

sd_adj = lapply(sort(unique(pred_mat$var))[c(1:3,5:6,4)], function(var){
        print(var)
        rast = pred_rast_list[[which(pred_mat$var == var & pred_mat$response == "delta_metrics" & pred_mat$type == "sd")]]
        rast = aggregate(rast, 2) # optional - aggregate to reduce file size

        rast_pts <- data.frame(rasterToPoints(rast))
        colnames(rast_pts) <- c('x','y','value')

        sites_var = sites[sites$variable == var,]
        sites_var = sites_var[! duplicated(sites_var$LOCATION_CODE),]
        print(nrow(sites_var))

        p = ggplot(data=rast_pts,aes(x=x,y=y,fill=value)) +
           geom_raster() +
           coord_equal() +
           geom_point(data=sites_var, aes(x=X,y=Y,fill=NULL), shape=20, size=0.5, stroke=0) +
           scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),
                                breaks=pretty_breaks(n=3),                                
                                guide=guide_colorbar(title=format_names(var))) +
           theme_bw() +
           map_theme +
           expand_limits(y=min(sites$Y) - 300)

        return(p)
})

top = annotate_figure(plot_grid(plotlist=sd_raw, nrow=2), fig.lab="A. unadjusted")
bot = annotate_figure(plot_grid(plotlist=sd_adj, nrow=2), fig.lab="B. Adjusted")

p = plot_grid(ggdraw()+draw_label("A. Unadjusted temperature standard deviation"), plot_grid(plotlist=sd_raw, nrow=2),
          ggdraw()+draw_label("B. Adjusted temperature standard deviation"), plot_grid(plotlist=sd_adj, nrow=2),
          ncol = 1, rel_heights=c(1,20,1,20))

png("output/sd_maps.png", width=6-0.2, height=8.8/6*2*2, units="in", res=400)
print(p)
dev.off()
