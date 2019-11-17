my_pred_ALE = function(X.model, newdata){
        predict(X.model,as.matrix(newdata))
}

##############################  set up data

scales = c(5, 10) # 5 is for elevation

rast_df = data.frame(read_excel("data_input/HJA variables_final.xlsx", na="NA"))
rast_df = rast_df[rast_df$Group != "Vegetation",]
rast_df$Variable_color = c(brewer.pal(4,"Greens")[c(3,4)],
                           brewer.pal(9,"Blues")[9],
                           brewer.pal(table(rast_df$Group)["Microtopography"],"YlOrBr"))

rast_reduce = readRDS("data_processed/rast_reduce.RDS")
rast_reduce = rast_reduce[rast_reduce$variable %in% rast_df$Variable & rast_reduce$scale %in% scales,]
rast_reduce$var_scale = paste0(rast_reduce$variable, "XX", rast_reduce$scale)
rast_spread = spread(rast_reduce[,c("site","var_scale","value"),],key="var_scale",value="value")

temp_merged = readRDS("data_processed/temperature_metrics.RDS")

data = merge(temp_merged, rast_spread, by="site")
data$site_year = paste(data$LOCATION_CODE, data$Year)

############################## fit models!

pred_mat = expand.grid(var=unique(data$variable), response=c("value","delta_metrics"), stringsAsFactors=FALSE)

partial_mat = data.frame(rbindlist(foreach(i=1:nrow(pred_mat)) %dopar%{
        print(i)
        var = pred_mat$var[i]; response=pred_mat$response[i];
        data_var = data[data$variable==var,]
        data_var$value = data_var[,response]        
        covar_data = data_var[! is.na(data_var$value),]

        mod = lightgbm(data = as.matrix(covar_data[,unique(rast_reduce$var_scale)]),
                label = covar_data$value,
                num_leaves = round(0.75*2^5),
                learning_rate = 1,
                objective = "quantile",
                alpha = 0.9,
                nrounds = 25,
                nthread = 1)
        df_train = na.omit(as.matrix(covar_data[,unique(rast_reduce$var_scale)])) # for ALEPlot (can't handle NAs apparently)

        out = data.frame(rbindlist(lapply(unique(rast_reduce$var_scale), function(feature){
                print(feature)
                ALE_data = suppressGraphics({ALEPlot(df_train, mod, pred.fun=my_pred_ALE,
                                   J=which(colnames(df_train) == feature), K=40, NA.plot = TRUE)})
                out = data.frame(pred_mat[i,,drop=FALSE], variable=feature, x=ALE_data$x.values, yhat=ALE_data$f.values)
                out$yhat = out$yhat - out$yhat[1]
                out$scaled_x = (out$x - min(out$x))/(max(out$x) - min(out$x)) # could instead use ecdf(..)(..) on training data
                return(out)
        })))
}))
partial_mat$variable = str_split_fixed(partial_mat$variable,"XX",2)[,1]
partial_mat$variable = factor(partial_mat$variable, levels=rast_df$Variable, labels=format_names(rast_df$Variable))
partial_mat$response_long = factor(partial_mat$response, levels=c("value","delta_metrics"),
                           labels=c("Unadjusted","Relative to free-air"))
partial_mat$var = format_names(partial_mat$var)

p = ggplot(partial_mat,
                aes(x=scaled_x,y=yhat,color=variable)) +
         geom_line(size=0.5) +
         facet_grid(var ~ response_long, scales="free") +
         scale_color_manual(values = rast_df$Variable_color,
                            guide=guide_legend(ncol=1,title=NULL)) +   
        geom_hline(yintercept=0, linetype="dashed") +
         theme_bw() +
         ylab("Change in average prediction") +
         scale_x_continuous(breaks=c(0,0.5,1)) +
         xlab("Scaled predictor") +
         theme(axis.ticks=element_line(color="black"),
               axis.text=element_text(color="black"),
               panel.border=element_rect(color="black"),
               panel.spacing.x=unit(0.6,"lines"))

png("output/quantile/ALE_plot.png", width=7.75, height=7.75, units="in", res=400)
print(p)
dev.off()

pdf("output/quantile/ALE_plot.pdf", width=7.75, height=7.75)
print(p)
dev.off()

