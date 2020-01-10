############################## set up data

rast_reduce = readRDS("data_processed/rast_reduce.RDS")
rast_reduce = rast_reduce[! rast_reduce$variable %in% c("PC1","PC2"),] # not using PC vars for main models
rast_reduce$var_scale = paste0(rast_reduce$variable, "XX", rast_reduce$scale)
rast_spread = spread(rast_reduce[,c("site","var_scale","value"),],key="var_scale",value="value")

temp_merged = readRDS("data_processed/temperature_metrics.RDS")

data = merge(temp_merged, rast_spread, by="site")
data$site_year = paste(data$LOCATION_CODE, data$Year)

rast_df = data.frame(read_excel("data_input/HJA variables_final.xlsx", na="NA"))
rast_df = rast_df[! rast_df$Variable %in% c("PC1","PC2"),]
rast_df$Variable_color = c(colorRampPalette(brewer.pal(9,"Greens"))(table(rast_df$Group)["Vegetation"]),
                           brewer.pal(9,"Blues")[9],
                           brewer.pal(6,"YlOrBr")[(6-table(rast_df$Group)["Microtopography"]+1):6])

############################## CV w/ plot (could do for some other variables too, e.g., sd across years)

args_mat = expand.grid(var=unique(data$variable),
                       response=c("value","delta_metrics"),
                       stringsAsFactors=FALSE)

pred = foreach(arg_row=1:nrow(args_mat), .combine="rbind") %dopar%{
        var=args_mat$var[arg_row]; response=args_mat$response[arg_row];
        data_var = data[data$variable==var,]
        data_var$value = data_var[,response]
        data_var = data_var[! is.na(data_var$value),]

        data_var_split = split(data_var, data_var$Year)

        pred_CV = foreach(i = 1:length(data_var_split), .combine="rbind") %do% {
                print(i)
                test_data = data_var_split[[i]] # for annual metrics, and spatial CV, test_data will just have one row
                train_data = data_var[data_var$Year != test_data$Year[1],]
                mod = xgboost(data = as.matrix(train_data[,unique(rast_reduce$var_scale)]), label = train_data$value,
                      booster = "gbtree",objective = "reg:linear", nrounds=25,max_depth=5, verbose=0, nthread=1)
                pred = predict(mod, newdata = as.matrix(test_data[,unique(rast_reduce$var_scale)]))
                out = data.frame(var=var, obs=test_data$value, pred=pred, year=test_data$Year, site=test_data$LOCATION_CODE,
                                 response=response, stringsAsFactors=FALSE)
        }
        return(pred_CV)
}

pred$response_long = factor(pred$response, levels=c("value","delta_metrics"),
                           labels=c("Unadjusted","Relative to free-air"))

############################## year-year correlations (for paper text)

pred %>%
  filter(response == "value") %>%
  dplyr::select(-obs, -response, -response_long) %>%
  spread("year", "pred") %>%
  dplyr::select(-site) %>%
  split(.$var) %>%
  map(dplyr::select, -var) %>%
  map(cor, use="pairwise.complete.obs", method="spearman") %>%
  map(min, na.rm=TRUE)

pred %>%
  filter(response == "value") %>%
  dplyr::select(-obs, -response, -response_long) %>%
  spread("year", "pred") %>%
  dplyr::select(-site) %>%
  split(.$var) %>%
  map(dplyr::select, -var) %>%
  map(cor, use="pairwise.complete.obs", method="spearman") %>%
  map(remove_diag) %>%
  melt %>%
  select(value) %>%
  summarise(mean = mean(value,na.rm=TRUE),
            median = median(value,na.rm=TRUE),
            total = sum(! is.na(value)))
  

pred %>%
  filter(response == "value") %>%
  dplyr::select(-obs, -response, -response_long) %>%
  spread("year", "pred") %>%
  dplyr::select(-site) %>%
  split(.$var) %>%
  map(dplyr::select, -var) %>%
  map(cor, use="pairwise.complete.obs", method="spearman") %>%
  map(remove_diag) %>%
  map(mean, na.rm=TRUE) %>%
  unlist %>%
  sort

# https://stats.stackexchange.com/questions/228540/how-to-calculate-out-of-sample-r-squared
stats = do.call("rbind", lapply(split(pred, paste(pred$response,pred$var)), function(df){
        R2 = max(0, round(1 - sum((df$obs - df$pred)^2) / sum((df$obs - mean(df$obs))^2),2))
        R2 = formatC(R2, format='f', digits=2)
        MAE = formatC(round(mean(abs(df$obs - df$pred)),2), format='f', digits=2)
        out = df[1,]
        out$lab = paste0("' '*R^2*': '*", R2, "*', MAE: '*", MAE)
        return(out)
}))
pred$year = factor(pred$year); stats$year = factor(stats$year);
stats$response = factor(stats$response, levels=c("value","delta_metrics"))
stats$response_long = factor(stats$response, levels=c("value","delta_metrics"),
                           labels=c("Unadjusted","Relative to free-air"))

dummy = pred; dummy$obs = pred$pred; dummy$pred = pred$obs;

plot_list = lapply(sort(unique(pred$var))[c(1:3,5:6,4)], function(var){
        p = ggplot(pred[pred$var == var,], aes(x=obs,y=pred,color=year)) +
        geom_point(data=dummy[dummy$var == var,], color="#FFFFFF00") +
        geom_point(size=0.6, stroke=0) +
        geom_abline(intercept=0,slope=1,linetype="dashed") +
        scale_x_continuous(breaks= pretty_breaks(4)) +
        scale_y_continuous(breaks= pretty_breaks(4)) +
        facet_wrap(~ response_long, scales="free", nrow=1) +
        scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(pred$year))),
                           drop=FALSE, # important!  (don't have all years for all metrics!!)
                           guide=guide_legend(title=NULL, nrow=1,
                                              override.aes = list(size = 1.5))) +
        theme_bw() +
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              axis.title=element_blank(),
              panel.border=element_rect(color="black"),
              legend.background=element_rect(color="black"),
              strip.background=element_blank(),
              strip.text=element_blank(),
              legend.position="none",
              aspect.ratio = 1,
              plot.title = element_text(hjust = 0.5)) +
        ggtitle(format_names(var)) +
        geom_text(data=stats[stats$var == var,], aes(label=lab), x=Inf, y=-Inf, color="black", vjust=-1.15, hjust=1.05, size=2.5, parse=TRUE) +
        geom_text(data=stats[stats$var == var,], aes(label=response_long), x=-Inf, y=Inf, color="black", vjust=2.15, hjust=-0.1, size=2.5)

        return(p)
})

leg = get_legend(plot_list[[1]] + theme(legend.position="bottom"))

p_top = plot_grid(plotlist=plot_list, nrow=3, align="hv")
y.grob <- textGrob("   Predicted", gp=gpar(fontsize=12), rot=90)
x.grob <- textGrob("     Observed", gp=gpar(fontsize=12))
p_top_all = arrangeGrob(p_top, left = y.grob, bottom = x.grob)

png("output/GRIDMET/GRIDMET_CV_year.png", width=7.2, height=6.5, units="in", res=300)
print(plot_grid(p_top_all, leg, rel_heights=c(1,0.1), nrow=2))
dev.off()

pdf("output/GRIDMET/GRIDMET_CV_year.pdf", width=7.2, height=6.5)
print(plot_grid(plot_grid(plotlist=plot_list, nrow=3, align="hv"), leg, rel_heights=c(1,0.1), nrow=2))
dev.off()

stats_year = do.call("rbind", lapply(split(pred, paste(pred$response,pred$var,pred$year)), function(df){
        MAE = mean(abs(df$obs - df$pred))
        out = data.frame(df[1,], MAE=MAE, n=nrow(df))
        out$year = as.numeric(as.character(out$year))
        out$var = format_names(out$var)
        return(out)
}))

p = ggplot(stats_year, aes(x=year, y=MAE, color=response_long)) +
        geom_point(aes(size=n)) +
        scale_size(range=c(0.5,3), guide=guide_legend(title="Sites")) +
        geom_line() +
        facet_wrap(~ var, scales="free_y") +
        scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(pred$year))),
                           guide=guide_legend(title="Response",
                                              override.aes = list(size = 1.5))) +
        theme_bw() +
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              panel.border=element_rect(color="black"),
              legend.position="bottom",
              legend.background=element_rect(color="black")) +
        scale_y_continuous(limits=c(0,NA)) +
        ylab("Mean absolute error")

png("output/GRIDMET/MAE_year.png", width=7.5, height=5, units="in", res=300)
print(p)
dev.off()

############################## quantile regression models for imp. plots etc.

pred_mat = expand.grid(var=unique(data$variable), quantile=seq(0.1,0.9,by=0.1), stringsAsFactors=FALSE)
pred_mat$mod_list = lapply(1:nrow(pred_mat), function(i){
        var = pred_mat$var[i]; quantile=pred_mat$quantile[i];
        data_var = data[data$variable==var,]
        data_var$geometry = NULL
        covar_data = data_var
        mod = lightgbm(data = as.matrix(covar_data[,unique(rast_reduce$var_scale)]),
                label = covar_data$value,
                num_leaves = round(0.75*2^5),
                learning_rate = 1,
                objective = "quantile",
                alpha=quantile,
                nrounds = 25,
                nthread = 8)
        return(mod)
})

############################## pseudo-R2 (paper text)

# https://stats.stackexchange.com/questions/129200/r-squared-in-quantile-regression
lapply(unique(pred_mat$var), function(var){
  data_var = data[data$variable==var,]
  covar_data = data_var
  mod_full = lightgbm(data = as.matrix(covar_data[,unique(rast_reduce$var_scale)]),
                      label = covar_data$value,
                      num_leaves = round(0.75*2^5),
                      learning_rate = 1,
                      objective = "quantile",
                      alpha=0.9,
                      nrounds = 25,
                      nthread = 8)
  y_hat = predict(mod_full, data=as.matrix(covar_data[,unique(rast_reduce$var_scale)]))
  y_bar = rep(quantile(covar_data$value, 0.9), nrow(covar_data))
  tau = 0.9
  y = covar_data$value
  out = sum(ifelse(y >= y_hat, tau*abs(y-y_hat), (1-tau)*abs(y-y_hat)))/sum(ifelse(y >= y_bar, tau*abs(y-y_bar), (1-tau)*abs(y-y_bar)))
  names(out) = var
  return(round(out,3))
})

############################## imp. plots for other vars...

rel_inf = do.call("rbind", lapply(1:nrow(pred_mat), function(i){
        mod = pred_mat$mod_list[[i]]
        imp = data.frame(lgb.importance(mod))
        imp = rbind.fill(imp, data.frame(Feature=setdiff(unique(rast_reduce$var_scale),imp$Feature),Gain=0)) # these are omitted
        out = data.frame(pred_mat[i,], imp)
        out$mod_list = NULL
        return(out)
}))
rel_inf = separate(rel_inf,"Feature",c("Variable","Scale"),"XX")
rel_inf$Scale = as.numeric(rel_inf$Scale)
rel_inf$value = 100 * rel_inf$Gain
rel_inf$temp_var = format_names(rel_inf$var)
rel_inf$Group = rast_df$Group[match(rel_inf$Variable, rast_df$Variable)]

rel_inf_agg = aggregate(value ~ temp_var + Scale + Group, FUN=sum, data=rel_inf[rel_inf$quantile == 0.9,])

p = ggplot(rel_inf_agg, aes(x=Scale,y=value,color=Group)) +
        geom_line() +
        geom_point(size=1) +
        facet_wrap(~ temp_var) +
        scale_color_manual(values = muted(c("blue","brown","green")),
                           breaks=c("Elevation","Microtopography","Vegetation"),
                           guide=guide_legend(title=NULL)) +
        scale_x_log10() +
        expand_limits(y=0) +
        annotation_logticks(sides="b") +
        xlab("Scale (m)") +
        ylab("Relative influence (%)") +
        theme_bw() +
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              panel.border=element_rect(color="black"),
              panel.grid.minor.y=element_blank(),
              legend.position="bottom",
              legend.key.width=unit(2,"lines"),
              legend.margin=margin(1,1,1,1))

dir.create("output/quantile")
png("output/quantile/influence_by_scale.png", width=7, height=5, units="in", res=300)
print(p)
dev.off()

rel_inf_agg_scale = aggregate(value ~ temp_var + Variable + Group + quantile, FUN=sum, data=rel_inf)
rel_inf_agg_scale$Variable = factor(rel_inf_agg_scale$Variable, levels=rast_df$Variable[c(1:15,17:21,16)])
rel_inf_agg_scale$Variable_text = as.character(rel_inf_agg_scale$Variable)
rel_inf_agg_scale$Variable_text[rel_inf_agg_scale$value < 2.5] = ""
rel_inf_agg_scale$temp_var = factor(rel_inf_agg_scale$temp_var, levels=unique(rel_inf_agg_scale$temp_var))

# https://stackoverflow.com/questions/34903368/how-to-center-stacked-percent-barchart-labels/34904604#34904604
p = ggplot(rel_inf_agg_scale, aes(x=quantile,y=value,fill=format_names(Variable))) +
        geom_bar(stat="identity") +
        facet_wrap(~ temp_var, nrow=6) +
        scale_fill_manual(values = unique(rast_df$Variable_color)[c(1:15,17:21,16)],
                          guide=guide_legend(title=NULL,nrow=21)) +
        xlab("Quantile") +
        ylab("Relative influence (%)") +
        theme_bw() +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0))) +
        scale_x_continuous(breaks=seq(0.1,0.9,by=0.1), labels=paste0(100*seq(0.1,0.9,by=0.1),"%")) +        
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              panel.border=element_rect(color="black"),
              panel.grid.minor.y=element_blank(),
              legend.position="right",
              legend.key.width=unit(2,"lines"),
              legend.margin=margin(1,1,1,1))

png("output/quantile/influence_by_quantile.png", width=8, height=6.5, units="in", res=300)
print(p)
dev.off()

# 90% quantile only
rel_inf_agg_scale = aggregate(value ~ temp_var + Variable + Group, FUN=sum, data=rel_inf[rel_inf$quantile == 0.9,])
rel_inf_agg_scale$Variable = factor(rel_inf_agg_scale$Variable, levels=rast_df$Variable)
rel_inf_agg_scale$Variable_text = as.character(rel_inf_agg_scale$Variable)
rel_inf_agg_scale$Variable_text = format_names(as.character(rel_inf_agg_scale$Variable))
rel_inf_agg_scale$Variable_text[rel_inf_agg_scale$value < 2.5] = ""
rel_inf_agg_scale$temp_var = factor(rel_inf_agg_scale$temp_var, levels=unique(rel_inf_agg_scale$temp_var))

p = ggplot(rel_inf_agg_scale, aes(x=temp_var,y=value,fill=Variable)) +
        geom_bar(stat="identity", show.legend=FALSE) +
        scale_fill_manual(values = unique(rast_df$Variable_color)) +
        xlab("Response variable") +
        ylab("Relative influence (%)") +
        theme_bw() +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0))) +        
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              panel.border=element_rect(color="black"),
              panel.grid.minor.y=element_blank(),
              legend.position="bottom",
              legend.key.width=unit(2,"lines"),
              legend.margin=margin(1,1,1,1)) +
        geom_text(aes(label=Variable_text),
                  position=position_stack(vjust=0.5), size=3)

png("output/quantile/influence_simple.png", width=8, height=6.5, units="in", res=300)
print(p)
dev.off()

pdf("output/quantile/influence_simple.pdf", width=8, height=6.5)
print(p)
dev.off()

df = aggregate(value ~ Group + temp_var, FUN=sum, data=rel_inf[rel_inf$quantile == 0.9,]) # paper text
df = df[order(df$Group, df$value, decreasing=TRUE),]
write.csv(df, "output/main_models.csv", row.names=FALSE)

