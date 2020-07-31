############################## set up data

rast_reduce = readRDS("data_processed/rast_reduce.RDS")
rast_reduce = rast_reduce[! rast_reduce$variable %in% c("PC1","PC2"),] # not using PC vars for main models
rast_reduce$var_scale = paste0(rast_reduce$variable, "XX", rast_reduce$scale)
rast_spread = spread(rast_reduce[,c("site","var_scale","value"),],key="var_scale",value="value")

temp_merged = readRDS("data_processed/temperature_metrics.RDS")
temp_merged = aggregate(cbind(value,delta_metrics) ~ site+variable,FUN=sd,temp_merged)

data = merge(temp_merged, rast_spread, by="site")

rast_df = data.frame(read_excel("data_input/HJA variables_final.xlsx", na="NA"))
rast_df = rast_df[! rast_df$Variable %in% c("PC1","PC2"),]
rast_df$Variable_color = c(colorRampPalette(brewer.pal(9,"Greens"))(table(rast_df$Group)["Vegetation"]),
                           brewer.pal(9,"Blues")[9],
                           brewer.pal(6,"YlOrBr")[(6-table(rast_df$Group)["Microtopography"]+1):6]) %>%
                           alpha(0.8)

############################## regression models for imp. plots etc.

pred_mat = expand.grid(var=unique(data$variable), stringsAsFactors=FALSE)
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
                alpha=0.5,
                nrounds = 25,
                nthread = 8)
        return(mod)
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

rel_inf_agg_scale = aggregate(value ~ temp_var + Variable + Group, FUN=sum, data=rel_inf)
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
              panel.grid.major.x=element_blank(),
              legend.position="bottom",
              legend.key.width=unit(2,"lines"),
              legend.margin=margin(1,1,1,1)) +
        geom_text(aes(label=Variable_text),
                  position=position_stack(vjust=0.5), size=3)

png("output/quantile/influence_sd.png", width=8, height=6.5, units="in", res=300)
print(p)
dev.off()

df = aggregate(value ~ Group + temp_var, FUN=sum, data=rel_inf) # paper text
df = df[order(df$Group, df$value, decreasing=TRUE),]
write.csv(df, "output/sd_models.csv", row.names=FALSE)

