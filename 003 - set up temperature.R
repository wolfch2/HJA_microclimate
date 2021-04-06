# days in each month
# https://rdrr.io/github/tsangyp/StreamThermal/src/R/T_magnitude.r https://en.wikipedia.org/wiki/Growing_degree-day
monthdays<-Vectorize(function(month){
  if(month==2){
    MonthDays<-28 # ignores leap years, I guess...
  }else if(month==1|month==3|month==5|month==7|month==8|month==10|month==12){
    MonthDays<-31
  }else{
    MonthDays<-30
  }
  return(MonthDays)
})

# compute mean/max/min provided at least missing_frac proportion of obs. are present
missing_mean = function(x, n_days, missing_frac=4/5){
        out = mean(x)
        if(length(x) < missing_frac * n_days) out = NA
        return(out)
}
missing_max = function(x, n_days, missing_frac=4/5){
        out = max(x)
        if(length(x) < missing_frac * n_days) out = NA
        return(out)
}
missing_min = function(x, n_days, missing_frac=4/5){
        out = min(x)
        if(length(x) < missing_frac * n_days) out = NA
        return(out)
}

T_custom_metrics = function(sitedata, missing_GDD = 1, missing_other = 4/5){ # see StreamThermal package for input format
        sitedata = separate(sitedata, "Date", c("Year", "Month", "Day"), sep = "-")
        is_leap = ((as.numeric(sitedata$Year)[1] %% 4) == 0) & # divis. by 4 and not "a centemnial not divis by 400" 
                ! ((as.numeric(sitedata$Year)[1] %% 100) == 0 & (as.numeric(sitedata$Year)[1] %% 400) != 0)
        winter_days = sum(monthdays(1:3)) + is_leap
        total_days = 365 + is_leap        
        # calculate raw metrics
        GDD_winter_5 = sum(pmax((sitedata$MaxT + sitedata$MinT)/2 - 5, 0)[sitedata$Month <= 3])
        # change to NA if we exceed missing_frag missing days
        if(sum(sitedata$Month <= 3) < missing_GDD * winter_days) GDD_winter_5=NA
        # rescale to 90 days (winter, no leapyear) or 365 days (annual, no leapyear)
        # GDD_winter_5 = GDD_winter_5 * 90 / sum(sitedata$Month <= 3)
        #
        out = data.frame(SiteInfo=sitedata$SiteInfo[1], # max daily mean looks better than mean daily max (maybe mean daily max is getting thrown off more by errors?)
                         Jul_Sep_mean_max=missing_mean(sitedata$MaxT[as.numeric(sitedata$Month) %in% 7:9],sum(monthdays(7:9)),missing_other),
                         Jul_Sep_mean_mean=missing_mean(sitedata$MeanT[as.numeric(sitedata$Month) %in% 7:9],sum(monthdays(7:9)),missing_other),
                         GDD_winter_5=GDD_winter_5,
                         Apr_Jun_mean_max=missing_mean(sitedata$MaxT[as.numeric(sitedata$Month) %in% 4:6],sum(monthdays(4:6)),missing_other),
                         Apr_Jun_mean_mean=missing_mean(sitedata$MeanT[as.numeric(sitedata$Month) %in% 4:6],sum(monthdays(4:6)),missing_other),                         
                         Apr_Jun_mean_min=missing_mean(sitedata$MinT[as.numeric(sitedata$Month) %in% 4:6],sum(monthdays(4:6)),missing_other),
                         stringsAsFactors=FALSE)
        return(out)
}

############################## set up temperature data

temp_files = list.files("data_input/temperature_cleaned",full.names=TRUE)
temperature = rbindlist(pblapply(temp_files, function(temp_file){
        name = as.numeric(gsub("015","",str_split(temp_file,"_")[[1]][4]))
        temp = fread(temp_file)        
        temp[, "LOCATION_CODE" := name] # https://stackoverflow.com/questions/19072053/adding-columns-to-a-data-table
        return(temp)
}))
date_time_sep = stri_split_fixed(temperature$DATE_TIME, " ", 2, simplify=TRUE)
temperature[, "Date" := date_time_sep[,1]]
temperature[, "Time" := date_time_sep[,2]]

# https://stackoverflow.com/questions/12064202/apply-several-summary-functions-on-several-variables-by-group-in-one-call
temperature = temperature[ , .(AIRTEMP_MIN_DAY = min(TEMP_C),
                                AIRTEMP_MEAN_DAY = mean(TEMP_C),
                                AIRTEMP_MAX_DAY = max(TEMP_C)), by = .(LOCATION_CODE, Date)] # aggregate to day level

temperature = separate(temperature, "Date", c("Year","Month","Day"), "-", convert=TRUE, remove=FALSE)
temperature$Date = ymd(temperature$Date)

############################## join with gridMET

GRIDMET = read.csv("data_input/GRIDMET_reduce.csv", as.is = TRUE)
GRIDMET$Date = mdy(GRIDMET$date)
GRIDMET$min = GRIDMET$tmmn - 273.15
GRIDMET$max = GRIDMET$tmmx - 273.15

temperature_merged = merge(temperature, GRIDMET[,c("Date","min","max")], all.x=TRUE, all.y=FALSE, by="Date")

############################## join with VANMET

weather = read.csv("data_input/MS00101_v9.csv", as.is=TRUE)
weather = weather[weather$HEIGHT == 450 & weather$SITECODE == "VANMET",]
weather = data.table(weather)[ , .(min_w = mean(AIRTEMP_MIN_DAY, na.rm=TRUE),
                                   mean_w = mean(AIRTEMP_MEAN_DAY, na.rm=TRUE),
                                   max_w = mean(AIRTEMP_MAX_DAY, na.rm=TRUE)), by = .(SITECODE, DATE)]
weather = na.omit(weather)
weather$Date = ymd(weather$DATE)
#range(year(weather$DATE)) # 1987-2016
#weather = weather[year(weather$DATE) >= 2009,] # optional -- match with our analysis time span

temperature_merged_again = merge(temperature_merged, weather[,c("Date","min_w","mean_w","max_w")], all.x=TRUE, all.y=FALSE, by="Date")

############################## compute annual metrics (three ways: unadjusted, relative to gridMET free-air, and relative to VANMET free-air)

temp_split = split(temperature_merged_again,paste(temperature_merged_again$Year,temperature_merged_again$LOCATION_CODE))
temperature_metrics = foreach(i=1:length(temp_split)) %dopar% {
        print(i)
        df = temp_split[[i]]
        ### 1. compute unadjusted metrics
        sitedata_micro = data.frame(SiteInfo=df$LOCATION_CODE,
                      Date=df$Date,
                      MaxT=df$AIRTEMP_MAX_DAY,
                      MinT=df$AIRTEMP_MIN_DAY,
                      MeanT=df$AIRTEMP_MEAN_DAY,
                      stringsAsFactors=FALSE)
        temps = t(T_custom_metrics(sitedata_micro)[,-1])
        out = data.frame(site=df$LOCATION_CODE[1],
                   Year=df$Year[1],
                   variable=rownames(temps),
                   value=temps,
                   stringsAsFactors=FALSE)
        ### 2. compute delta metrics (difference between fine-scale and gridMET free-air)
        sitedata_GRIDMET = data.frame(SiteInfo=df$LOCATION_CODE,
                      Date=df$Date,
                      MaxT=df$max, # just use GRIDMET data here...
                      MinT=df$min,
                      MeanT=(df$max + df$min) / 2,
                      stringsAsFactors=FALSE)
        out$GRIDMET_value = t(T_custom_metrics(sitedata_GRIDMET)[,-1])
        out$delta_metrics = out$value - out$GRIDMET_value
        ### 3. compute VANMET delta metrics (difference between fine-scale and VANMET free-air)
        sitedata_VANMET = data.frame(SiteInfo=df$LOCATION_CODE,
                      Date=df$Date,
                      MaxT=df$max_w, # just use GRIDMET data here...
                      MinT=df$min_w,
                      MeanT=df$mean_w,
                      stringsAsFactors=FALSE) %>% na.omit # VANMET has a few missing obs. here and there (which appear as NA)
        out$VANMET_value = t(T_custom_metrics(sitedata_VANMET)[,-1])
        out$delta_metrics_VANMET = out$value - out$VANMET_value
        return(out)
}
temperature_metrics = data.frame(rbindlist(temperature_metrics))
temperature_metrics = temperature_metrics[! is.na(temperature_metrics$value),]
table(temperature_metrics$variable, temperature_metrics$Year) # looks fine..

temperature_metrics$POINT = temperature_metrics$LOCATION_CODE = temperature_metrics$SITECODE = temperature_metrics$site

saveRDS(temperature_metrics, "data_processed/temperature_metrics.RDS")
temperature_metrics = readRDS("data_processed/temperature_metrics.RDS")

############################## map # days

boundaries = read_sf("data_input/boundaries/boundary83.shp")

sites = readRDS("data_processed/sites.RDS")
temperature_sites = merge(temperature, sites, by="LOCATION_CODE")
temp_days = data.table(temperature_sites)[ , .(count=.N), by = .(POINT, Year, X, Y)]

p = ggplot(boundaries) +
        geom_sf(fill="lightgray", color="black", size=0.3) +
        geom_point(data=temp_days, aes(x=X,y=Y,fill=count),size=0.75,pch=21,color="black") +
        facet_wrap(~ Year, nrow=4) +
        scale_size_manual(values=c(0.25,0.25,0.25,1.75), guide = guide_legend(title = "Site type")) +        
        scale_fill_gradientn(colours=brewer.pal(11,"Spectral"),
                                      limits=c(0,366),
                                      guide = guide_colorbar(title = "Days")) +        
        theme_bw() +
        theme(legend.direction="horizontal",
                      legend.position=1:0,
                      legend.justification=1:0,
                      axis.ticks=element_blank(),
                      axis.title=element_blank(),
                      axis.text=element_blank(),
                      panel.grid.major = element_line(colour="transparent"),
                      panel.grid.minor = element_line(colour="transparent"),
                      legend.background=element_rect(fill=NA),
                      legend.box="horizontal")

png("output/days_map.png", width=6, height=7, units="in", res=300)
print(p)
dev.off()

############################## plot example time series

site_obs_count = temperature_metrics %>%
        group_by(site) %>%
        summarize(count = n()) %>%
        arrange(desc(count))

example_sites = data.frame(read_excel("data_input/site_locations.xlsx")) %>%
        arrange(Plantation) %>%
        merge(site_obs_count, by.x="POINT", by.y="site") %>%
        arrange(desc(count)) %>%
        group_by(Plantation) %>%
        slice_head(n=4) %>%
        inner_join(temperature_metrics) %>%
        mutate(var = format_names(variable),
               Plantation = factor(Plantation,
                                   levels=c(1,0),
                                   labels=c("Plantation","Mature forest/old growth"))) %>%
        arrange(Plantation) %>%
        mutate(POINT = factor(POINT,levels=unique(.$POINT)))

p = ggplot(example_sites, aes(x=Year,y=delta_metrics,color=Plantation,group=POINT)) +
        geom_point() +
        geom_line() +
        facet_grid(var ~ POINT, scales="free_y") +
                scale_color_manual(values=brewer.pal(3,"Set1")[3:2],
                           guide=guide_legend(title=NULL, order=1)) +
        theme_bw() +
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              panel.border=element_rect(color="black"),
              legend.position="bottom",
              panel.spacing.x=unit(1.5,"lines"),
              legend.background=element_rect(color="black")) +
        ylab("Value relative to free-air") +
        geom_hline(yintercept=0, linetype="dashed") +
        scale_x_continuous(breaks=seq(2009,2018,by=3)) +
        xlab("Year")

png("output/example_ts.png", width=12, height=9, units="in", res=300)
print(p)
dev.off()

############################## compare gridMET and VANMET for our six temperature metrics

df = temperature_metrics %>%
        filter(! duplicated(paste(Year,variable))) %>% # ignore duplicates (same across sites)
        na.omit %>%
        mutate(variable=format_names(variable))

# https://stats.stackexchange.com/questions/228540/how-to-calculate-out-of-sample-r-squared
stats = do.call("rbind", lapply(split(df, df$variable), function(df_var){
        r = cor(df_var$GRIDMET_value, df_var$VANMET_value) %>% formatC(format='f', digits=2)
        b = lm(VANMET_value ~ GRIDMET_value, data=df_var)$coef[2] %>% formatC(format='f', digits=2)
        out = data.frame(variable=df_var$variable[1],
                         Year=df_var$Year[1],
                         lab=paste0("' '*r*': '*", r, "*', '*beta[1]*': '*", b))
        return(out)
}))
df$year = factor(df$Year); stats$year = factor(stats$Year);

dummy = df; dummy$GRIDMET_value = df$VANMET_value; dummy$VANMET_value = df$GRIDMET_value;

###

p = ggplot(df, aes(x=GRIDMET_value, y=VANMET_value)) +
        geom_smooth(method="lm") +
        geom_point(aes(color=year)) +
        facet_wrap(~ variable, scales="free") +
        geom_point(data=dummy, color="#FFFFFF00") +
        geom_abline(intercept=0,slope=1,linetype="dashed") +
        scale_x_continuous(breaks= pretty_breaks(4)) +
        scale_y_continuous(breaks= pretty_breaks(4)) +
        scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(df$year))),
                           guide=guide_legend(title=NULL, nrow=1,
                                              override.aes = list(size = 1.5))) +
        geom_text(data=stats, aes(label=lab),
                  x=Inf, y=-Inf, color="black", vjust=-1.15, hjust=1.05, size=2.5, parse=TRUE) +
        theme_bw() +
        theme(axis.text=element_text(color="black"),
              axis.ticks=element_line(color="black"),
              axis.title=element_blank(),
              panel.border=element_rect(color="black"),
              legend.background=element_rect(color="black"),
              legend.position="bottom",
              aspect.ratio = 1,
              plot.title = element_text(hjust = 0.5)) 

png("output/gridMET_VANMET_metrics.png", width=8.5, height=6, units="in", res=300)
print(p)
dev.off()
