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

T_custom_metrics = function(sitedata, missing_frac = 4/5){ # see StreamThermal package for input format
        sitedata = separate(sitedata, "Date", c("Year", "Month", "Day"), sep = "-")
        is_leap = ((as.numeric(sitedata$Year)[1] %% 4) == 0) & # divis. by 4 and not "a centemnial not divis by 400" 
                ! ((as.numeric(sitedata$Year)[1] %% 100) == 0 & (as.numeric(sitedata$Year)[1] %% 400) != 0)
        winter_days = sum(monthdays(1:3)) + is_leap
        total_days = 365 + is_leap        
        # calculate raw metrics
        GDD_winter_5 = sum(pmax((sitedata$MaxT + sitedata$MinT)/2 - 5, 0)[sitedata$Month <= 3])
        # change to NA if we exceed missing_frag missing days
        if(sum(sitedata$Month <= 3) < missing_frac * winter_days) GDD_winter_5=NA
        # rescale to 90 days (winter, no leapyear) or 365 days (annual, no leapyear)
        GDD_winter_5 = GDD_winter_5 * 90 / sum(sitedata$Month <= 3)
        #
        out = data.frame(SiteInfo=sitedata$SiteInfo[1], # max daily mean looks better than mean daily max (maybe mean daily max is getting thrown off more by errors?)
                         Jul_Sep_mean_max=missing_mean(sitedata$MaxT[as.numeric(sitedata$Month) %in% 7:9],sum(monthdays(7:9)),missing_frac),
                         Jul_Sep_mean_mean=missing_mean(sitedata$MeanT[as.numeric(sitedata$Month) %in% 7:9],sum(monthdays(7:9)),missing_frac),
                         GDD_winter_5=GDD_winter_5,
                         Apr_Jun_mean_max=missing_mean(sitedata$MaxT[as.numeric(sitedata$Month) %in% 4:6],sum(monthdays(4:6)),missing_frac),
                         Apr_Jun_mean_mean=missing_mean(sitedata$MeanT[as.numeric(sitedata$Month) %in% 4:6],sum(monthdays(4:6)),missing_frac),                         
                         Apr_Jun_mean_min=missing_mean(sitedata$MinT[as.numeric(sitedata$Month) %in% 4:6],sum(monthdays(4:6)),missing_frac),
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

############################## compute annual metrics (two ways: unadjusted and relative to free-air)

temp_split = split(temperature_merged,paste(temperature_merged$Year,temperature_merged$SITECODE))
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
        out = data.frame(site=df$SITECODE[1],
                   Year=df$Year[1],
                   variable=rownames(temps),
                   value=temps,
                   stringsAsFactors=FALSE)
        ### 2. compute delta metrics (difference between fine-scale and free-air)
        sitedata_GRIDMET = data.frame(SiteInfo=df$LOCATION_CODE,
                      Date=df$Date,
                      MaxT=df$max, # just use GRIDMET data here...
                      MinT=df$min,
                      MeanT=(df$max + df$min) / 2,
                      stringsAsFactors=FALSE)
        out$GRIDMET_value = t(T_custom_metrics(sitedata_GRIDMET)[,-1])       
        out$delta_metrics = out$value - out$GRIDMET_value
        return(out)
}
temperature_metrics = data.frame(rbindlist(temperature_metrics))
temperature_metrics = temperature_metrics[! is.na(temperature_metrics$value),]
table(temperature_metrics$variable, temperature_metrics$Year) # looks fine..

temperature_metrics$POINT = temperature_metrics$LOCATION_CODE = temperature_metrics$SITECODE = temperature_metrics$site

saveRDS(temperature_metrics, "data_processed/temperature_metrics.RDS")

############################## map # days

boundaries = read_sf("data_processed/boundaries/boundary83.shp")

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
p
dev.off()

