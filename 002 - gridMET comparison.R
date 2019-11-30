# https://stackoverflow.com/questions/743812/calculating-moving-average
ma <- function(x, n = 5){as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))}

################################ gridMET

temp = read.csv("data_input/GRIDMET_reduce.csv", as.is = TRUE)
temp$tmmx = temp$tmmx - 273.15
temp$tmmn = temp$tmmn - 273.15
temp$Date = mdy(temp$date)
temp$yday = yday(temp$Date)
temp$Year = year(temp$Date)
temp$max = as.numeric(ma(temp$tmmx, 30))
temp$min = as.numeric(ma(temp$tmmn, 30))

temp_long = melt(temp, measure.vars=c("max","min"))
temp_long$gp = paste(temp_long$Year, temp_long$variable)

p = ggplot(temp_long, aes(x=yday,y=value,color=Year,group=gp)) +
  geom_line() +
  scale_color_gradientn(colors=rev(brewer.pal(11,"Spectral"))) +
  theme_bw() +
  theme(axis.text=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        panel.border=element_rect(color="black"),
        legend.position=c(0.05,0.95),
        legend.justification=0:1,
        legend.background=element_rect(color="black")) +
  ylab("Temperature (°C)") +
  xlab("Day of year")

################################ join w/ VANMET

weather = read.csv("data_input/MS00101_v8.csv", as.is=TRUE)
table(weather$SITECODE, weather$HEIGHT)
weather = weather[weather$HEIGHT == 450 & weather$SITECODE == "VANMET",]
weather = data.table(weather)[ , .(AIRTEMP_MIN_DAY = mean(AIRTEMP_MIN_DAY, na.rm=TRUE),
                                AIRTEMP_MEAN_DAY = mean(AIRTEMP_MEAN_DAY, na.rm=TRUE),
                                AIRTEMP_MAX_DAY = mean(AIRTEMP_MAX_DAY, na.rm=TRUE)), by = .(SITECODE, DATE)]
weather = na.omit(weather)
weather$Date = ymd(weather$DATE)
weather$max_w = as.numeric(weather$AIRTEMP_MAX_DAY)
weather$min_w = as.numeric(weather$AIRTEMP_MIN_DAY)
range(year(weather$DATE)) # 1987-2016
weather = weather[year(weather$DATE) >= 2009,] # optional -- match with our analysis time span

temp_long = melt(temp, measure.vars=c("max","min"))
temp_long$gp = paste(temp_long$Year, temp_long$variable)

df = merge(temp, weather, by="Date")

max_cor = round(cor(na.omit(cbind(df$tmmx, df$max_w)))[1,2],3)
min_cor = round(cor(na.omit(cbind(df$tmmn, df$min_w)))[1,2],3)

df_max = df
df_max$GRIDMET = df_max$tmmx
df_max$VANMET = df_max$max_w
df_min = df
df_min$GRIDMET = df_min$tmmn
df_min$VANMET = df_min$min_w
df_long = rbind(data.frame(df_max,cat=paste0("Maximum (r = ", max_cor, ")")),
                data.frame(df_min,cat=paste0("Minimum (r = ", min_cor, ")"))) # 2009-2016

q = ggplot(df_long, aes(x=GRIDMET,y=VANMET)) +
  geom_point(alpha=0.1) +
  theme_bw() +
  theme(axis.text=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        panel.border=element_rect(color="black"),
        legend.position=c(0.05,0.95),
        legend.justification=0:1,
        legend.background=element_rect(color="black")) +
  facet_wrap(~ cat) +
  ylab("VANMET temperature (°C)") +
  xlab("gridMET temperature (°C)") +
  geom_abline(intercept=0,slope=1,color="red",lwd=1,linetype="dashed")

dir.create("output/GRIDMET")
png("output/GRIDMET/daily_GRIDMET_VANMET.png", width=9, height=9, units="in", res=400)
print(plot_grid(p, q, labels = c('A', 'B'), label_size = 12, nrow=2))
dev.off()

