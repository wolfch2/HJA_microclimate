x = seq(20,28,length=1e3)
y1 = dnorm(x,mean=23,sd=0.85)
y2 = dnorm(x,mean=27,sd=0.15)

df = data.frame(x=c(x,x),y=c(y1,y2),gp=rep(1:2,each=length(x)))
df$gp_name = factor(df$gp, levels=2:1,
                    labels=c("Stable microclimate","Stable microrefugium"))

p_1 = ggplot(df, aes(x=x,y=y,group=gp_name,linetype=gp_name)) +
  geom_line() +
  theme_classic() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position=c(0.05,0.95),
        legend.justification=c(0,1)) +
  geom_vline(xintercept=26, lwd=1, color="red") +
  xlab("Annual statistic") +
  ylab("Density") +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  scale_linetype(guide=guide_legend(title=NULL))

pdf("output/conceptual_1.pdf", width=4.5, height=4.5) # for paper
print(p_1)
dev.off()

png("output/conceptual_1.png", width=4.5, height=4.5, res=450, units="in")
print(p_1)
dev.off()

