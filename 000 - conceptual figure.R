set.seed(2)
x1 = rnorm(500,mean=23,sd=0.85)
x2 = rnorm(500,mean=27,sd=0.15)

df = data.frame(x=c(x1,x2),gp=rep(1:2,each=length(x1)))
df$gp_name = factor(df$gp, levels=2:1,
                    labels=c("Stable microclimate","Stable microrefugium"))

p = ggplot(df, aes(x=x,group=gp_name,fill=gp_name)) +
  geom_histogram(position="identity",binwidth=0.2,boundary=0,color="black") +
  theme_classic() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position=c(0.05,0.95),
        legend.justification=c(0,1)) +
  geom_vline(xintercept=26, lwd=1, color="red") +
  xlab("Maximum summer temperature") +
  ylab("Frequency") +
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  scale_x_continuous(limits=c(20,28)) +
  scale_fill_manual(values=brewer.pal(11,"Paired")[c(1,4)],
                    guide=guide_legend(title=NULL))

pdf("output/conceptual_1.pdf", width=4.5, height=4.5) # for paper
print(p)
dev.off()

png("output/conceptual_1.png", width=4.5, height=4.5, res=450, units="in")
print(p)
dev.off()

df = expand.grid(mean=seq(1,2,length=100),sd=seq(1,2,length=100))
df$quan = qnorm(0.9,mean=df$mean,sd=df$sd)

p = ggplot(df, aes(x=mean,y=sd,fill=quan)) +
    geom_raster() +
    scale_fill_gradientn(colors=rev(brewer.pal(11,"Spectral")),
                         guide=guide_colorbar(title="90% quantile")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    xlab("Mean") +
    ylab("Standard deviation") +
    theme(axis.ticks=element_line(color="black"),
          axis.text=element_text(color="black"),
          panel.border=element_rect(color="black",fill=NA)) +
      coord_fixed()

png("output/conceptual_2.png", width=1.2*4.5, height=1.2*3.5, res=450, units="in")
print(p)
dev.off()
