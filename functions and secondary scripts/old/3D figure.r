library(tidyverse)
library(mvtnorm)
library(plotly)
library(MASS)
library(cowplot)
library(gridExtra)

sigma1 <- matrix(c(3,0,0,3), ncol = 2)
sigma1
means1 <- c(1, 5)

sigma2 <- matrix(c(6,0,0,5), ncol = 2)
sigma2
means2 <- c(2, -2)

n <- 1000
set.seed(42)
x1 <- rmvnorm(n = n, mean = means1, sigma = sigma1)
x2 <- rmvnorm(n = n, mean = means2, sigma = sigma2)
d <- rbind(data.frame(x1),data.frame(x2))

p1 <- ggplot(d, aes(x = X1, y = X2)) +
  stat_density_2d(aes(fill = ..density..),geom = "raster", contour = FALSE) +
  scale_fill_viridis_c()+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0.5),axis.text=element_blank(),axis.ticks=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+coord_fixed(ratio=1,expand=F)+labs(fill="fitness")+xlab("Trait value species 1")+ylab("Trait value species 2")+xlim(c(-5,8))+ylim(c(-5,8))+
geom_point(shape=21,x=2.5,y=0.75,fill="black",col="white",size=4)+geom_point(shape=21,x=1,y=5,fill="orange",col="white",size=4)+geom_point(shape=21,x=1.5,y=-3,fill="dodgerblue",col="white",size=4)+ggtitle("Fitness landscape of species 1")


sigma1 <- matrix(c(2,0,0,3), ncol = 2)
sigma1
means1 <- c(-1, 4)

sigma2 <- matrix(c(1,0,0,1.2), ncol = 2)
sigma2
means2 <- c(1.5, -3)

sigma3 <- matrix(c(3,0,0,2), ncol = 2)
sigma3
means3 <- c(2, 5)

n <- 1000
set.seed(42)
x1 <- rmvnorm(n = n, mean = means1, sigma = sigma1)
x2 <- rmvnorm(n = n, mean = means2, sigma = sigma2)
x3 <- rmvnorm(n = n, mean = means3, sigma = sigma3)
d <- rbind(data.frame(x1),data.frame(x2),data.frame(x3))

p2 <- ggplot(d, aes(x = X1, y = X2)) +
  stat_density_2d(aes(fill = ..density..),geom = "raster", contour = FALSE) +
  scale_fill_viridis_c()+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0.5),legend.position="none",axis.text=element_blank(),axis.ticks=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+coord_fixed(ratio=1,expand=F)+labs(fill="fitness")+xlab("Trait value species 1")+ylab("Trait value species 2")+xlim(c(-5,8))+ylim(c(-5,8))+
geom_point(shape=21,x=2.5,y=0.75,fill="black",col="white",size=4)+geom_point(shape=21,x=1,y=5,fill="orange",col="white",size=4)+geom_point(shape=21,x=1.5,y=-3,fill="dodgerblue",col="white",size=4)+ggtitle("Fitness landscape of species 2")



leg <- ggpubr::as_ggplot(cowplot::get_legend(p1))
p1=p1+theme(legend.position="none")

top=plot_grid(p1,p2,ncol=2,align="hv",rel_widths=c(1,1,0.2))

top

sigma1 <- matrix(c(6,0,0,3), ncol = 2)
sigma1
means1 <- c(-3, 3)

sigma2 <- matrix(c(5,0,0,2), ncol = 2)
sigma2
means2 <- c(3.5, -2)

sigma3 <- matrix(c(6,0,0,8), ncol = 2)
sigma3
means3 <- c(2, 3)

n <- 1000
set.seed(42)
x1 <- rmvnorm(n = n, mean = means1, sigma = sigma1)
x2 <- rmvnorm(n = n, mean = means2, sigma = sigma2)
x3 <- rmvnorm(n = n, mean = means3, sigma = sigma3)
d <- rbind(data.frame(x1),data.frame(x2),data.frame(x3)) 
 
p3 <- ggplot(d, aes(x = X1, y = X2)) +
  stat_density_2d(aes(fill = ..density..),geom = "raster", contour = FALSE) +
  scale_fill_viridis_c(option = "A")+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0.5),legend.position="none",axis.text=element_blank(),axis.ticks=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+coord_fixed(ratio=1,expand=F)+labs(fill="fitness")+xlab("Trait value species 1")+ylab("Trait value species 2")+xlim(c(-5,8))+ylim(c(-5,8))+
geom_point(shape=21,x=1,y=5,fill="orange",col="white",size=4)+geom_point(shape=21,x=1.5,y=-3,fill="dodgerblue",col="white",size=4)+ggtitle("Network structure (e.g. connectance)")


x1=data.frame(x1)
x2=data.frame(x2)
x3=data.frame(x3)
mat1=kde2d(x1[,1],x1[,2],lims = c(c(-5,8), c(-5,8)),n =100)
mat2=kde2d(x2[,1],x2[,2],lims = c(c(-5,8), c(-5,8)),n =100)
mat3=kde2d(x3[,1],x3[,2],lims = c(c(-5,8), c(-5,8)),n =100)

mat=mat1
mat$z[,]="1"
mat$z[mat1$z>mat2$z]="2"
mat$z[mat1$z>mat3$z]="3"
mat$z[mat1$z>mat3$z & mat3$z>mat2$z]="4"
mat$z[mat1$z<mat3$z & mat3$z>mat2$z]="5"

d=rcosmosis::vmat2df(mat)


p4 <- ggplot(data=d, aes(x = x, y = y,fill=z,group=z)) +
  geom_raster() +
  scale_fill_viridis_c(option = "E")+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0.5),legend.position="none",axis.text=element_blank(),axis.ticks=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+coord_fixed(ratio=1,expand=F)+labs(fill="fitness")+xlab("Trait value species 1")+ylab("Trait value species 2")+xlim(c(-5,8))+ylim(c(-5,8))+
geom_point(shape=21,x=1,y=5,fill="orange",col="white",size=4)+geom_point(shape=21,x=1.5,y=-3,fill="dodgerblue",col="white",size=4)+ggtitle("Structural stability")


grid.arrange(top,p3,p4,ncol=1)
 
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
pdf("Fig3D.pdf",width=6,height=10)
grid.arrange(top,p3,p4,ncol=1)
dev.off();