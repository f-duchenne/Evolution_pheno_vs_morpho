###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","viridis") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
empf=fread("empirical_networks.csv")
names(empf)=gsub("weighted ","",names(empf))
names(empf)[3]="inter. evenness"
fwrite(unique(empf[,c("site","na","np","nround","nyear")]),"table1.csv")
empf2=subset(empf,competition==5 & rho==0.05)
colo=c("dodgerblue4","chartreuse3")

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info_empir.csv")
names(indf)=gsub("weighted ","",names(indf))
names(indf)[3]="inter. evenness"

indf2=subset(indf,rho==0.01)
acp=PCA(indf2[,1:6],graph=FALSE)
indf2=cbind(indf2,acp$ind$coord)

empf=cbind(empf,predict(acp,newdata=empf[,1:6])$coord)

plot_acp=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf2,time==2000),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait))+
geom_point(data=subset(indf2,time==0),aes(x=Dim.1,y=Dim.2,alpha=time),color="grey",alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("a")+geom_point(data=empf,aes(x=Dim.1,y=Dim.2),color="black")

pl2=fviz_pca_var(acp, col.var = "black", repel =TRUE)+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold",hjust = 0))

png("Fig.S5.png",width=1200,height=800,res=160)
plot_grid(plot_acp,pl2,ncol=2,align="hv")
dev.off();


pl1=ggplot()+
geom_boxplot(data=empf2,aes(x="Empirical",y=motifn/motifntot))+
geom_boxplot(data=subset(indf,time==2000 & rho==0.05),aes(y=motifn/motifntot,x="Simulations",color=trait))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
axis.title.x=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+ggtitle("a")+
ylab("Percentage of 'V' motifs promoting facilitation")+scale_y_continuous(labels=scales::percent)


pl2=ggplot(data=empf2,aes(x=feas_without_pheno,y=feas_with_pheno))+
geom_point(aes(color=na+np),size=2)+
geom_abline(intercept=0,slope=1)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position = c(0.75, 0.3))+
scale_color_viridis(option="E")+ggtitle("b")+coord_fixed(ratio=1,expand=F)+ylim(c(0,0.5))+xlim(c(0,0.5))+
ylab("Feasibility when including seasonal structure")+
xlab("Feasibility when neglecting seasonal structure")+labs(color=expression(N[sp]))


pl3=ggplot(data=empf,aes(x=as.factor(competition),y=feas_without_pheno/feas_with_pheno-1,color=as.factor(rho)))+
geom_hline(linetype="dashed",yintercept=0)+
geom_hline(linetype="dashed",yintercept=-1)+
geom_boxplot()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+ggtitle("c")+
scale_y_continuous(labels=scales::percent)+ylab("Decrease in feasibility\nwhen neglecting seasonal structure")+
xlab("Competition strength")+scale_colour_viridis(option="D",discrete=TRUE)+labs(color=expression(rho))


plot_grid(pl1,pl2,pl3,ncol=3,align="hv",rel_widths=c(0.75,2,1.5))


pdf("Fig.5.pdf",width=10,height=5)
plot_grid(pl1,pl2,pl3,ncol=3,align="hv",rel_widths=c(1,2.5,1.5))
dev.off();



