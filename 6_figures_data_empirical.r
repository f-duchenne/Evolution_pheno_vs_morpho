###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","viridis","glmmTMB","lme4") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

path_folder="C:/Users/Duchenne/Documents/evolution_pheno_morpho/data_zenodo/"
setwd(dir=path_folder)

#STRUCTURE EMPIRICAL:
empf=fread("data/empirical/stability_of_empirical_networks.csv")
names(empf)[1:2]=paste0("w",names(empf)[1:2])
names(empf)[3]="inter. evenness"
names(empf)[5]="wNODF"
setcolorder(empf,c(names(empf)[1:5],"modularity","mut.strength",names(empf)[8:ncol(empf)]))
colo=c("chartreuse3","dodgerblue4","gold3")

empf$dive=empf$na+empf$np 
empf$prop_motif=(empf$motifn/empf$motifntot)
empf$comp=apply(empf[,c("comp_a","comp_p")],1,mean)
#model=lmer(feas_with_pheno~rho*(prop_motif+dive)+(1|site),data=subset(empf,competition==5 & rho %in% c(0.01,0.05,0.2)))

#STRUCTURE SIMULATIONS:
indf=fread("data/empirical/outputs_simulations/networks_info_empir.csv")
names(indf)=gsub("weighted ","w",names(indf))
names(indf)[3]="inter. evenness"
names(indf)[5]="wNODF"
names(indf)[7]="mut.strength"
names(indf)[11:12]=c("na","np")
indf$dive=indf$nbsp_a+indf$nbsp_p 
indf$prop_motif=(indf$motifn/indf$motifntot)


bidon=subset(indf,time==2000)
mean(bidon$comp/bidon$mut.strength)

### PROPORTION OF V+ MOTIFS
pl1=ggplot()+
geom_boxplot(data=subset(empf,competition==4),aes(x="Empirical",y=motifn/motifntot))+
geom_boxplot(data=subset(indf,time==2000 & trait=="both" & competition==4),aes(y=motifn/motifntot,x="Simulated\nCoevolved"),color=colo[1])+
geom_boxplot(data=subset(indf,time==0 & trait=="both" & competition==4),aes(y=motifn/motifntot,x="Simulated\nInitial"),color="grey")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
axis.title.x=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
ggtitle("a")+
ylab('Percentage of "V+" motifs')+scale_y_continuous(labels=scales::percent)

### STRUCTURAL STABILITY
empf$dive=empf$na+empf$np
b=melt(empf,id.vars=c("dive","competition"),measure.vars=c("feas_with_pheno","feas_without_pheno"))
b$variable=as.character(b$variable)
b$variable[b$variable=="feas_without_pheno"]="without seasonal structure"
b$variable[b$variable=="feas_with_pheno"]="with seasonal structure"

pl2=ggplot(data=subset(b,competition==4),aes(x=dive,y=value,color=variable,fill=variable))+
geom_point(size=2,alpha=0.9)+
stat_smooth(method="glm",alpha=0.2,,method.args = list(family = quasibinomial(link = 'logit')))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position = c(0.25, 0.15))+scale_fill_manual(values=c("black","grey"))+
scale_color_manual(values=c("black","grey"))+ggtitle("b")+coord_cartesian()+
ylab(expression(paste("Structural stability  ",(omega))))+
xlab("Diversity")+labs(color="",fill="")


pl3=ggplot(data=empf,aes(x=as.factor(competition),y=feas_without_pheno/feas_with_pheno-1))+
geom_hline(linetype="dashed",yintercept=0)+
geom_boxplot()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+ggtitle("c")+
scale_y_continuous(labels=scales::percent)+ylab("Decrease in structural stability\nwhen neglecting seasonal structure")+
xlab("Competition importance")


plot_grid(pl1,pl2,pl3,ncol=3,align="hv",rel_widths=c(1,2,1))


pdf("Figures/Fig.5.pdf",width=10,height=5)
plot_grid(pl1,pl2,pl3,ncol=3,align="hv",rel_widths=c(0.75,2,0.75))
dev.off();


empf2$dive=empf2$na+empf2$np
b=melt(empf2,id.vars="dive",measure.vars=c("feas_with_pheno","feas_without_pheno"))
b$variable=as.character(b$variable)
b$variable[b$variable=="feas_without_pheno"]="without seasonal structure"
b$variable[b$variable=="feas_with_pheno"]="with seasonal structure"

ggplot(data=b,aes(x=dive,y=value^(1/dive),color=variable,fill=variable))+
geom_point(size=2,alpha=0.9)+
stat_smooth(method="glm",alpha=0.2,,method.args = list(family = quasibinomial(link = 'logit')))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position = c(0.75, 0.3))+scale_fill_manual(values=c("black","grey"))+
scale_color_manual(values=c("black","grey"))+ggtitle("b")+coord_cartesian()+
ylab("Size of the feasibility domain")+
xlab("Diversity")+labs(color="",fill="")






##### ACP
acp=PCA(rbind(indf[,c(1:6)],empf[,1:6]),graph=FALSE)

indf=cbind(indf,acp$ind$coord[1:nrow(indf),])
indf$trait2=indf$trait
indf$trait2[indf$time==0]="initial"

empf=cbind(empf,acp$ind$coord[(nrow(indf)+1):(nrow(indf)+nrow(empf)),])
empf$trait2="empirical"

######################################FIGURE

plot_acp=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf,time==2000),aes(x=Dim.1,y=Dim.2,color=trait2))+
geom_point(data=subset(indf,time==0),aes(x=Dim.1,y=Dim.2,color=trait2),alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=c(colo[1],"black","grey",colo[2:3]))+ggtitle("a")+geom_point(data=empf,aes(x=Dim.1,y=Dim.2,color=trait2))+
labs(color="")

pl2=fviz_pca_var(acp, col.var = "black", repel =TRUE)+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold",hjust = 0))

dists=NULL
for(i in unique(indf$site)){
for(tr in c("both","pheno")){
bidon1=subset(indf,site==i & trait==tr & time==2000 & competition==4)
bidon2=subset(empf,site==i & competition==4)

dists=rbind(dists,data.frame(site=i, distance_to_simu=mean(sqrt((bidon1$Dim.1-bidon2$Dim.1)^2+(bidon1$Dim.2-bidon2$Dim.2)^2+((bidon1$Dim.3-bidon2$Dim.3)^2))),trait=tr))
}}

dist_plot=ggplot(data=dists,aes(x=trait,y=distance_to_simu,color=trait))+
geom_boxplot()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
axis.title.x=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
xlab("Simulation with")+
ylab("Structural distance to empirical networks")+
scale_color_manual(values=colo)+ggtitle("c")+
labs(color="")

nodf=ggplot()+
geom_boxplot(data=empf2,aes(x="Empirical",y=wNODF))+
geom_boxplot(data=subset(indf2,time==2000),aes(y=wNODF,x="Simulated\nCoevolved",color=trait))+
geom_boxplot(data=subset(indf2,time==0),aes(y=wNODF,x="Simulated\nInitial"),color="grey")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
axis.title.x=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
ggtitle("d")+scale_color_manual(values=colo)+
ylab('wNODF')

connectance=ggplot()+
geom_boxplot(data=empf2,aes(x="Empirical",y=wconnectance))+
geom_boxplot(data=subset(indf2,time==2000),aes(y=wconnectance,x="Simulated\nCoevolved",color=trait))+
geom_boxplot(data=subset(indf2,time==0),aes(y=wconnectance,x="Simulated\nInitial"),color="grey")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
axis.title.x=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
ggtitle("e")+scale_color_manual(values=colo)+
ylab('Weighted connectance')

modularity=ggplot()+
geom_boxplot(data=empf2,aes(x="Empirical",y=wconnectance))+
geom_boxplot(data=subset(indf2,time==2000),aes(y=modularity,x="Simulated\nCoevolved",color=trait))+
geom_boxplot(data=subset(indf2,time==0),aes(y=modularity,x="Simulated\nInitial"),color="grey")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
axis.title.x=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
ggtitle("f")+scale_color_manual(values=colo)+
ylab('Modularity')

top=plot_grid(plot_acp,pl2,ncol=2,align="hv")
bottom=plot_grid(dist_plot,nodf,connectance,modularity,ncol=4,align="hv")

grid.arrange(top,bottom,ncol=1)

png("Fig.S5.png",width=1200,height=800,res=160)

dev.off();
