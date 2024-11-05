###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","viridis","glmmTMB","lme4") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
empf=fread("empirical_networks.csv")
names(empf)[1:2]=paste0("w",names(empf)[1:2])
names(empf)[3]="inter. evenness"
names(empf)[5]="wNODF"
colo=c("chartreuse3","dodgerblue4","gold3")

empf$dive=empf$na+empf$np 
empf$prop_motif=(empf$motifn/empf$motifntot)
#model=lmer(feas_with_pheno~rho*(prop_motif+dive)+(1|site),data=subset(empf,competition==5 & rho %in% c(0.01,0.05,0.2)))

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info_empir.csv")
names(indf)=gsub("weighted ","w",names(indf))
names(indf)[3]="inter. evenness"
names(indf)[5]="wNODF"
names(indf)[7]="mut.strength"

indf$dive=indf$nbsp_a+indf$nbsp_p 
indf$prop_motif=(indf$motifn/indf$motifntot)
#model=glmmTMB(feas~(prop_motif+dive)*trait+(1|site),family=beta_family(link="logit"),data=subset(indf,competition==10 & rho==0.05 & time==2000))

indf2=subset(indf,rho==0.01)
acp=PCA(indf2[,1:6],graph=FALSE)
indf2=cbind(indf2,acp$ind$coord)

empf=cbind(empf,predict(acp,newdata=empf[,1:6])$coord)
empf$trait2="empirical"

indf2$trait2=indf2$trait
indf2$trait2[indf2$time==0]="initial"

plot_acp=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf2,time==2000),aes(x=Dim.1,y=Dim.2,color=trait2))+
geom_point(data=subset(indf2,time==0),aes(x=Dim.1,y=Dim.2,color=trait2),alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="bottom",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=c("black","grey",colo))+ggtitle("a")+geom_point(data=empf,aes(x=Dim.1,y=Dim.2,color=trait2))+
labs(color="")

pl2=fviz_pca_var(acp, col.var = "black", repel =TRUE)+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold",hjust = 0))

png("Fig.S5.png",width=1200,height=800,res=160)
plot_grid(plot_acp,pl2,ncol=2,align="hv")
dev.off();


empf2=subset(empf,competition==5 & rho==0.05 & trait=="both")
pl1=ggplot()+
geom_boxplot(data=empf2,aes(x="Empirical",y=motifn/motifntot))+
geom_boxplot(data=subset(indf,rho==0.05 & time==2000),aes(y=motifn/motifntot,x="Simulated\nCoevolved"),color=colo[1])+
geom_boxplot(data=subset(indf,rho==0.05 & time==0),aes(y=motifn/motifntot,x="Simulated\nInitial"),color="grey")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
axis.title.x=element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
ggtitle("a")+
ylab('Percentage of "V+" motifs')+scale_y_continuous(labels=scales::percent)


empf2=subset(empf,competition==5)
empf2$dive=empf2$na+empf2$np
b=melt(empf2,id.vars=c("dive","rho"),measure.vars=c("feas_with_pheno","feas_without_pheno"))
b$variable=as.character(b$variable)
b$variable[b$variable=="feas_without_pheno"]="without seasonal structure"
b$variable[b$variable=="feas_with_pheno"]="with seasonal structure"

pl2=ggplot(data=subset(b,rho==0.05),aes(x=dive,y=value,color=variable,fill=variable))+
geom_point(size=2,alpha=0.9)+
stat_smooth(method="glm",alpha=0.2,,method.args = list(family = quasibinomial(link = 'logit')))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position = c(0.75, 0.9))+scale_fill_manual(values=c("black","grey"))+
scale_color_manual(values=c("black","grey"))+ggtitle("b")+coord_cartesian()+
ylab(expression(paste("Structural stability  ",(omega))))+
xlab("Diversity")+labs(color="",fill="")


pl3=ggplot(data=subset(empf,competition %in% c(2,5)),aes(x=as.factor(competition),y=feas_without_pheno/feas_with_pheno-1,color=as.factor(rho)))+
geom_hline(linetype="dashed",yintercept=0)+
geom_hline(linetype="dashed",yintercept=-1)+
geom_boxplot()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border=element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+ggtitle("c")+
scale_y_continuous(labels=scales::percent)+ylab("Decrease in structural stability\nwhen neglecting seasonal structure")+
xlab("Competition strength")+scale_colour_manual(values=rev(scales::viridis_pal(option="rocket")(6))[1:3])+labs(color=expression(rho))


plot_grid(pl1,pl2,pl3,ncol=3,align="hv",rel_widths=c(0.75,2,1.5))


pdf("Fig.5.pdf",width=10,height=5)
plot_grid(pl1,pl2,pl3,ncol=3,align="hv",rel_widths=c(1.4,2.5,1.5))
dev.off();


empf2=subset(empf,competition==5 & rho==0.2)

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