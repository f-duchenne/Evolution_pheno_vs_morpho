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
axis.text.x = element_text(angle = 45, hjust=1),
strip.background=element_rect(fill=NA,color=NA))+
ggtitle("a")+
ylab('Percentage of "V+" motifs')+scale_y_continuous(labels=scales::percent)

### STRUCTURAL STABILITY
empf$dive=empf$na+empf$np
b=melt(empf,id.vars=c("dive","competition"),measure.vars=c("feas_with_pheno","feas_without_pheno"))
b$variable=as.character(b$variable)
b$variable[b$variable=="feas_without_pheno"]="without seasonal structure"
b$variable[b$variable=="feas_with_pheno"]="with seasonal structure"


model=glm(value~dive*variable*as.factor(competition),family=quasibinomial,data=b)
car::Anova(model)

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
xlab(expression(paste("Competition\nimportance ",(italic(c)))))


plot_grid(pl1,pl2,pl3,ncol=3,align="hv",rel_widths=c(1,2,1))


pdf("Figures/Fig.5.pdf",width=10,height=5)
plot_grid(pl1,pl2,pl3,ncol=3,align="hv",rel_widths=c(0.75,2,0.75))
dev.off();
