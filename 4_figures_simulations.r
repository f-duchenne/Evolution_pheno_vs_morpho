###################################################################################
###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2","bipartite","FactoMineR","factoextra","gridExtra","cowplot","ggpubr","scales","piecewiseSEM",
"igraph","qgraph","car","nlme","viridis") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

colo=c("dodgerblue4","chartreuse3")
invlogit=function(x){return(1+49*exp(x)/(1+exp(x)))}
invlogit1=function(x){return(80+205*exp(x)/(1+exp(x)))}

##################################### SPECIES LEVEL ########################
setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
datf=fread("species_level_data.csv")
endpo=subset(datf,time==2000 & comp==10)


sp1=ggplot()+
geom_density(data=subset(endpo,type=="mu"),aes(x=invlogit1(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("density")+
xlab("Mean value of the trait")+ggtitle("a")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_grid(. ~rich, labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")

sp2=ggplot()+
geom_density(data=subset(endpo,type=="sd"),aes(x=invlogit(value),color=trait,fill=trait),alpha=0.2)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),
legend.key = element_rect(fill = "white"))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+ylab("density")+
xlab("Standard deviation of the trait")+ggtitle("b")+labs(linetype=expression(paste(n[sp]," / guild")))+
guides(linetype=guide_legend(override.aes = list(fill = "white",colour="black")))+facet_wrap( ~rich, labeller = label_bquote(cols=n[sp] == .(rich)), scales = "free_y")

plot_grid(sp1,sp2,ncol=1,align="hv")

png("Fig.S2.png",width=1200,height=800,res=150)
plot_grid(sp1,sp2,ncol=1,align="hv")
dev.off();

##################################### COMMUNITY LEVEL ########################

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho")
indf=fread("C:/Users/Duchenne/Documents/evolution_pheno_morpho/networks_info.csv")
names(indf)=gsub("weighted ","",names(indf))
names(indf)[3]="inter. evenness"

indf2=subset(indf,competition==10 & rho==0.01)
acp=PCA(indf2[,1:6],graph=FALSE)
indf2=cbind(indf2,acp$ind$coord)

pl1=ggplot(data=subset(indf2,rich==10),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,group=paste(essai,trait)))+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point()+
geom_line(show_guide =F)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("a")+coord_fixed(ratio=1)

pl2=fviz_pca_var(acp, col.var = "black", repel =TRUE)+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold",hjust = 0))

pl3=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf2,time==2000),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,shape=as.factor(rich)))+
geom_point(data=subset(indf2,time==0),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("c")

pl4=ggplot(data=subset(indf2,time==2000),aes(x=as.factor(rich),y=NODF,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Nestedness (NODF)")+xlab(expression(paste(n[sp], " / guild")))+
ggtitle("d")

pl5=ggplot(data=subset(indf2,time==2000),aes(x=as.factor(rich),y=connectance,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Connectance")+xlab(expression(paste(n[sp], " / guild")))+ggtitle("e")

top=plot_grid(pl1,pl2,align="hv",ncol=2,rel_widths=c(1.5,1))
bottom=plot_grid(pl3,pl4,pl5,ncol=3,rel_widths=c(2,1,1))

pdf("Fig.2.pdf",width=7,height=7)
grid.arrange(top,bottom,ncol=1)
dev.off();

#################### FIG. S3
indf3=subset(indf,competition==5 & rho==0.01)
indf3=cbind(indf3,predict(acp,newdata=indf3[,1:6])$coord)

pl3=ggplot()+
geom_vline(xintercept=0,linetype="dashed")+
geom_hline(yintercept=0,linetype="dashed")+
geom_point(data=subset(indf3,time==2000),aes(x=Dim.1,y=Dim.2,alpha=time,color=trait,shape=as.factor(rich)))+
geom_point(data=subset(indf3,time==0),aes(x=Dim.1,y=Dim.2,alpha=time,shape=as.factor(rich)),color="grey",alpha=0.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+
xlab(paste0("Dimension 1 (",round(acp$eig[1,2],digits=1)," %)"))+
ylab(paste0("Dimension 2 (",round(acp$eig[2,2],digits=1)," %)"))+
scale_color_manual(values=colo)+ggtitle("a")


pl4=ggplot(data=subset(indf3,time==2000),aes(x=as.factor(rich),y=NODF,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Nestedness (NODF)")+xlab(expression(paste(n[sp], " / guild")))+
ggtitle("b")


pl5=ggplot(data=subset(indf3,time==2000),aes(x=as.factor(rich),y=connectance,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width")+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(N[A] == N[P]))+ylab("Connectance")+xlab(expression(paste(n[sp], " / guild")))+ggtitle("c")

bottom=plot_grid(pl3,pl4,pl5,ncol=3,rel_widths=c(2,1,1))

png("Fig.S3.png",width=1200,height=600,res=150)
bottom
dev.off();

#########################################

indf4=subset(indf,competition==10)

pl6=ggplot(data=subset(indf4,rho==0.05),aes(x=as.factor(time),y=feas,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width",show.legend=FALSE)+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab("Size of the feasibility domain")+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_discrete(breaks=c(0,2000))+ggtitle("a",subtitle="low interspecific competition")

pl6b=ggplot(data=subset(indf4,rho==0.4),aes(x=as.factor(time),y=feas,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width",show.legend=FALSE)+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab("Size of the feasibility domain")+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_discrete(breaks=c(0,2000))+ggtitle("b",subtitle="high interspecific competition")

plot_grid(pl6,pl6b,ncol=1,align="hv")

pdf("Fig.3.pdf",width=6,height=6)
plot_grid(pl6,pl6b,ncol=1,align="hv")
dev.off();

############################## FIGURE S4
indf4=subset(indf,competition==5)

pl6=ggplot(data=subset(indf4,rho==0.05),aes(x=as.factor(time),y=feas,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width",show.legend=FALSE)+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.background = element_blank(),
strip.background=element_rect(fill=NA,color=NA))+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab("Size of the feasibility domain")+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_discrete(breaks=c(0,2000))+ggtitle("a",subtitle="low interspecific competition")

pl6b=ggplot(data=subset(indf4,rho==0.4),aes(x=as.factor(time),y=feas,color=trait))+geom_violin(position=position_dodge(width=0.5),width=1,scale="width",show.legend=FALSE)+
geom_boxplot(width=0.1,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
scale_color_manual(values=colo)+scale_fill_manual(values=colo)+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab("Size of the feasibility domain")+xlab("Time")+facet_grid(cols=vars(rich), labeller = label_bquote(cols=n[sp] == .(rich)))+scale_x_discrete(breaks=c(0,2000))+ggtitle("b",subtitle="high interspecific competition")

plot_grid(pl6,pl6b,ncol=1,align="hv")
png("Fig.S4.png",width=1200,height=1000,res=150)
plot_grid(pl6,pl6b,ncol=1,align="hv")
dev.off();

############################# FIGURE 5
indf5=subset(indf4,time %in% c(0,2000))

pl7=ggplot(data=subset(indf5,time==0),aes(x=rich,y=feas,color=as.factor(rho)))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
stat_summary(fun.y=mean, geom="line", size = 0.5,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab("Size of the feasibility domain")+xlab("Number of species per guild")+facet_wrap(~trait)+scale_colour_viridis(option="D",discrete=TRUE)+
labs(colour=expression(rho))+ggtitle("a",subtitle="Initial networks")

pl8=ggplot(data=subset(indf5,time==2000),aes(x=rich,y=feas,color=as.factor(rho)))+
stat_summary(fun.y = mean, fun.ymin = function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom = "pointrange")+
stat_summary(fun.y=mean, geom="line", size = 0.5,position=position_dodge(width=0.5))+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank())+
labs(linetype=expression(paste(n[sp]," / guild")))+ylab("Size of the feasibility domain")+xlab("Number of species per guild")+facet_wrap(~trait)+scale_colour_viridis(option="D",discrete=TRUE)+
labs(colour=expression(rho))+ggtitle("b",subtitle="Coevolved networks")


pdf("Fig.4.pdf",width=6,height=6)
plot_grid(pl7,pl8,ncol=1,align="hv")
dev.off();

indf %>% group_by(rich,trait,time) %>% summarise(mean(feas))




indf4$NODFl=scale(logit(indf4$NODF/100,adjust=0.0001))
indf4$connectancel=scale(logit(indf4$connectance,adjust=0.0001))
indf4$modularityl=scale(logit(indf4$modularity,adjust=0.0001))
indf4$feasl=scale(logit(indf4$feas,adjust=0.0001))
indf4$compl=scale(logit(indf4$comp,adjust=0.0001))
indf4$richl=scale(indf4$rich)


pdf("Fig.4.pdf",width=12,height=7)
par(mfrow=c(2,2))
for(ti in c(0,2000)){
for(tr in c("morpho","pheno")){
indfb=as.data.frame(subset(indf4,trait==tr & time==ti & rho==0.05))

modelco=lme(connectancel~richl,data=indfb,random=~1|essai,control = lmeControl(opt = "optim"))
suma1=as.data.frame(summary(modelco)$tTable)[-1,]
suma1$Predictor=rownames(suma1)
suma1$Response="connectancel"

modelnest=lme(NODFl~richl+connectancel,data=indfb,random=~1|essai,control = lmeControl(opt = "optim"))
print(car::vif(modelnest))
suma2=as.data.frame(summary(modelnest)$tTable)[-1,]
suma2$Predictor=rownames(suma2)
suma2$Response="NODFl"

modelcomp=lme(modularityl~richl+connectancel,data=indfb,random=~1|essai,control = lmeControl(opt = "optim"))
suma3=as.data.frame(summary(modelcomp)$tTable)[-1,]
suma3$Predictor=rownames(suma3)
suma3$Response="modularityl"

print(car::vif(modelcomp))
modelstab=lme(feasl~(NODFl+connectancel+richl+modularityl),data=indfb,random=~1|essai,control = lmeControl(opt = "optim"))
suma4=as.data.frame(summary(modelstab)$tTable)[-1,]
suma4$Predictor=rownames(suma4)
suma4$Response="feasl"
print(car::vif(modelstab))

#obj=piecewiseSEM::psem(modelco,modelnest,modelcomp,modelstab,data=indf)
objb=rbind(suma1,suma2,suma3,suma4) #summary(obj)
names(objb)[1]="Std.Estimate"
left=0
right=30
center=(right-left)/2
haut=20
bas=0
cex=1.3
echelle_fleche=3
cex.arrows=0.7
l=objb[,c("Predictor","Response","Std.Estimate","p-value")]

l$colo=ifelse(l$Std.Estimate<0,"firebrick4","dodgerblue4")
l$lty=1
l$lty=ifelse(l$"p-value">0.05,2,l$lty)
l$curv=NA
l$curv[l$Predictor=="richl" & l$Response=="feasl"]=2.7
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
vertex_attr(g, "name") 
g= g %>% set_vertex_attr("name", value =
c("Diversity",
paste0("Connectance\n(R2 = ",round(rsq::rsq.lmm(modelco)$fixed,digits=2),")"),
paste0("Nestedness\n(R2 = ",round(rsq::rsq.lmm(modelnest)$fixed,digits=2),")"),
paste0("Modularity\n(R2 = ",round(rsq::rsq.lmm(modelcomp)$fixed,digits=2),")"),
paste0("Feasibility\n(R2 = ",round(rsq::rsq.lmm(modelstab)$fixed,digits=2),")")))
coord=data.frame(label=vertex_attr(g, "name"),
x=c(center,center,right,left,center),
y=c(haut,bas+11,bas+11,bas+11,bas),vsize=25)

EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])*echelle_fleche
asi[asi<5]=5
asi[asi>=15]=15

title=paste0(ifelse(tr=="morpho","a - Morphological trait","b - Phenological trait"),", time = ",ti)

qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=cex,label.scale=F,lty=l$lty,
edge.label.cex = cex.arrows*2,edge.label.position=0.60,vsize2=9,
vsize=coord$vsize,curve=l$curv,
shape="rectangle",edge.labels=T,fade=F,asize=asi,
mar=c(5,5,5,5),knot.border.color="white",curveShape=-1,title=title)
}

}
dev.off();
#



