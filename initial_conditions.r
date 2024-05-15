#' **Generate initial conditions for models**
#'
#' François Duchenne
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)


nb_rep=500

for(i in 1:nb_rep){
nbsp_a=30
nbsp_p=30

duree=70
#pheno
mu_phen=rnorm(nbsp_a+nbsp_p,190,duree)
while(length(mu_phen[mu_phen>365 | mu_phen<0])>0){
mu_phen[mu_phen>365]=mu_phen[mu_phen>365]-365
mu_phen[mu_phen<0]=365+mu_phen[mu_phen<0]
}
mu_phen=runif(nbsp_a+nbsp_p,-2.5,2.5)
sd_phen=runif(nbsp_a+nbsp_p,-2.5,2.5)

#morpho
mu_morpho=rnorm(nbsp_a+nbsp_p,190,duree)
while(length(mu_morpho[mu_morpho>365 | mu_morpho<0])>0){
mu_morpho[mu_morpho>365]=mu_morpho[mu_morpho>365]-365
mu_morpho[mu_morpho<0]=365+mu_morpho[mu_morpho<0]
}
sd_morpho=runif(nbsp_a+nbsp_p,-2.5,2.5)


#vecteurs parametres
rmax <- runif(nbsp_a+nbsp_p,-0.5,-0.5)
K <- runif(nbsp_a+nbsp_p,100,100)
Nini <- runif(nbsp_a+nbsp_p,1,1)


#####SPECIES LEVEL INFORMATIONS TO EXPORT

final=data.frame(sp=c(paste0("a",1:nbsp_a),paste0("f",1:nbsp_p)),type=c(rep("poll",nbsp_a),rep("plant",nbsp_p)),mu_phen=mu_phen,
sd_phen=sd_phen,mu_morpho=mu_morpho,sd_morpho=sd_morpho,rmax=rmax,K=K,Nini=Nini)
final$random=i

setwd(dir=paste0("C:/Users/Duchenne/Documents/evolution_pheno_morpho/initial"))
fwrite(final,paste0("pops_ini_",i,".csv"))
}




