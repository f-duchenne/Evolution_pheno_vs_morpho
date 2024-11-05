ENV["JULIA_NUM_PRECOMPILE_TASKS"]=250

using Pkg; Pkg.instantiate()
using CSV
using DelimitedFiles
using DataFrames
using Infinity
using LinearAlgebra
using Statistics 
using ForwardDiff
using Distributions
using QuadGK


const jj=ARGS[1]
println("essai",jj)

include("C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts/functions and secondary scripts/derivatives_function_m_pheno_or_morpho_discrete.jl")
const alpha=1.0
const r=-0.5
const epsilon=0.01
const tf=3000

function simue(pini)
    alpha,r,epsilon,df,dive,competition,tf = pini
    nbsp_a=dive
    nbsp_p=dive
    traits=["morpho";"pheno"]
    uinit=[ones(Float64,nbsp_a+nbsp_p);df.mu_phen[1:dive*2];df.sd_phen[1:dive*2];]
    final=missing
    for trait in traits
        p=nbsp_a,nbsp_p,epsilon,alpha,competition,trait,r
        sol = @inbounds mDerivative(uinit,p,tf)
        if trait=="morpho"
            final=DataFrame(sol[:,(dive*2+1):end],:auto)
            final[!,"trait"] .=trait
        else
            bidon=DataFrame(sol[:,(dive*2+1):end],:auto)
            bidon[!,"trait"] .=trait
            final=vcat(final,bidon)
        end
    end
    return(final)
end

for competition in [2;5;]
  for dive in [10;20;30;]
      println("essai",jj)
      df = DataFrame(CSV.File(join(["C:/Users/Duchenne/Documents/evolution_pheno_morpho/initial/pops_ini_",jj,".csv"])))
      pini= alpha,r,epsilon,df,dive,competition,tf
      finalf=simue(pini)
      GC.gc()
      CSV.write(join(["/home/duchenne/pheno/results/ueq_",jj,"_",dive,"_",competition,".csv"]),finalf)
  end
end
 

##############
nbsp_a=dive
nbsp_p=dive
bidon=finalf[finalf.trait.=="morpho",:] 
sol=Matrix(bidon[:,1:(end-1)])
plot(sol[:,end],sol[:,(nbsp_a*1+nbsp_p*1+1):(nbsp_a*2+nbsp_p*2)], legend = false) #plot the mean traits


ENV["JULIA_NUM_PRECOMPILE_TASKS"]=250

using Pkg; Pkg.instantiate()
using CSV
using DelimitedFiles
using DataFrames
using Infinity
using LinearAlgebra
using Statistics 
using ForwardDiff
using Distributions
using QuadGK

const jj=parse(Int,ARGS[1])
println("essai",jj)

const liste = DataFrame(CSV.File(join(["/home/duchenne/pheno/initial_empir/liste.csv"])))
const tab=DataFrame(site=repeat(liste.site,10),essai=repeat([1:1:10;],17))
const sitef=tab[jj,"site"]
const essai=tab.essai[jj]
const nbsp_a=liste.na[liste.site .== sitef][1]
const nbsp_p=liste.np[liste.site .== sitef][1]
const tf=3000

include("/home/duchenne/pheno/derivatives_function_m_pheno_or_morpho_discrete_symmetric.jl")
include("/home/duchenne/pheno/derivatives_function_m_pheno_and_morpho_discrete_symmetric.jl")
const alpha=1
const r=-0.5
const epsilon=0.01

function simue(pini)
    alpha,r,epsilon,df,nbsp_a,nbsp_p,competition,tf = pini
    traits=["both";"pheno";]
    final=missing
    for trait in traits
        if trait=="both"
            uinit=[ones(Float64,nbsp_a+nbsp_p);df.mu_phen[1:dive*2];df.sd_phen[1:dive*2];reverse(df.mu_phen)[1:dive*2];reverse(df.sd_phen)[1:dive*2];]
            funcdev=mDerivative2
            p=nbsp_a,nbsp_p,epsilon,alpha,competition,trait,r
        else
            uinit=[ones(Float64,nbsp_a+nbsp_p);df.mu_phen[1:(nbsp_a+nbsp_p)];df.sd_phen[1:(nbsp_a+nbsp_p)];]
            funcdev=mDerivative
            p=nbsp_a,nbsp_p,epsilon,alpha,competition,r
        end
        sol = @inbounds funcdev(uinit,p,tf)
        if trait=="both"
            final=DataFrame(sol[:,((nbsp_a+nbsp_p)+1):end],:auto)
            final[!,"trait"] .=trait
        else
            bidon=DataFrame(sol[:,((nbsp_a+nbsp_p)+1):end],:auto)
            bidon[!,"trait"] .=trait
            final=vcat(final,bidon)
        end
    end
    return(final)
end

for competition in [5;]
      println("essai",jj)
      df = DataFrame(CSV.File(join(["/home/duchenne/pheno/initial_empir/pops_ini_",sitef,"_",essai,".csv"])))
      pini= alpha,r,epsilon,df,nbsp_a,nbsp_p,competition,tf
      finalf=simue(pini)
      GC.gc()
      CSV.write(join(["/home/duchenne/pheno/results_empir/ueq_",sitef,"_",essai,".csv"]),finalf)
end
  














