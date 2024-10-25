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
const tf=2000

function simue(pini)
    alpha,r,epsilon,df,dive,competition,tf = pini
    nbsp_a=dive
    nbsp_p=dive
    traits=["morpho";"pheno"]
    uinit=[ones(Float64,nbsp_a+nbsp_p);df.mu_phen[1:dive*2];df.sd_phen[1:dive*2];]
    final=missing
    for trait in traits
        p=nbsp_a,nbsp_p,epsilon,alpha,competition,trait,r
        sol = @inbounds mDerivative2(uinit,p,tf)
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
















