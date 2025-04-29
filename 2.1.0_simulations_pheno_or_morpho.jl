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

include("/home/duchenne/pheno/derivatives_function_pheno_or_morpho.jl")
const alpha=1.0
const r=0
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


for competition in [2;4;6;]
  for dive in [10;20;30;]
      println("essai",jj)
      df = DataFrame(CSV.File(join(["/home/duchenne/pheno/initial/pops_ini_",jj,".csv"])))
      pini= alpha,r,epsilon,df,dive,competition,tf
      finalf=simue(pini)
      GC.gc()
      CSV.write(join(["/home/duchenne/pheno/results_symmetric/ueq_",jj,"_",dive,"_",competition,".csv"]),finalf)
  end
end
  