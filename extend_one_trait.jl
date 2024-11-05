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
const tf=1000


function simue(pini)
    alpha,r,epsilon,df,dive,competition,tf = pini
    nbsp_a=dive
    nbsp_p=dive
    traits=["morpho";"pheno"]
    final=df
    for trait in traits
        df2=filter(:trait => ==(trait), df)
        uinit=[ones(Float64,nbsp_a+nbsp_p);collect(df2[end,1:(end-2)]);]
        p=nbsp_a,nbsp_p,epsilon,alpha,competition,trait,r
        sol = @inbounds mDerivative2(uinit,p,tf)
        bidon=DataFrame(sol[2:end,(dive*2+1):end],:auto)
        bidon[:,end] = bidon[:,end] .+ maximum(df2[:,(end-1)])
        bidon[!,"trait"] .=trait
        final=vcat(final,bidon)
    end
    return(final)
end


for competition in [2;5;]
  for dive in [10;20;30;]
      println("essai",jj)
      df = DataFrame(CSV.File(join(["C:/Users/Duchenne/Documents/evolution_pheno_morpho/results_new/ueq_",jj,"_",dive,"_",competition,".csv"])))
      pini= alpha,r,epsilon,df,dive,competition,tf
      finalf=simue(pini)
      GC.gc()
      CSV.write(join(["/home/duchenne/pheno/results/ueq_",jj,"_",dive,"_",competition,".csv"]),finalf)
  end
end
  












