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

const liste = DataFrame(CSV.File(join(["/home/duchenne/pheno/initial_empir/liste.csv"])))
const tab=DataFrame(site=repeat(liste.site,10),essai=repeat([1:1:10;],17))
const site=tab.site[jj]
const essai=tab.essai[jj]
const nbsp_a=liste.na[liste.site .== site][1]
const nbsp_p=liste.np[liste.site .== site][1]

include("/home/duchenne/pheno/derivatives_function_m_pheno_or_morpho_discrete.jl")
const alpha=1
const r=-0.5
const epsilon=0.001

function simue(pini)
    alpha,r,epsilon,df,nbsp_a,nbsp_p,competition,tf = pini
    traits=["morpho";"pheno";]
    uinit=[ones(Float64,nbsp_a+nbsp_p);df.mu_phen[1:(nbsp_a+nbsp_p)];df.sd_phen[1:(nbsp_a+nbsp_p)];]
    final=missing
    for trait in traits
        p=nbsp_a,nbsp_p,epsilon,alpha,competition,trait,r
        sol = @inbounds mDerivative2(uinit,p,tf)
        if trait=="morpho"
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

for competition in [10;]
      println("essai",jj)
      tf=2000
      df = DataFrame(CSV.File(join(["/home/duchenne/pheno/initial_empir/pops_ini_",site,"_",essai,".csv"])))
      pini= alpha,r,epsilon,df,nbsp_a,nbsp_p,competition,tf
      finalf=simue(pini)
      GC.gc()
      CSV.write(join(["/home/duchenne/pheno/results_empir/ueq_",site,"_",essai,".csv"]),finalf)
end
  












