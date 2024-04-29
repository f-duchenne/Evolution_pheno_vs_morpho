ENV["JULIA_NUM_PRECOMPILE_TASKS"]=250

using Pkg; Pkg.instantiate()
using CSV
using DelimitedFiles
using DataFrames
using Infinity
using OrdinaryDiffEq
using LinearAlgebra
using Statistics 
using ForwardDiff
using Distributions

const θ = range(0.0, 365.0, length=180)
const conv=1e-10
const seuil=1e-5
const maxiter=30000
const alg = AutoTsit5(Rosenbrock23())
const tspan=(0.0,0.1)
const maxdiv=30
const threat=0

const jj=1 #ARGS[1]
println("essai",jj)

const df = DataFrame(CSV.File(join(["C:/Users/Duchenne/Documents/evolution_pheno_morpho/initial/pops_ini_",jj,".csv"])))

include("C:/Users/Duchenne/Documents/evolution_pheno_morpho/scripts/derivatives_function_m_pheno_or_morpho.jl")

dive=30
alpha=0.5
competition=0.5
nbsp_a=dive
nbsp_p=dive
Kp=1
r=df.rmax
uinit=[ones(Float64,nbsp_a+nbsp_p);df.mu_phen[1:dive*2];df.sd_phen[1:dive*2];df.mu_morpho[1:dive*2];df.sd_morpho[1:dive*2];]
trait="morpho"
p=Kp,nbsp_a,nbsp_p,epsilon,seuil,θ,alpha,competition,morpho
prob = ODEProblem(mDerivative,uinit,(0.0,0.1),p)

sol = @inbounds solve(prob, alg,saveat=1)
plot(sol.t,transpose(sol[(1):(nbsp_a+nbsp_p),1:end]))
plot(sol.t,transpose(sol[(nbsp_a*2+nbsp_p*2+1):(nbsp_a*3+nbsp_p*3),1:end]))
plot(sol.t,sol[(nbsp_a+nbsp_p+1),1:end])

function solve_model(problem,algo,tspan,p)
    ac=1
    b=0
    sol = @inbounds solve(problem, algo,saveat=1)
    b=b+maximum(tspan)
    ac=maximum(var.(eachrow(sol[:,(end-9):end])))
    while ac > conv && b <= maxiter
        problem = ODEProblem(mDerivative,sol[:,end],tspan,p)
        sol = @inbounds solve(problem, algo, saveat = 1)
        b=b+maximum(tspan)
        ac=maximum(var.(eachrow(sol[:,(end-9):end])))
    end
    return sol,ac
end

obj= @inbounds solve_model(prob,alg,tspan,p)


@inbounds function simulations(jj) 
    GC.gc()

    for epsilon=[1.0;]
        println(epsilon)
        for dive=[10:10:maxdiv;]
            println(dive)
            for kind=["mutualism"]
                println(kind)

                alpha=1.0
                competition=1.0
                nbsp_a=30
                nbsp_p=30
                Kp=100
                r=df.rmax
                uinit=[ones(Float64,nbsp_a+nbsp_p)*10;df.mu_phen;df.sd_phen;df.mu_morpho;df.sd_morpho;]
                p=Kp,nbsp_a,nbsp_p,epsilon,seuil,θ,alpha,competition
                prob = ODEProblem(mDerivative,uinit,tspan,p)

                function solve_model(problem,algo,tspan,p)
                    ac=1
                    b=0
                    sol = @inbounds solve(problem, algo,saveat=1)
                    b=b+maximum(tspan)
                    ac=maximum(var.(eachrow(sol[:,(end-9):end])))
                    while ac > conv && b <= maxiter
                        problem = ODEProblem(mDerivative,sol[:,end],tspan,p)
                        sol = @inbounds solve(problem, algo, saveat = 1)
                        b=b+maximum(tspan)
                        ac=maximum(var.(eachrow(sol[:,(end-9):end])))
                    end
                    return sol,ac
                end

                obj= @inbounds solve_model(prob,alg,tspan,p)
                sol=obj[1]
                ac=obj[2]

                Species=Vector{Union{Missing, String}}(missing, (nbsp_p+nbsp_f))
                for i=1:(nbsp_p+nbsp_f)
                    Species[i]=join(["sp",i])
                end

                final1=DataFrame(Species=[repeat(Species,3);Species[(nbsp_p+1):(nbsp_p+nbsp_f)]],guild=[repeat(vcat(repeat(["I"],nbsp_p),repeat(["F"],nbsp_f)),3);repeat(["F"],nbsp_f)],
                type=[repeat(["abundance";"trait_mean";"trait_sd"],inner=(nbsp_p+nbsp_f));repeat(["resistance"],nbsp_f)],
                Nini=uinit,N_eq1=sol[:,end],rzero=[repeat(r,3);r[(nbsp_p+1):(nbsp_p+nbsp_f)]],
                alpha=repeat([alpha],(dive*3*2+nbsp_f)),essai=repeat([jj],(dive*3*2+nbsp_f)),conv=repeat([ac],(dive*3*2+nbsp_f)),Kp=repeat([Kp],(dive*3*2+nbsp_f)),epsilon=repeat([epsilon],(dive*3*2+nbsp_f)),
                kind=repeat([kind],(dive*3*2+nbsp_f)),dive=repeat([dive],(dive*3*2+nbsp_f)))

                final=vcat(final,final1)
                final1=0
                GC.gc()
            end
        end
    end
    return final
end

finalf= @inbounds simulations(jj)

CSV.write(join(["/home/duchenne/ER/eq1/initial_",jj,".csv"]),finalf[(2:size(finalf,1)),:])
