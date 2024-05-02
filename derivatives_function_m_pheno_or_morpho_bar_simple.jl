@inbounds function mDerivative2(u,p,t)
    Kp,nbsp_a,nbsp_p,epsilon,seuil,θ,alpha,competition,trait,r,mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p = p    
    du = similar(u)

    function invlogit(x)
        # if x>=30
        #     return(50)
        # end

        # if x<=(-30)
        #     return(1)
        # end

        # if abs(x)<30
        #     return(1+50*exp(x)/(1+exp(x)))
        # end

        return(1+50*exp(x)/(1+exp(x)))

    end
    
    #Defining some shorcuts
    sd_phen_a=invlogit.(sd_phen_a)
    sd_phen_p=invlogit.(sd_phen_p)
    #abundances
    abund_poll= @view u[1:nbsp_a]
    abund_flower= @view u[(1+nbsp_a):(nbsp_a+nbsp_p)]
    #interaction matrix
    p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a
    @inbounds function interaction_matrix(p2)
        mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a = p2
        m = Array{Float64}(undef, nbsp_p, nbsp_a)
        for j in 1:nbsp_a
            mu1= @view mu_phen_a[j]
            sd1 = @view sd_phen_a[j]
            @inbounds function mut_inter(v)
                mu2= @view mu_phen_p[v]
                sd2 = @view sd_phen_p[v]
                phen1= mu1 .- 2.0 .* sd1..mu1 .+ 2.0 .* sd1
                phen2= mu2 .- 2.0 .* sd2..mu2 .+ 2.0 .* sd2
                phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                phen
            end
            m[:,j].=mut_inter.(1:nbsp_p)
        end
        return(m)
    end
    m=interaction_matrix(p2)

    #Loop defining derivatives for each species according to their guild
    for i=1:(nbsp_a+nbsp_p)
        af= @view u[i]
            
        if u[i] > seuil
            if i< nbsp_a+1 #if pollinator
                mu_phenf= @view mu_phen_a[i]
                sd_phenf= @view sd_phen_a[i]
                #interactions and competition vectors
                mut_interactions= @view m[:,i]
                @inbounds function comp_inter_a(v)
                    if trait=="pheno"
                        phen1= mu_phenf .- 2.0 .* sd_phenf..mu_phenf .+ 2.0 .* sd_phenf
                        phen2= mu_phen_a[v] .- 2.0 .* sd_phen_a[v]..mu_phen_a[v] .+ 2.0 .* sd_phen_a[v]
                        phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                    else
                        phen=1.0
                    end
                    competitor= @view m[:,v]
                    similarity=sum(mut_interactions .* competitor) / sum(mut_interactions)
                    phen.*similarity
                end
                comp_interactions= comp_inter_a.(1:nbsp_a)
                comp_interactions[i,] = 0.0
                ### DERIVATIVES FOR ABUNDANCE
                du[i]= af .* (r .- af./Kp .+ alpha*sum(mut_interactions .* abund_flower) /(1+sum(mut_interactions .* abund_flower)) .- competition*sum(comp_interactions .* abund_poll))

            else#if plant
                mu_phenf= @view mu_phen_p[i-nbsp_a]
                sd_phenf= @view sd_phen_p[i-nbsp_a]
                #interactions and competition vectors
                mut_interactions= @view m[(i-nbsp_a),:]
                @inbounds function comp_inter_p(v)
                    if trait=="pheno"
                        phen1= mu_phenf .- 2.0 .* sd_phenf..mu_phenf .+ 2.0 .* sd_phenf
                        phen2= mu_phen_p[v] .- 2.0 .* sd_phen_p[v]..mu_phen_p[v] .+ 2.0 .* sd_phen_p[v]
                        phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                    else
                        phen=1.0
                    end
                    competitor = @view m[v,:]
                    similarity=sum(mut_interactions .* competitor) / sum(mut_interactions)
                    phen.*similarity
                end
                comp_interactions= comp_inter_p.(1:nbsp_p)
                comp_interactions[i-nbsp_a,] = 0.0  
                ### DERIVATIVES FOR ABUNDANCE
                du[i]= af .* (r .- af./Kp .+ alpha*sum(mut_interactions .* abund_poll) /(1+sum(mut_interactions .* abund_poll))  .- competition*sum(comp_interactions .* abund_flower))
                
            end
        else
            du[i] = 0
        end

    end
    return du
end
