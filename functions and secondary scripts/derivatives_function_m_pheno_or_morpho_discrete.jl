@inbounds function mDerivative2(uinit,p,t)
    nbsp_a,nbsp_p,epsilon,alpha,competition,trait,r = p
    result = Array{Float64}(undef, t+1, length(uinit)+1)
    result[1,:]=[uinit;0;]
    u=copy(uinit)
    
    function invlogit(x)
        if x>=30.0
            return(50.0)
        end

        if x<=(-30.0)
            return(1.0)
        end

        if abs(x)<30.0
            return(1.0+49.0*exp(x)/(1.0+exp(x)))
        end
    end

    function invlogit1(x)
        return(80.0+205.0*exp(x)/(1.0+exp(x)))

    end

    function inte(theta,m1,s1,m2,s2)
        n1= pdf.(Normal.(m1,s1),theta)
        n2= pdf.(Normal.(m2,s2),theta)
        return(min.(n1,n2))
    end

    for ti in 1:t
        u2=copy(u)
        #Defining some shorcuts
        #pheno
        mu_phen_a = invlogit1.(u[(1+nbsp_a+nbsp_p):(nbsp_a*2+nbsp_p)])
        sd_phen_a = invlogit.(u[(1+nbsp_a*2+nbsp_p*2):(nbsp_a*3+nbsp_p*2)])
        mu_phen_p =invlogit1.(u[(1+nbsp_a*2+nbsp_p):(nbsp_a*2+nbsp_p*2)])
        sd_phen_p = invlogit.(u[(1+nbsp_a*3+nbsp_p*2):(nbsp_a*3+nbsp_p*3)])
        #abundances
        abund_poll= @view u[1:nbsp_a]
        abund_flower= @view u[(1+nbsp_a):(nbsp_a+nbsp_p)]
        #interaction matrix   
        p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a
        @inbounds function interaction_matrix(p2)
            mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a = p2
            m = Array{Float64}(undef, nbsp_p, nbsp_a)
            for j in 1:nbsp_a
                mu1= @view mu_phen_a[j]
                sd1 = @view sd_phen_a[j]
                for v in 1:nbsp_p
                    mu2= @view mu_phen_p[v]
                    sd2 = @view sd_phen_p[v]
                    phen=quadgk(theta -> inte.(theta,mu1,sd1,mu2,sd2), 0, 365, rtol=1e-6)[1]
                    m[v,j]=phen
                end  
            end
            return(m)
        end
        m = interaction_matrix(p2)

        #Loop defining derivatives for each species according to their guild
        for i=1:(nbsp_a+nbsp_p)

            af= @view u[i]

            @inbounds function pop_derivative_a(x,y)
                #functions for derivation
                p2b = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m
                @inbounds function mut_inter_a(p2b)
                    mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m = p2b
                    mut_interactions = Vector{}(undef, nbsp_p)
                    for v in 1:nbsp_p
                        mu2= @view mu_phen_p[v]
                        sd2= @view sd_phen_p[v]
                        phen=quadgk(theta -> inte.(theta,invlogit1.(x),invlogit.(y),mu2,sd2), 0, 365, rtol=1e-6)[1]
                        mut_interactions[v] = phen
                    end
                    return(mut_interactions)
                end
                mut_interactions = mut_inter_a(p2b)

                p3 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m,mut_interactions
                @inbounds function comp_inter_a(p3)
                    mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m,mut_interactions = p3
                    comp_interactions = Vector{}(undef, nbsp_a)
                    for v in 1:nbsp_a
                        if trait=="pheno"
                            mu2= @view mu_phen_a[v]
                            sd2= @view sd_phen_a[v]
                            phen=quadgk(theta -> inte.(theta,invlogit1.(x),invlogit.(y),mu2,sd2), 0, 365, rtol=1e-6)[1]
                        else
                            phen=1.0
                        end
                        competitor= @view m[:,v]
                        similarity=mean(sqrt.(mut_interactions .* competitor))
                        comp_interactions[v]=phen.*similarity
                    end
                    return(comp_interactions)
                end
                comp_interactions = comp_inter_a(p3)
                comp_interactions[i,] = 1.0
                d_pop= r .+ alpha*sum(mut_interactions .* abund_flower) .- competition*sum(comp_interactions .* abund_poll)
                return d_pop
            end

            @inbounds function pop_derivative_p(x,y)
                #functions for derivation
                p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m
                @inbounds function mut_inter_p(p2)
                    mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m = p2
                    mut_interactions = Vector{}(undef, nbsp_a)
                    for v in 1:nbsp_a
                        mu2= @view mu_phen_a[v]
                        sd2= @view sd_phen_a[v]
                        phen=quadgk(theta -> inte.(theta,invlogit1.(x),invlogit.(y),mu2,sd2), 0, 365, rtol=1e-6)[1]
                        mut_interactions[v] = phen
                    end
                    return(mut_interactions)
                end
                mut_interactions = mut_inter_p(p2)

                p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m,mut_interactions
                @inbounds function comp_inter_p(p2)
                    mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m,mut_interactions = p2
                    comp_interactions = Vector{}(undef, nbsp_p)
                    for v in 1:nbsp_p
                        if trait=="pheno"
                            mu2= @view mu_phen_p[v]
                            sd2= @view sd_phen_p[v]
                            phen=quadgk(theta -> inte.(theta,invlogit1.(x),invlogit.(y),mu2,sd2), 0, 365, rtol=1e-6)[1]
                        else
                            phen=1.0
                        end
                        competitor = @view m[v,:]
                        similarity=mean(sqrt.(mut_interactions .* competitor))
                        comp_interactions[v]=phen.*similarity
                    end
                    return(comp_interactions)
                end
                comp_interactions= comp_inter_p(p2)
                comp_interactions[(i-nbsp_a),] = 1.0
                d_pop= r .+ alpha*sum(mut_interactions .* abund_poll) .- competition*sum(comp_interactions .* abund_flower)
                return d_pop
            end

            ###################################################################################################

            if i< nbsp_a+1 #if pollinator
                funcdev=pop_derivative_a
            else#if plant
                funcdev=pop_derivative_p
            end
            ### EVOLUTION
            ### PARTIAL DERIVATIVES FOR EVOLUTION
            fit=funcdev(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
            fitval=[funcdev(u[(i+nbsp_a+nbsp_p)]+epsilon,u[(i+nbsp_a*2+nbsp_p*2)]);funcdev(u[(i+nbsp_a+nbsp_p)]-epsilon,u[(i+nbsp_a*2+nbsp_p*2)])]
            if(maximum(fitval)<=fit)
                partial_dev_mu_phen =0
            else
                if fitval[1] .> fitval[2]
                    partial_dev_mu_phen = epsilon
                else
                    partial_dev_mu_phen = -1 * epsilon
                end
            end
            u2[(nbsp_a+nbsp_p+i)]=u[(nbsp_a+nbsp_p+i)] + partial_dev_mu_phen

            fit=funcdev(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
            fitval=[funcdev(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)]+epsilon);funcdev(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)]-epsilon)]
            if(maximum(fitval)<=fit)
                partial_dev_sd_phen =0
            else
                if fitval[1] .> fitval[2]
                    partial_dev_sd_phen = epsilon
                else
                    partial_dev_sd_phen = -1 * epsilon
                end
            end
            u2[(nbsp_a*2+nbsp_p*2+i)]=u[(nbsp_a*2+nbsp_p*2+i)] + partial_dev_sd_phen

        end
        u=copy(u2)
        result[(ti+1),:]=[u;ti;]
    end
    return result[1:10:end,1:end]
end
