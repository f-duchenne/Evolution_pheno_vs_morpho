@inbounds function mDerivative2(uinit,p,t)
    nbsp_a,nbsp_p,epsilon,alpha,competition,trait,r = p
    result = Array{Float64}(undef, t+1, length(uinit)+1)
    result[1,:]=[uinit;0;]
    u=copy(uinit)
    preci=0.000001
    max_gen_var=5
    for ti in 1:t
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
                    function inte(theta)
                        phen1= pdf.(Normal.(mu1,sd1),theta)
                        phen2= pdf.(Normal.(mu2,sd2),theta)
                        return(min.(phen1,phen2))
                    end
                    phen=quadgk(inte, 0, 365, rtol=1e-6)[1]
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
                        mu1= @view mu_phen_p[v]
                        sd1= @view sd_phen_p[v]
                        function inte(theta)
                            phen1= pdf.(Normal.(invlogit1.(x),invlogit.(y)),theta)
                            phen2= pdf.(Normal.(mu1,sd1),theta)
                            return(min.(phen1,phen2))
                        end
                        phen=quadgk(inte, 0, 365, rtol=1e-6)[1]
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
                            mu1= @view mu_phen_a[v]
                            sd1= @view sd_phen_a[v]
                            function inte(theta)
                                phen1= pdf.(Normal.(invlogit1.(x),invlogit.(y)),theta)
                                phen2= pdf.(Normal.(mu1,sd1),theta)
                                return(min.(phen1,phen2))
                            end
                            phen=quadgk(inte, 0, 365, rtol=1e-6)[1]
                        else
                            phen=1.0
                        end
                        competitor= @view m[:,v]
                        similarity=sum(mut_interactions .* competitor) / (sum(mut_interactions))
                        comp_interactions[v]=phen.*similarity
                    end
                    return(comp_interactions)
                end
                comp_interactions = comp_inter_a(p3)
                comp_interactions[i,] = 0.0
                d_pop= r .- af .+ alpha*sum(mut_interactions .* abund_flower) .- competition*sum(comp_interactions .* abund_poll)
                return d_pop
            end
            ∂f_∂x_a(x, y) = @inbounds ForwardDiff.derivative(x -> pop_derivative_a(x,y),x)
            ∂f_∂y_a(x, y) = @inbounds ForwardDiff.derivative(y -> pop_derivative_a(x,y),y)

            @inbounds function pop_derivative_p(x,y)
                #functions for derivation
                p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m
                @inbounds function mut_inter_p(p2)
                    mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,nbsp_p,nbsp_a,m = p2
                    mut_interactions = Vector{}(undef, nbsp_a)
                    for v in 1:nbsp_a
                        mu1= @view mu_phen_a[v]
                        sd1= @view sd_phen_a[v]
                        function inte(theta)
                            phen1= pdf.(Normal.(invlogit1.(x),invlogit.(y)),theta)
                            phen2= pdf.(Normal.(mu1,sd1),theta)
                            return(min.(phen1,phen2))
                        end
                        phen=quadgk(inte, 0, 365, rtol=1e-6)[1]
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
                            mu1= @view mu_phen_p[v]
                            sd1= @view sd_phen_p[v]
                            function inte(theta)
                                phen1= pdf.(Normal.(invlogit1.(x),invlogit.(y)),theta)
                                phen2= pdf.(Normal.(mu1,sd1),theta)
                                return(min.(phen1,phen2))
                            end
                            phen=quadgk(inte, 0, 365, rtol=1e-6)[1]
                        else
                            phen=1.0
                        end
                        competitor = @view m[v,:]
                        similarity=sum(mut_interactions .* competitor) / (sum(mut_interactions))
                        comp_interactions[v]=phen.*similarity
                    end
                    return(comp_interactions)
                end
                comp_interactions= comp_inter_p(p2)
                comp_interactions[(i-nbsp_a),] = 0.0
                d_pop= r .- af .+ alpha*sum(mut_interactions .* abund_poll) .- competition*sum(comp_interactions .* abund_flower)
                return d_pop
            end
            ∂f_∂x_p(x, y) = @inbounds ForwardDiff.derivative(x -> pop_derivative_p(x,y),x)
            ∂f_∂y_p(x, y) = @inbounds ForwardDiff.derivative(y -> pop_derivative_p(x,y),y)

            ###################################################################################################

            if i< nbsp_a+1 #if pollinator
                ### EVOLUTION
                ### PARTIAL DERIVATIVES FOR EVOLUTION
                fit=pop_derivative_a(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                vec=[pop_derivative_a(u[(i+nbsp_a+nbsp_p)]+preci,u[(i+nbsp_a*2+nbsp_p*2)]);pop_derivative_a(u[(i+nbsp_a+nbsp_p)]-preci,u[(i+nbsp_a*2+nbsp_p*2)])]
                if(maximum(vec)<=fit)
                    partial_dev_mu_phen =0
                else
                    partial_dev_mu_phen = vec[1].-vec[2] #@inbounds ∂f_∂x_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                    partial_dev_mu_phen = (1/preci)*partial_dev_mu_phen
                    if abs(partial_dev_mu_phen)>max_gen_var
                        partial_dev_mu_phen = sign(partial_dev_mu_phen)*max_gen_var
                    end
                end
                u[(nbsp_a+nbsp_p+i)]=u[(nbsp_a+nbsp_p+i)]+(sqrt(af .+ 1) .* epsilon .* partial_dev_mu_phen)

                fit=pop_derivative_a(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                vec=[pop_derivative_a(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)]+preci);pop_derivative_a(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)]-preci)]
                if(maximum(vec)<=fit)
                    partial_dev_sd_phen =0
                else
                    partial_dev_sd_phen = vec[1].-vec[2]#@inbounds ∂f_∂y_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                    partial_dev_sd_phen = (1/preci)*partial_dev_sd_phen
                    if abs(partial_dev_sd_phen)>max_gen_var
                        partial_dev_sd_phen = sign(partial_dev_sd_phen)*max_gen_var
                    end
                end
                u[(nbsp_a*2+nbsp_p*2+i)]=u[(nbsp_a*2+nbsp_p*2+i)]+(sqrt(af .+ 1) .* epsilon .* partial_dev_sd_phen)

            else#if plant
                #### EVOLUTION                
                ### PARTIAL DERIVATIVES FOR EVOLUTION
                fit=pop_derivative_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                vec=[pop_derivative_p(u[(i+nbsp_a+nbsp_p)]+preci,u[(i+nbsp_a*2+nbsp_p*2)]);pop_derivative_p(u[(i+nbsp_a+nbsp_p)]-preci,u[(i+nbsp_a*2+nbsp_p*2)])]
                if(maximum(vec)<=fit)
                    partial_dev_mu_phen =0
                else
                    partial_dev_mu_phen = vec[1].-vec[2] #@inbounds ∂f_∂x_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                    partial_dev_mu_phen = (1/preci)*partial_dev_mu_phen
                    if abs(partial_dev_mu_phen)>max_gen_var
                        partial_dev_mu_phen = sign(partial_dev_mu_phen)*max_gen_var
                    end
                end
                u[(nbsp_a+nbsp_p+i)]=u[(nbsp_a+nbsp_p+i)]+(sqrt(af .+ 1) .* epsilon .* partial_dev_mu_phen)

                fit=pop_derivative_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                vec=[pop_derivative_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)]+preci);pop_derivative_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)]-preci)]
                if(maximum(vec)<=fit)
                    partial_dev_sd_phen =0
                else
                    partial_dev_sd_phen = vec[1].-vec[2]#@inbounds ∂f_∂y_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                    partial_dev_sd_phen = (1/preci)*partial_dev_sd_phen
                    if abs(partial_dev_sd_phen)>max_gen_var
                        partial_dev_sd_phen = sign(partial_dev_sd_phen)*max_gen_var
                    end
                end
                u[(nbsp_a*2+nbsp_p*2+i)]=u[(nbsp_a*2+nbsp_p*2+i)]+(sqrt(af .+ 1) .* epsilon .* partial_dev_sd_phen)
            end
        end
        result[(ti+1),:]=[u;ti;]
    end
    return result[1:10:end,1:end]
end
