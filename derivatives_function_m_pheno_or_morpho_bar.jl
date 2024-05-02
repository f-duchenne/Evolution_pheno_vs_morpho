@inbounds function mDerivative2(u,p,t)
    Kp,nbsp_a,nbsp_p,epsilon,seuil,θ,alpha,competition,trait,r,h,eco = p    
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

    function invlogit1(x)
        return(80+285*exp(x)/(1+exp(x)))

    end
    
    
    #Defining some shorcuts
    #pheno
    mu_phen_a = invlogit1.(u[(1+nbsp_a+nbsp_p):(nbsp_a*2+nbsp_p)])
    sd_phen_a = invlogit.(u[(1+nbsp_a*2+nbsp_p*2):(nbsp_a*3+nbsp_p*2)])
    mu_phen_p =invlogit1.(u[(1+nbsp_a*2+nbsp_p):(nbsp_a*2+nbsp_p*2)])
    sd_phen_p = invlogit.(u[(1+nbsp_a*2+nbsp_p*2):(nbsp_a*3+nbsp_p*3)])
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
            for v in 1:nbsp_p
                mu2= @view mu_phen_p[v]
                sd2 = @view sd_phen_p[v]
                phen1= mu1 .- 2.0 .* sd1..mu1 .+ 2.0 .* sd1
                phen2= mu2 .- 2.0 .* sd2..mu2 .+ 2.0 .* sd2
                phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
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
            p2b = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m
            @inbounds function mut_inter_a(p2b)
                mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m = p2b
                mut_interactions = Vector{}(undef, nbsp_p)
                for v in 1:nbsp_p
                    mu1= @view mu_phen_p[v]
                    sd1= @view sd_phen_p[v]
                    phen1= invlogit1.(x) .- 2.0 .* invlogit.(y)..invlogit1.(x) .+ 2.0 .* invlogit.(y)
                    phen2= mu1 .- 2.0 .* sd1..mu1 .+ 2.0 .* sd1
                    phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                    mut_interactions[v] = phen
                end
                return(mut_interactions)
            end
            mut_interactions = mut_inter_a(p2b)

            p3 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m,mut_interactions
            @inbounds function comp_inter_a(p3)
                mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m,mut_interactions = p3
                comp_interactions = Vector{}(undef, nbsp_a)
                for v in 1:nbsp_a
                    if trait=="pheno"
                        phen1= invlogit1.(x) .- 2.0 .* invlogit.(y)..invlogit1.(x) .+ 2.0 .* invlogit.(y)
                        phen2= mu_phen_a[v] .- 2.0 .* sd_phen_a[v]..mu_phen_a[v] .+ 2.0 .* sd_phen_a[v]
                        phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                    else
                        phen=1.0
                    end
                    competitor= @view m[:,v]
                    similarity=sum(mut_interactions .* competitor) / sum(mut_interactions)
                    comp_interactions[v]=phen.*similarity
                end
                return(comp_interactions)
            end
            comp_interactions = comp_inter_a(p3)
            comp_interactions[i,] = 0.0
            d_pop= r .- af./Kp .+ alpha*sum(mut_interactions .* abund_flower) / (1+h*sum(mut_interactions .* abund_flower)) .- competition*sum(comp_interactions .* abund_poll)
            return d_pop
        end
        ∂f_∂x_a(x, y) = @inbounds ForwardDiff.derivative(x -> pop_derivative_a(x,y),x)
        ∂f_∂y_a(x, y) = @inbounds ForwardDiff.derivative(y -> pop_derivative_a(x,y),y)

    
        @inbounds function pop_derivative_p(x,y)
            #functions for derivation
            p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m
            @inbounds function mut_inter_p(p2)
                mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m = p2
                mut_interactions = Vector{}(undef, nbsp_a)
                for v in 1:nbsp_a
                    mu1= @view mu_phen_a[v]
                    sd1= @view sd_phen_a[v]
                    phen1= invlogit1.(x) .- 2.0 .* invlogit.(y)..invlogit1.(x) .+ 2.0 .* invlogit.(y)
                    phen2= mu1 .- 2.0 .* sd1..mu1 .+ 2.0 .* sd1
                    phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                    mut_interactions[v] = phen
                end
                return(mut_interactions)
            end
            mut_interactions = mut_inter_p(p2)

            p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m,mut_interactions
            @inbounds function comp_inter_p(p2)
                mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m,mut_interactions = p2
                comp_interactions = Vector{}(undef, nbsp_p)
                for v in 1:nbsp_p
                    if trait=="pheno"
                        phen1= invlogit1.(x) .- 2.0 .* invlogit.(y)..invlogit1.(x) .+ 2.0 .* invlogit.(y)
                        phen2= mu_phen_p[v] .- 2.0 .* sd_phen_p[v]..mu_phen_p[v] .+ 2.0 .* sd_phen_p[v]
                        phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                    else
                        phen=1.0
                    end
                    competitor = @view m[v,:]
                    similarity=sum(mut_interactions .* competitor) / sum(mut_interactions)
                    comp_interactions[v]=phen.*similarity
                end
                return(comp_interactions)
            end
            comp_interactions= comp_inter_p(p2)
            comp_interactions[(i-nbsp_a),] = 0.0
            d_pop= r .- af./Kp .+ alpha*sum(mut_interactions .* abund_poll) / (1+h*sum(mut_interactions .* abund_poll)) .- competition*sum(comp_interactions .* abund_flower)
            return d_pop
        end
        ∂f_∂x_p(x, y) = @inbounds ForwardDiff.derivative(x -> pop_derivative_p(x,y),x)
        ∂f_∂y_p(x, y) = @inbounds ForwardDiff.derivative(y -> pop_derivative_p(x,y),y)

        ###################################################################################################

        if u[i] > seuil
            if i< nbsp_a+1 #if pollinator
                if eco > 0.0
                    mu_phenf= @view mu_phen_a[i]
                    sd_phenf= @view sd_phen_a[i]
                    #interactions and competition vectors
                    mut_interactions= @view m[:,i]
                    p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m,mut_interactions
                    @inbounds function comp_inter_a(p2)
                        mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m,mut_interactions = p2
                        comp_interactions = Array{Float64}(undef, nbsp_a, 1)
                        for v in 1:nbsp_a
                            if trait=="pheno"
                                phen1= mu_phenf .- 2.0 .* sd_phenf..mu_phenf .+ 2.0 .* sd_phenf
                                phen2= mu_phen_a[v] .- 2.0 .* sd_phen_a[v]..mu_phen_a[v] .+ 2.0 .* sd_phen_a[v]
                                phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                            else
                                phen=1.0
                            end
                            competitor= @view m[:,v]
                            similarity=sum(mut_interactions .* competitor) / sum(mut_interactions)
                            comp_interactions[v]=phen.*similarity
                        end
                        return(comp_interactions)
                    end
                    comp_interactions= comp_inter_a(p2)
                    comp_interactions[i,] = 0.0
                    ### DERIVATIVES FOR ABUNDANCE
                    du[i]= eco .* af .* (r .- af./Kp .+ alpha *sum(mut_interactions .* abund_flower) / (1+h*sum(mut_interactions .* abund_flower)) .- competition*sum(comp_interactions .* abund_poll))
                else
                    du[i]=0.0
                end
                ### PARTIAL DERIVATIVES FOR EVOLUTION
                partial_dev_mu_phen = @inbounds ∂f_∂x_a.(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                partial_dev_sd_phen = @inbounds ∂f_∂y_a.(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])


                ### EVOLUTION
                du[(nbsp_a+nbsp_p+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_mu_phen)
                du[(nbsp_a*2+nbsp_p*2+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_sd_phen)


            else#if plant
                if eco > 0
                    mu_phenf= @view mu_phen_p[i-nbsp_a]
                    sd_phenf= @view sd_phen_p[i-nbsp_a]
                    #interactions and competition vectors
                    mut_interactions= @view m[(i-nbsp_a),:]
                    p2 = mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m,mut_interactions
                    @inbounds function comp_inter_p(p2)
                        mu_phen_a,sd_phen_a,mu_phen_p,sd_phen_p,θ,nbsp_p,nbsp_a,m,mut_interactions = p2
                        comp_interactions = Array{Float64}(undef, nbsp_a, 1)
                        for v in 1:nbsp_p
                            if trait=="pheno"
                                phen1= mu_phenf .- 2.0 .* sd_phenf..mu_phenf .+ 2.0 .* sd_phenf
                                phen2= mu_phen_p[v] .- 2.0 .* sd_phen_p[v]..mu_phen_p[v] .+ 2.0 .* sd_phen_p[v]
                                phen=Intervals.span(intersect(phen1,phen2))/max(Intervals.span(phen1),Intervals.span(phen2))
                            else
                                phen=1.0
                            end
                            competitor = @view m[v,:]
                            similarity=sum(mut_interactions .* competitor) / sum(mut_interactions)
                            comp_interactions[v]=phen.*similarity
                        end
                        return(comp_interactions)
                    end
                    comp_interactions= comp_inter_p(p2)
                    comp_interactions[i-nbsp_a,] = 0.0  
                    ### DERIVATIVES FOR ABUNDANCE
                    phen1= mu_phenf .- 2.0 .* sd_phenf..mu_phenf .+ 2.0 .* sd_phenf
                    du[i]= eco .* af .* (r .- af./Kp .+ alpha*sum(mut_interactions .* abund_poll) /(1+h*sum(mut_interactions .* abund_poll)) .- competition*sum(comp_interactions .* abund_flower))
                else
                    du[i]=0.0
                end
                                
                ### PARTIAL DERIVATIVES FOR EVOLUTION

                partial_dev_mu_phen = @inbounds ∂f_∂x_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])
                partial_dev_sd_phen = @inbounds ∂f_∂y_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)])

                #### EVOLUTION
                du[(nbsp_a+nbsp_p+i)]= (sqrt(af .+ 1) .* epsilon .* partial_dev_mu_phen)
                du[(nbsp_a*2+nbsp_p*2+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_sd_phen)
            end
        else
            du[i] = 0
            du[(nbsp_a+nbsp_p+i)] = 0
            du[(nbsp_a*2+nbsp_p*2+i)] = 0
        end

    end
    return du
end
