@inbounds function mDerivative(u,p,t)
    Kp,nbsp_a,nbsp_p,epsilon,seuil,θ,alpha,competition,morpho = p    
    du = similar(u)

    function invlogit(x)
        if x>=30
            return(50)
        end

        if x<=(-30)
            return(1)
        end

        if abs(x)<30
            return(1+50*exp(x)/(1+exp(x)))
        end

    end
    
    #Defining some shorcuts
    #pheno
    mu_phen_a = @view u[(1+nbsp_a+nbsp_p):(nbsp_a*2+nbsp_p)]
    sd_phen_a = invlogit.(u[(1+nbsp_a*2+nbsp_p*2):(nbsp_a*3+nbsp_p*2)])
    mu_phen_p = @view u[(1+nbsp_a*2+nbsp_p):(nbsp_a*2+nbsp_p*2)]
    sd_phen_p = invlogit.(u[(1+nbsp_a*2+nbsp_p*2):(nbsp_a*3+nbsp_p*3)])
    #morpho
    mu_morpho_a = @view u[(1+nbsp_a*3+nbsp_p*3):(nbsp_a*4+nbsp_p*3)]
    sd_morpho_a = invlogit.(u[(1+nbsp_a*4+nbsp_p*4):(nbsp_a*5+nbsp_p*4)])
    mu_morpho_p = @view u[(1+nbsp_a*4+nbsp_p*3):(nbsp_a*4+nbsp_p*4)]
    sd_morpho_p = invlogit.(u[(1+nbsp_a*5+nbsp_p*4):(nbsp_a*5+nbsp_p*5)])
    #abundances
    abund_poll= @view u[1:nbsp_a]
    abund_flower= @view u[(1+nbsp_a):(nbsp_a+nbsp_p)]
    #interaction matrix
    p2 = mu_phen_a,sd_phen_a,mu_morpho_a,sd_morpho_a,mu_phen_p,sd_phen_p,mu_morpho_p,sd_morpho_p,θ,nbsp_p,nbsp_a
    @inbounds function interaction_matrix(p2)
        mu_phen_a,sd_phen_a,mu_morpho_a,sd_morpho_a,mu_phen_p,sd_phen_p,mu_morpho_p,sd_morpho_p,θ,nbsp_p,nbsp_a = p2
        m = Array{Float64}(undef, nbsp_p, nbsp_a)
        for j in 1:nbsp_a
            @inbounds function mut_inter(v) 
                phen=sum(min.(pdf.(Normal(mu_phen_a[j],sd_phen_a[j]),θ), pdf.(Normal.(mu_phen_p[v],sd_phen_p[v]),θ)))*step(θ)
                morpho=sum(min.(pdf.(Normal(mu_morpho_a[j],sd_morpho_a[j]),θ), pdf.(Normal.(mu_morpho_p[v],sd_morpho_p[v]),θ)))*step(θ)
                phen.*morpho
            end
            m[:,j].=mut_inter.(1:nbsp_p)
        end
        return(m)
    end
    m=interaction_matrix(p2)

    #Loop defining derivatives for each species according to their guild
    for i=1:(nbsp_a+nbsp_p)
        af= @view u[i]
        @inbounds function pop_derivative_a(x,y,w,z)
            #functions for derivation
            @inbounds function mut_inter_a(v)
                mu1= @view mu_phen_p[v]
                sd1= @view sd_phen_p[v]
                mu2= @view mu_morpho_p[v]
                sd2= @view sd_morpho_p[v]
                phen=sum(pdf.(Normal(x,invlogit.(y)),θ) .* pdf.(Normal.(mu1,sd1),θ))*step(θ)
                morpho=sum(pdf.(Normal(w,invlogit.(z)),θ) .* pdf.(Normal.(mu2,sd2),θ))*step(θ)
                phen.*morpho
            end
            mut_interactions = mut_inter_a.(1:nbsp_p)

            @inbounds function comp_inter_a(v)
                mu1= @view mu_phen_a[v]
                sd1= @view sd_phen_a[v]
                phen=sum(min.(pdf.(Normal(x,invlogit.(y)),θ), pdf.(Normal.(mu1,sd1),θ)))*step(θ)
                competitor = @view m[:,v]
                similarity=sum(m[:,i] .* competitor)
                phen*similarity
            end
            comp_interactions = comp_inter_a.(1:nbsp_a)
            comp_interactions[i,] = 1.0
            r=-0.6 #+(exp(-(abs(x-182)/180)^3)+exp(-(abs(w-182)/180)^3))/4
            d_pop= r .+ alpha*sum(mut_interactions .* abund_flower) .- competition*sum(comp_interactions .* abund_poll)
            return d_pop
        end
        ∂f_∂x_a(x, y, w, z) = @inbounds ForwardDiff.derivative(x -> pop_derivative_a(x,y,w,z),x)
        ∂f_∂y_a(x, y, w, z) = @inbounds ForwardDiff.derivative(y -> pop_derivative_a(x,y,w,z),y)
        ∂f_∂w_a(x, y, w, z) = @inbounds ForwardDiff.derivative(w -> pop_derivative_a(x,y,w,z),w)
        ∂f_∂z_a(x, y, w, z) = @inbounds ForwardDiff.derivative(z -> pop_derivative_a(x,y,w,z),z)
    
        @inbounds function pop_derivative_p(x,y,w,z)
            #functions for derivation
            @inbounds function mut_inter_p(v)
                mu1= @view mu_phen_a[v]
                sd1= @view sd_phen_a[v]
                mu2= @view mu_morpho_a[v]
                sd2= @view sd_morpho_a[v] 
                phen=sum(min.(pdf.(Normal.(x,invlogit.(y)),θ), pdf.(Normal.(mu1,sd1),θ)))*step(θ)
                morpho=sum(min.(pdf.(Normal.(w,invlogit.(z)),θ), pdf.(Normal.(mu2,sd2),θ)))*step(θ)
                phen.*morpho
            end
            mut_interactions = mut_inter_p.(1:nbsp_a)

            @inbounds function comp_inter_p(v)
                mu1= @view mu_phen_p[v]
                sd1= @view sd_phen_p[v]
                phen=sum(min.(pdf.(Normal.(x,invlogit.(y)),θ), pdf.(Normal.(mu1,sd1),θ)))*step(θ)
                competitor = @view m[v,:]
                similarity=sum(m[(i-nbsp_a),:] .* competitor)
                phen.*similarity
            end
            comp_interactions = comp_inter_p.(1:nbsp_p)
            comp_interactions[(i-nbsp_a),] = 1.0
            r=-0.6 #+(exp(-(abs(x-182)/180)^3)+exp(-(abs(w-182)/180)^3))/4
            d_pop= r .+ alpha*sum(mut_interactions .* abund_poll) .- competition*sum(comp_interactions .* abund_flower)
            return d_pop
        end
        ∂f_∂x_p(x, y, w, z) = @inbounds ForwardDiff.derivative(x -> pop_derivative_p(x,y,w,z),x)
        ∂f_∂y_p(x, y, w, z) = @inbounds ForwardDiff.derivative(y -> pop_derivative_p(x,y,w,z),y)
        ∂f_∂w_p(x, y, w, z) = @inbounds ForwardDiff.derivative(w -> pop_derivative_p(x,y,w,z),w)
        ∂f_∂z_p(x, y, w, z) = @inbounds ForwardDiff.derivative(z -> pop_derivative_p(x,y,w,z),z)

        if u[i] > seuil
            if i< nbsp_a+1 #if pollinator
                mu_phenf= @view mu_phen_a[i]
                sd_phenf= invlogit.(sd_phen_a[i])
                #interactions and competition vectors
                mut_interactions= @view m[:,i]
                @inbounds function comp_inter_a(v)
                    phen=sum(min.(pdf.(Normal.(mu_phenf,sd_phenf),θ), pdf.(Normal.(mu_phen_a[v],sd_phen_a[v]),θ)))*step(θ)
                    competitor= @view m[:,v]
                    similarity=sum(mut_interactions .* competitor)
                    phen.*similarity
                end
                comp_interactions= comp_inter_a.(1:nbsp_a)
                r=-0.6 #+(exp(-(abs(mu_phenf.-182)/180)^3)+exp(-(abs(mu_morphof.-182)/180)^3))/4
                ### DERIVATIVES FOR ABUNDANCE
                du[i]= af .* (r .+ alpha*sum(mut_interactions .* abund_flower) .- competition*sum(comp_interactions .* abund_poll))

                ### PARTIAL DERIVATIVES FOR EVOLUTION
                partial_dev_mu_phen = @inbounds ∂f_∂x_a.(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)],u[(i+nbsp_a*3+nbsp_p*3)],u[(i+nbsp_a*4+nbsp_p*4)])
                partial_dev_sd_phen = @inbounds ∂f_∂y_a.(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)],u[(i+nbsp_a*3+nbsp_p*3)],u[(i+nbsp_a*4+nbsp_p*4)])
                if morpho==1
                    partial_dev_mu_morpho = @inbounds ∂f_∂w_a.(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)],u[(i+nbsp_a*3+nbsp_p*3)],u[(i+nbsp_a*4+nbsp_p*4)])
                    partial_dev_sd_morpho = @inbounds ∂f_∂z_a.(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)],u[(i+nbsp_a*3+nbsp_p*3)],u[(i+nbsp_a*4+nbsp_p*4)])
                else
                    partial_dev_mu_morpho = 0.0
                    partial_dev_sd_morpho = 0.0
                end

                ### EVOLUTION
                du[(nbsp_a+nbsp_p+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_mu_phen)
                du[(nbsp_a*2+nbsp_p*2+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_sd_phen)
                du[(nbsp_a*3+nbsp_p*3+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_mu_morpho)
                du[(nbsp_a*4+nbsp_p*4+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_sd_morpho)

            else#if plant
                mu_phenf= @view mu_phen_a[i-nbsp_a]
                sd_phenf= invlogit.(sd_phen_a[i-nbsp_a])
                #interactions and competition vectors
                mut_interactions= @view m[(i-nbsp_a),:]
                @inbounds function comp_inter_p(v)
                    phen=sum(min.(pdf.(Normal.(mu_phenf,sd_phenf),θ), pdf.(Normal.(mu_phen_a[v],sd_phen_a[v]),θ)))*step(θ)
                    competitor = @view m[v,:]
                    similarity=sum(mut_interactions .* competitor)
                    phen.*similarity
                end
                comp_interactions= comp_inter_p.(1:nbsp_a)
                r=-0.6 #+(exp(-(abs(mu_phenf.-182)/180)^3)+exp(-(abs(mu_morphof.-182)/180)^3))/4
                ### DERIVATIVES FOR ABUNDANCE
                du[i]= af .* (r .+ alpha*sum(mut_interactions .* abund_poll) .- competition*sum(comp_interactions .* abund_flower))
                
                ### PARTIAL DERIVATIVES FOR EVOLUTION

                partial_dev_mu_phen = @inbounds ∂f_∂x_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)],u[(i+nbsp_a*3+nbsp_p*3)],u[(i+nbsp_a*4+nbsp_p*4)])
                partial_dev_sd_phen = @inbounds ∂f_∂y_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)],u[(i+nbsp_a*3+nbsp_p*3)],u[(i+nbsp_a*4+nbsp_p*4)])
                if morpho==1
                    partial_dev_mu_morpho = @inbounds ∂f_∂w_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)],u[(i+nbsp_a*3+nbsp_p*3)],u[(i+nbsp_a*4+nbsp_p*4)])
                    partial_dev_sd_morpho = @inbounds ∂f_∂z_p(u[(i+nbsp_a+nbsp_p)],u[(i+nbsp_a*2+nbsp_p*2)],u[(i+nbsp_a*3+nbsp_p*3)],u[(i+nbsp_a*4+nbsp_p*4)])
                else
                    partial_dev_mu_morpho = 0.0
                    partial_dev_sd_morpho = 0.0
                end
                #### EVOLUTION
                du[(nbsp_a+nbsp_p+i)]= (sqrt(af .+ 1) .* epsilon .* partial_dev_mu_phen)
                du[(nbsp_a*2+nbsp_p*2+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_sd_phen)
                du[(nbsp_a*3+nbsp_p*3+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_mu_morpho)
                du[(nbsp_a*4+nbsp_p*4+i)]=(sqrt(af .+ 1) .* epsilon .* partial_dev_sd_morpho)

            end
        else
            du[i] = 0
            du[(nbsp_a+nbsp_p+i)] = 0
            du[(nbsp_a*2+nbsp_p*2+i)] = 0
            du[(nbsp_a*3+nbsp_p*3+i)] = 0
            du[(nbsp_a*4+nbsp_p*4+i)] = 0
        end

    end
    return du
end
