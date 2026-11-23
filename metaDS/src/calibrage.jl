using Statistics

# calibrage automatique de t0
function calibrer_T0(X::Matrix{Float64}, k::Int, z_init::Vector{Int}, cfg::Dict, rng::AbstractRNG)
    deltas = Float64[]      # gardeer les degradations
    z = copy(z_init)        # pour pas toucher solut de depart
    n = length(z)
    # mouv aleatoire 
    for iter in 1:500
        i=rand(rng,1:n)     # un pt au hasard
        ancien_c =z[i]
        nb_membres=0
        # comptage nb de points dans ce cluster
        for x in z
            if x==ancien_c
                nb_membres +=1
            end
        end
        if nb_membres==1
            continue
        end
        
        nouveau_c =rand(rng, 1:k)
        while nouveau_c == ancien_c       # on retire tant que c'est le meme cluster
            nouveau_c =rand(rng, 1:k)
        end
        
        z[i] =nouveau_c
        I_av, _=calculer_inertie(X,z_init,k)   # inertie avant
        I_ap, _=calculer_inertie(X,z,k)        # apres le mouvement
        delta =I_ap-I_av
        # on garde que les mauvais mouvements (ceux qui augmentent l'inertie)
        if delta >0                  
            push!(deltas,delta)  
        end 
        z[i] = ancien_c     
    end
    # aucune degradation trouve ->met 1.0 pour pas planter
    delta_moyen =isempty(deltas) ? 1.0 : mean(deltas)
    taux = get(cfg["sa"],"taux_cible",0.50)
    T0= -delta_moyen/log(taux)    
    N =cfg["sa"]["n_iter"]
    T_min =cfg["sa"]["T_min"]
    alpha=(T_min/T0)^(1 /N)

    return T0, alpha
end