# un seul run du recuit simule
function sa_un_run(X::Matrix{Float64}, k::Int, cfg::Dict, rng::AbstractRNG)
    n=size(X, 1)              # nb de points
    z =Vector{Int}(undef, n)  # solution de depart
    for i in 1:n
        z[i] = mod(i-1,k) + 1   # calcul baseline cyclique
    end
    inertie_courante, centres =calculer_inertie(X, z, k)  # inertie de depart

    f_best           =inertie_courante
    meilleurs_centres=copy(centres)
    T0, alpha= calibrer_T0(X,k,z,cfg,rng)   

    # si T0_auto est false on prend celui du config 
    if !get(cfg["sa"],"T0_auto",true)
        T0    =cfg["sa"]["T_init"]
        alpha =(cfg["sa"]["T_min"] / T0)^(1 / cfg["sa"]["n_iter"]) #recalcule avec nouv T0
    end
    T =T0
    # param lus dans config
    n_iter            =cfg["sa"]["n_iter"]             #nb max d iterations
    T_min             =cfg["sa"]["T_min"]              #temperature d arret
    patience          =cfg["convergence"]["patience"]  #iterations sans amelioration avant arret

    sans_amelioration =0              # compteur de stagnation
    premier_app       =false        
    # historiques: traçage courbes
    hist_best =Float64[]
    hist_cur  =Float64[]
    hist_T    =Float64[]

    for iter in 1:n_iter
        i =rand(rng, 1:n)  # point au hasard
        ancien_c =z[i]
        # calcul nb pts dans son cluster
        nb_membres =0
        for x in z
            if x==ancien_c
                nb_membres +=1
            end
        end
        # si cluster vide: on annule mais on note quand meme l'iteration
        if nb_membres ==1
            push!(hist_best, f_best)
            push!(hist_cur, inertie_courante)
            push!(hist_T, T)
            T =T*alpha      # temperature baisse
            continue
        end
        # recherce d'un nv cluster et retire tq le meme
        nouveau_c =rand(rng, 1:k)       
        while nouveau_c == ancien_c
            nouveau_c=rand(rng, 1:k)
        end

        z[i] = nouveau_c         # deplacement temporaire du pt
        inertie_apres, centres_apres = calculer_inertie(X, z, k) 
        delta = inertie_apres - inertie_courante   
        # critere de metropolis (accept à amelioration)
        if delta < 0||rand(rng) < exp(-delta/T)
            inertie_courante =inertie_apres  
            centres=centres_apres
            # nouveau meilleur resultat global
            if inertie_courante<f_best
                f_best            =inertie_courante
                meilleurs_centres =copy(centres)
                sans_amelioration =0
                premier_app=true
            elseif premier_app
                sans_amelioration+=1   # accepte mais pas de nouveau meilleur
            end
        else # mouvement refuse => remet l ancien cluster
            z[i] =ancien_c
            if premier_app
                sans_amelioration +=1
            end
        end
        # enregistrement des historiques 
        push!(hist_best, f_best)
        push!(hist_cur, inertie_courante)
        push!(hist_T, T)
        # arret anticipe
        if cfg["convergence"]["enabled"] && premier_app &&(sans_amelioration >=patience)
            break
        end
        T =T*alpha      # refroidissement geometrique
        if T<=T_min     
            break
        end
    end
    return meilleurs_centres, f_best, hist_best, hist_cur, hist_T
end

# algorithme du recuit simule avec runs independants
function sa_kmeans(X::Matrix{Float64}, k::Int, cfg::Dict; rng::AbstractRNG = Random.default_rng())
    n_runs = get(cfg["sa"], "n_runs", 1)  # nb run indep
    resultat = nothing    # meilleur run pour l instant
    for r in 1:n_runs
        res = sa_un_run(X, k, cfg, rng)
        # garde le run qui a la plus petite inertie 
        if resultat === nothing || res[2] < resultat[2]
            resultat = res
        end
    end
    return resultat
end