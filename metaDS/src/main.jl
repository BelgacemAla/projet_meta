using Random
using TOML
include("charger.jl")
include("inertie.jl")
include("calibrage.jl")
include("recuit.jl")
include("affichage.jl")

function main()
    cfg =TOML.parsefile("config.toml")    # chargement config
    # liste des fichiers test et leur k
    fichiers =["data/donnees_didactique.csv", "data/instance2_R2_grand.csv", "data/instance2b_R2_densites.csv", "data/instance3_iris.csv", "data/instance4_wine.csv", "data/instance5_haute_dim.csv"]
    ks       =[5, 6, 5, 3, 3, 7]

    # affichage du menu
    for i in 1:length(fichiers)
        println(i, " : ", fichiers[i])
    end
    print("choix : ")
    choix =parse(Int,readline())
    
    k         =ks[choix]
    data_file =fichiers[choix]

    X= charger_donnees(data_file)
    rng=MersenneTwister(42)

    # baseline cyclique, pour comparer avec le SA
    z_base =[mod(i-1,k) + 1 for i in 1:size(X,1)]
    I_base, _ =calculer_inertie(X, z_base, k)
    println("Baseline : $I_base")

    # lancement et affichage
    centres, best_inertie, hist_best, hist_cur,hist_T = sa_kmeans(X, k, cfg; rng=rng)
    println("Fichier : $data_file (k = $k)")
    println("Meilleure inertie trouvee : $best_inertie")
    println("Iterations : ",length(hist_best))

    # on s'assure que le dossier output existe
    mkpath("output")
    afficher_convergence(hist_cur,hist_best,I_base)
    afficher_temperature(hist_T)
    println("Graphes sauvegardes dans dossier output")
end

main()