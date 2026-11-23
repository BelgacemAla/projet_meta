using Plots

# fonction pour tracer et sauvegarder ou afficher le graphique d'evolution de l'inertie
function afficher_graphique(hist_inertie::Vector{Float64})
    p = plot(hist_inertie, xlabel="iterations", ylabel="inertie", title="evolution du recuit simule",)
    display(p)
    # sauvegarde 
    savefig("output/evolution_inertie.png")
    println("Graphique sauvegarde dans output/evolution_inertie.png")
end

function afficher_convergence(hist_cur, hist_best, I_base)
    p = plot(hist_cur, label="f_cur", xlabel="iterations", ylabel="inertie")
    plot!(p, hist_best, label="f_best", color=:pink, linewidth=2)
    hline!(p, [I_base], label="baseline", color=:red)
    savefig(p,"output/convergence.png")
end

function afficher_temperature(hist_T)
    p = plot(hist_T, yscale=:log10, label="T", xlabel="iterations", ylabel="temperature")
    savefig(p,"output/temperature.png")
end