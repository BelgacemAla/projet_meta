using DelimitedFiles
# chargement des donnees
function charger_donnees(chemin::String)
    M =readdlm(chemin,',',comments=true) #ignorer commentaires #...
    # nleve colonne = label_vrai
    return Float64.(M[:,1:end-1])
end