# calcul de l inertie et des centroides
function calculer_inertie(X::Matrix{Float64}, z::Vector{Int}, k::Int)
    n=size(X, 1)                # nb points
    d=size(X, 2)                # nb de dim
    centres  = zeros(k, d)      # une ligne par cluster
    compteurs=zeros(Int, k)     # nb de pts dans chaque cluster
    # calcul centroides ligne par ligne
    for i in 1:n
        c =z[i]
        compteurs[c] += 1
        for j in 1:d
            centres[c,j] +=X[i,j]
        end
    end
    # moy pour chaque centre
    for c in 1:k
        if compteurs[c]>0   #evite division par 0 si cluster vide
            for j in 1:d
                centres[c,j] =centres[c,j]/compteurs[c]
            end
        end
    end
    # calcule inertie totale
    inertie=0.0
    for i in 1:n
        c =z[i]
        somme_carres=0.0
        for j in 1:d          # distance au carre entre le point et son centre
            diff =X[i,j]-centres[c,j]
            somme_carres += diff^2
        end
        inertie +=somme_carres
    end
    return inertie,centres
end