function NBMembedding(g; matrix = ollin_matrix)
    edgeidmap, m, aristas = mapa(g)
    NBM, edgeidmap = matrix(g)
    the_values = eigvals(NBM)
    the_real = real.(the_values)
    the_imag = imag.(the_values)
    the_norm = sqrt.((the_real.^2)+(the_imag.^2))
    both = hcat(the_imag,the_norm)
    sorted = sortslices(both, dims=1,by=x->x[2],rev=true)
    threshold = sorted[findall(x->x≠0,sorted[:,1]),2][1]
    hm, index = num_com(the_values,threshold)
    if length(hm) > 1
        vect = eigvecs(NBM)
        index = index[2:length(index)]
        matriz_embedded = real(vect[:,index])
        contracted = contraction(g,index,matriz_embedded,edgeidmap)
    else
        vect= eigvecs(NBM)
        index = index[2:length(index)]
        matriz_embedded = real(vect[:,index])
        contracted = contraction(g,index,matriz_embedded,edgeidmap)
    end
    contracted
end # module
