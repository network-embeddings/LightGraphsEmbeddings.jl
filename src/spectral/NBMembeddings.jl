include("matrices.jl")
include("../utils.jl")
function NBMembedding(g; matrix = ihara_matrix)
    edgeidmap, m, aristas = mapa(g)
    NBM, edgeidmap = matrix(g)
    vect = eigvecs(NBM)
    the_values = eigvals(NBM)
    the_real = real.(the_values)
    the_imag = imag.(the_values)
    the_norm = sqrt.((the_real.^2)+(the_imag.^2))
    both = hcat(the_imag,the_norm)
    sorted = sortslices(both, dims=1,by=x->x[2],rev=true)
    threshold = sorted[findall(x->x≠0,sorted[:,1]),2][1]
    hm, index = num_com(the_values,threshold)
    vectores = eigvecs(NBM)
    index = index[1:length(index)]
    matrix_embedded = real(vectores[:,index])
    cont = contraction(g,index,matrix_embedded,edgeidmap)
    M = fit(PCA, cont'; maxoutdim=2)
    contracted = transform(M, cont')
    contracted
end
