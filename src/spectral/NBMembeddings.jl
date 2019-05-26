include(matrices.jl)
function NBMembedding(g; matrix = ihara_matrix)
    edgeidmap, m, aristas = mapa(g)
    NBM, edgeidmap = matrix(g)
    vect = eigvecs(NBM)
    index = 1:nv(g) |> collect
    contracted = contraction(g,index,real(vect[:,index]),edgeidmap)
    contracted
end
