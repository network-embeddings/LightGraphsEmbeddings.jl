#!/usr/bin/env julia

function kron_δ(i,j)
    if i == j
        return 1
    else
        return 0
    end
end

function harmean(i,j)
    2 / ((1/i)+(1/j))
end

function inv_degree(g,i)
    1/degree(g,i)
end

function strength(g,i)
    neigh = neighbors(g,i)
    strengthi = sum(map(x-> g.weights[x,i],neigh))
    strengtho = sum(map(x-> g.weights[i,x],neigh))
    strength = strengthi + strengtho
    strength
end

function strengthin(g,i)
    neigh = inneighbors(g,i)
    strength = sum(map(x-> g.weights[x,i],neigh))
    strength
end

function strengthout(g,i)
    neigh = outneighbors(g,i)
    strength = sum(map(x-> g.weights[x,i],neigh))
    strength
end

function inv_degree_ollin(g,i,v)
    try
        1/strength(g,i,v)
    catch
        0
    end
end

function rel_single(g,i)
    if outdegree(g,i) == 1
        1
    else
        0
    end
end

function mglobal_clustering_coefficient(g::AbstractGraph)
    c = 0
    ntriangles = 0
    for v in vertices(g)
        neighs = neighbors(g, v)
        for i in neighs, j in neighs
            i == j && continue
            if has_edge(g, i, j)
                c += 1
            end
        end
        k = degree(g, v)
        ntriangles += k * (k - 1)
    end
    ntriangles == 0 && return 0.
    return c / ntriangles
end

function inv_degree_N(g,i)
    if degree(g,i) != 1
        1/(degree(g,i)-1)
    else
        1/degree(g,i)
    end
end

function nonBM_embedding(g)
  NB = la_NBM(g)
  λ, eigve = eigs(NB, nev=20)
  los_reales = Array{Float64}(0)
  los_indices_reales = Array{Int64}(0)
  for i in 1:length(λ)
        if imag(λ[i]) == 0 && real(λ[i]) != 0 && i*real(λ[i]) >= sqrt(real(λ[1]))
          push!(los_reales,λ[i])
          push!(los_indices_reales,i)
      end
  end
  matriz_embedded = sub(real(eigve),:,los_indices_reales)
  length(los_reales), los_indices_reales
  ϕ = zeros(Float64, nv(g), length(los_reales))
  for n=1:(length(los_reales))
      v= matriz_embedded[:,n]
      ϕ[:,n] = contract(g, v)
  end
  return length(los_reales), ϕ
end

function grados(g)
	grados = Array{Int64}(0)
	gradosdiv = Array{Int64}(0)
	for i in vertices(g)
		push!(grados,degree(g,i))
	end
	grados
end

function avg_grados(g)
    grado = grados(g)
    avg_grad = mean(grado)
    avg_grad
end

function max_grados(g)
    grado = grados(g)
    maxi = maximum(grado)
    maxi
end

function mat_an(g,v,NBM)
	unos = ones(2ne(g))/sqrt(2ne(g))
	unnos = unos*unos'
	R = NBM - unnos
	R
end


function eigens(M,threshold)
	valores, vectores = eigs(M,which=:LR,nev=3)
    ind = findall(x->imag(x)==0.0,valores)
    valores = valores[ind]
    vectores = vectores[:,ind]
    if real(last(valores)) > threshold
		valores, vectores = eigs(M,which=:LR,nev=4)
        ind = findall(x->imag(x)==0.0,valores)
        valores = valores[ind]
        vectores = vectores[:,ind]
		if real(last(valores)) > threshold
		    valores, vectores = eigs(M,which=:LR,nev=5)
            ind = findall(x->imag(x)==0.0,valores)
            valores = valores[ind]
            vectores = vectores[:,ind]
		    if real(last(valores)) > threshold
		        valores, vectores = eigs(M,which=:LR,nev=10)
                ind = findall(x->imag(x)==0.0,valores)
                valores = valores[ind]
                vectores = vectores[:,ind]
		        if real(last(valores)) > threshold
		            valores, vectores = eigs(M,which=:LR,nev=nv(g))
                    ind = findall(x->imag(x)==0.0,valores)
                    valores = valores[ind]
                    vectores = vectores[:,ind]
		        end
		    end
		end
	end
	valores, vectores
end



function num_com(v,threshold)
	cuantos = Array{Float64}(undef,0)
	index = Array{Int64}(undef,0)
	for i in 1:length(v)
		if imag(v[i]) == 0 && real(v[i]) > threshold
		    push!(cuantos,real(v[i]))
		    push!(index,i)
		end
	end
    cuantos,index
end


function embe_nbm(g,ind,vec)
	reales = real(vectores[:,index])
	arriba = reales[1:nv(g),:]
	abajo = reales[nv(g)+1:2*nv(g),:]
	matriz_embedded = arriba + abajo
	matriz_embedded
end

function contract(g,v,ei)
    y = zeros(Float64, nv(g))
    for i in 1:nv(g)
        for j in inneighbors(g,i)
            u = ei[Edge(j,i)]
            y[i] += v[u]
        end
    end
    y
end

function contraction(g,vec,emb,ei)
	contracted = zeros(Float64, nv(g), length(vec))
	for n in 1:length(vec)
    	contracted[:,n] = contract(g,emb,ei)
	end
	contracted
end
