
function getbasisarcs(mcfinstance, x_lp)

    commArcSet = [[] for k in mcfinstance.commodities]
	A_minus_k, A_plus_k = Dict(), Dict()
	for i in mcfinstance.nodes, k in mcfinstance.commodities
		A_minus_k[k, i] = []
		A_plus_k[k, i] = []
	end

	#Find basis
	for k in mcfinstance.commodities, a in 1:mcfinstance.numarcs
        if value(x_lp[k,a]) > 1e-5
            i, j = mcfinstance.arcs[a]
            push!(A_minus_k[k,j], a)
            push!(A_plus_k[k,i], a)
            push!(commArcSet[k], a)
        end
	end

    basisarcs = (A=commArcSet, A_minus=A_minus_k, A_plus=A_plus_k)
    
    return basisarcs

end
