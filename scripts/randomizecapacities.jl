
function randomizeperturbations(numarcs, nodes, randomseedval)

	Random.seed!(randomseedval)

	#Arc and node capacity randomization
	arcperturbation, nodeperturbation = Dict(), Dict()
	for a in 1:numarcs
		arcperturbation[a] = rand() * (maxcapacityperturb - mincapacityperturb) + mincapacityperturb
	end
	for n in nodes
		nodeperturbation[n] = rand() * (maxcapacityperturb - mincapacityperturb) + mincapacityperturb
	end

	return arcperturbation, nodeperturbation

end

#--------------------------------------------------------------------------------------------#

function setcapacities(gamma_arc, gamma_node, q, numarcs, nodes, arcperturbation, nodeperturbation)
    
	#Arc and node capacities
	d, p = [], []
	avgcapacity = sum(q) * numarcs / gamma_arc
	for a in 1:numarcs
		arccapacity = max(1,round(avgcapacity * arcperturbation[a]))
		push!(d, arccapacity)
	end
	avgnodecapacity = sum(q) * numnodes / gamma_node
	for n in nodes
		nodecapacity = max(1,round(avgnodecapacity * nodeperturbation[n]))
		push!(p, nodecapacity)
	end

    return d, p

end

#--------------------------------------------------------------------------------------------#

function randomizecapacities(gamma_arc, gamma_node, q, numarcs, nodes, randomseedval, maxcapacityperturb, mincapacityperturb)

    arcperturbation, nodeperturbation = randomizeperturbations(numarcs, nodes, maxcapacityperturb, mincapacityperturb, randomseedval)
	d,p = setcapacities(gamma_arc, gamma_node, q, numarcs, nodes, arcperturbation, nodeperturbation)

	return d, p

end

