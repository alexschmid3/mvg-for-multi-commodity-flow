
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

#--------------------------------------------------------------------------------------------#

function findpossibledestinations(coordinates, nodes)

    possibleDestinations = Dict()
	numdestinations = convert(Int64, ceil(destdistpercentile*(numnodes - 1)))
	for i in nodes
		destlist = []
		distancedict = Dict()
		for j in nodes
			if i != j
				dist = norm(coordinates[i,:] - coordinates[j,:],2)
				distancedict[j] = dist
			end
		end
		sorteddict = sort(collect(distancedict), by=x->x[2])[numnodes - numdestinations : numnodes - 1]
		for n in sorteddict
			push!(destlist, n[1])
		end
		possibleDestinations[i] = destlist
	end

    return possibleDestinations, numdestinations

end

#--------------------------------------------------------------------------------------------#

function checktooclose(tempcoords, coordinates)

	for j in 1:size(coordinates)[1]
		if sqrt((tempcoords[1] - coordinates[j,1])^2 + (tempcoords[2] - coordinates[j,2])^2) <= 0.01
			return true
		end
	end
	return false
end

#--------------------------------------------------------------------------------------------#

function randomizeinstance(randomseedval)

	Random.seed!(randomseedval)
	
	commodities = 1:numcom
	nodes = 1:numnodes

	#Randomize coordinates of nodes
	coordinates = zeros(numnodes, 2)
	for i in nodes
		while 1==1
			tempcoords = rand(), rand()
			if !(checktooclose(tempcoords, coordinates))
				coordinates[i,1] = tempcoords[1]
				coordinates[i,2] = tempcoords[2]
				break
			end
		end
	end

	#Find possible destinations for each origin node
	possibleDestinations, numdestinations = findpossibledestinations(coordinates, nodes)

	#Randomize commodities 
	b, q = zeros(numnodes, numcom), []
	Origin, Destination = [], []
	for k in commodities
		orig = randperm(numnodes)[1]
		dest = possibleDestinations[orig][randperm(numdestinations)[1]]
		b[orig, k] = 1
		b[dest, k] = -1
		orderqty = rand(minorder:maxorder)
		#orderqty = round(rand() * (maxorder - minorder) + minorder)
		push!(q, orderqty)
		push!(Origin, orig)
		push!(Destination, dest)
	end

	#Create arcs between all nodes + arc costs
	arcindex = 1
	arcs, arcLookup = Dict(), Dict()
	A_minus, A_plus = Dict(), Dict()
	for i in nodes
		A_minus[i] = []
		A_plus[i] = []
	end
	arcdist = Dict()
	for i in nodes, j in nodes
		if i != j
			dist = norm(coordinates[i,:] - coordinates[j,:],2)
			if dist <= radius
				arcs[arcindex] = (i, j)
				arcLookup[i, j] = arcindex
				push!(A_minus[j], arcindex)
				push!(A_plus[i], arcindex)
				arcdist[arcindex] = dist
				arcindex += 1
			end
		end
	end
	numarcs = arcindex - 1

	c = zeros(numarcs, length(commodities))
	for a in 1:numarcs
		dist = arcdist[a]
		perturbation = rand() * (maxdistanceperturb - mindistanceperturb) + mindistanceperturb
		arccost = perturbation * dist
		for k in commodities
			c[a, k] = round(arccost, digits = 2)
		end
	end

	return coordinates, commodities, nodes, arcs, arcLookup, numarcs, A_minus, A_plus, c, b, q, Origin, Destination

end

