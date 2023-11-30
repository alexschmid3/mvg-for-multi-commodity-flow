
function createinstance_mcf(radius, destdistpercentile, maxdistanceperturb, mindistanceperturb)

	commodities = 1:numcom
	nodes = 1:numnodes

	#Randomize nodes
	coordinates = zeros(numnodes, 2)
	for i in nodes
		coordinates[i,2] = rand()
		coordinates[i,1] = rand()
	end

	#Find possible destinations for each origin node
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
	c = Dict()
	for i in nodes
		A_minus[i] = []
		A_plus[i] = []
	end
	for i in nodes, j in nodes
		if i != j
			dist = norm(coordinates[i,:] - coordinates[j,:],2)
			if dist <= radius
				arcs[arcindex] = (i, j)
				arcLookup[i, j] = arcindex
				push!(A_minus[j], arcindex)
				push!(A_plus[i], arcindex)
				dist = norm(coordinates[i,:] - coordinates[j,:],2)
				perturbation = rand() * (maxdistanceperturb - mindistanceperturb) + mindistanceperturb
				arccost = perturbation * dist
		 		for k in commodities
					c[arcindex, k] = round(arccost, digits = 2)
				end
				arcindex += 1
			end
		end
	end
	numarcs = arcindex - 1

	return coordinates, commodities, nodes, arcs, arcLookup, numarcs, A_minus, A_plus, c, b, q, Origin, Destination

end

function randomizecapacities_mcf(numarcs, nodes, maxcapacityperturb, mincapacityperturb)

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

function setcapacities_mcf(gamma_arc, gamma_node, q, numarcs, nodes, arcperturbation, nodeperturbation)

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

function findshortestpath_mcf(k, arcredcosts, numnodes, numarcs, arcs, arcLookup, Origin, Destination)

	#Find origin and destination nodes 
	orig = Origin[k]
	dest = Destination[k]

	#Initialize shortest path algorithm (Bellman-Ford)
	currdistance = repeat([999999999.0],outer=[numnodes])
	currdistance[orig] = 0
	prevnode, prevarc = zeros(numnodes), zeros(numnodes)
	
	for iteration in 1:numnodes #Max path length
		for a in 1:numarcs
			n_end, n_start = arcs[a][2], arcs[a][1]
			if currdistance[n_end] > currdistance[n_start] + arcredcosts[k,a] + 0.00001
				currdistance[n_end] = currdistance[n_start] + arcredcosts[k,a]
				prevnode[n_end] = n_start
				prevarc[n_end] = a
			end
		end
	end

	#Format the shortest path output
	shortestpathnodes_rev = [dest]
	shortestpatharcs_rev = []
	node = dest
	while node != orig
		push!(shortestpatharcs_rev, Int(prevarc[node]))
		node = Int(prevnode[node])
		push!(shortestpathnodes_rev, node)
	end
	shortestpathnodes = reverse(shortestpathnodes_rev) 
	shortestpatharcs = reverse(shortestpatharcs_rev) 

	return currdistance[dest], shortestpathnodes, shortestpatharcs

end

function findgoodinstance_arctuning(gamma_arc_init, gamma_node, optgap, maxiter, timegoal, numcom, numnodes, maxorder, minorder, radius, destdistpercentile, maxdistanceperturb, mindistanceperturb, maxcapacityperturb, mincapacityperturb)
	
	currgamma = gamma_arc_init
	bestfeasiblegamma = 0
	bestinfeasiblegamma = 50*gamma_arc_init
	goodinstance_flag = 0
	iter = 1

	#Create instance
	coordinates, commodities, nodes, arcs, arcLookup, numarcs, A_minus, A_plus, c, b, q, Origin, Destination = createinstance_mcf(radius, destdistpercentile, maxdistanceperturb, mindistanceperturb)
	arcperturbation, nodeperturbation = randomizecapacities_mcf(numarcs, nodes, maxcapacityperturb, mincapacityperturb)

	while (goodinstance_flag == 0) & (iter <= maxiter)
		println("Iteration $iter arc gamma = ", currgamma)
		
		#Randomize capacities
		d, p = setcapacities_mcf(currgamma, gamma_node, q, numarcs, nodes, arcperturbation, nodeperturbation)

		#Write corresponding MCF IP
		m = Model(Gurobi.Optimizer)
		set_optimizer_attribute(m, "TimeLimit", timegoal+1)
		set_optimizer_attribute(m, "OutputFlag", 0)
		set_optimizer_attribute(m, "MIPGap", optgap)
		@variable(m, x_ip[k in commodities, a in 1:numarcs], Bin)
		@objective(m, Min, sum(sum(c[a,k] * q[k] * x_ip[k,a] for a in 1:numarcs) for k in commodities) )
		@constraint(m, arccapacitycon_ip[a = 1:numarcs], sum(q[k] * x_ip[k, a] for k in commodities) <= d[a])
		#@constraint(m, nodecapacitycon_ip[n in nodes], sum(sum(q[k] * x_ip[k, a] for a in union(A_plus[n], A_minus[n])) for k in commodities) <= p[n])
		@constraint(m, flowbalance_ip[k in commodities, i in nodes], sum(x_ip[k, a] for a in A_plus[i]) - sum(x_ip[k, a] for a in A_minus[i]) == b[i,k])

		#Solve IP
		status_ip = optimize!(m)

		#Check status
		if (termination_status(m) == MOI.OPTIMAL) & (solve_time(m) < timegoal)
			total_time = solve_time(m)
			println("Too short = ", total_time)
			bestfeasiblegamma = currgamma
			currgamma = round((currgamma + bestinfeasiblegamma) / 2, digits=0)
		elseif (termination_status(m) == MOI.OPTIMAL) & (solve_time(m) >= timegoal)
			total_time = solve_time(m)
			bestfeasiblegamma = currgamma
			println("Good instance = ", total_time)
			goodinstance_flag = 1
		elseif termination_status(m) == MOI.TIME_LIMIT
			total_time = solve_time(m)
			bestfeasiblegamma = currgamma
			println("Good instance = ", total_time)
			goodinstance_flag = 1
		elseif (termination_status(m) == MOI.INFEASIBLE_OR_UNBOUNDED) || (termination_status(m) == MOI.INFEASIBLE)
			println("Infeasible")
			total_time = "Inf"
			bestinfeasiblegamma = currgamma
			currgamma = round((currgamma + bestfeasiblegamma) / 2, digits=0)
		end

		df = DataFrame(experiment_id = [runid],
			randomseed = [randomseedval],
			numcom = [numcom], 
			numnodes = [numnodes], 
			radius = [radius], 
			gamma_init = [gamma_arc_init], 
			gamma_arc = [currgamma], 
			gamma_node = [gamma_node],
			optgap = [optgap], 
			numarcs = [numarcs], 
			total_time = [total_time],
			goodinstance_flag = [goodinstance_flag]
           )

		if iter == 1
			CSV.write(outputfilename, df)
		else
			CSV.write(outputfilename, df, append=true)
		end

		iter += 1

	end

	return bestfeasiblegamma, goodinstance_flag

end

function findgoodinstance_nodetuning(gamma_arc, gamma_node_init, optgap, maxiter, timegoal, numcom, numnodes, maxorder, minorder, radius, destdistpercentile, maxdistanceperturb, mindistanceperturb, maxcapacityperturb, mincapacityperturb)
	
	#gamma_arc, gamma_node_init, optgap, maxiter, timegoal, numcom, numnodes, maxorder, minorder, radius, destdistpercentile, maxdistanceperturb, mindistanceperturb, maxcapacityperturb, mincapacityperturb = gamma_arc, gamma_node_init, opt_gap, maxiter, timegoal2, numcom, numnodes, maxorder, minorder, radius, destdistpercentile, maxdistanceperturb, mindistanceperturb, maxcapacityperturb, mincapacityperturb

	currgamma = gamma_node_init
	bestfeasiblegamma = 0
	bestinfeasiblegamma = max(4*gamma_node_init, 10000)
	goodinstance_flag = 0
	iter = 1

	#Create instance
	coordinates, commodities, nodes, arcs, arcLookup, numarcs, A_minus, A_plus, c, b, q, Origin, Destination = createinstance_mcf(radius, destdistpercentile, maxdistanceperturb, mindistanceperturb)
	arcperturbation, nodeperturbation = randomizecapacities_mcf(numarcs, nodes, maxcapacityperturb, mincapacityperturb)

	while (goodinstance_flag == 0) & (iter <= maxiter)
		println("Iteration $iter node gamma = ", currgamma)
		
		#Randomize capacities
		d, p = setcapacities_mcf(gamma_arc, currgamma, q, numarcs, nodes, arcperturbation, nodeperturbation)

		#Write corresponding MCF IP
		m = Model(Gurobi.Optimizer)
		set_optimizer_attribute(m, "TimeLimit", timegoal+1)
		set_optimizer_attribute(m, "OutputFlag", 0)
		set_optimizer_attribute(m, "MIPGap", optgap)
		@variable(m, x_ip[k in commodities, a in 1:numarcs], Bin)
		@objective(m, Min, sum(sum(c[a,k] * q[k] * x_ip[k,a] for a in 1:numarcs) for k in commodities) )
		@constraint(m, arccapacitycon_ip[a = 1:numarcs], sum(q[k] * x_ip[k, a] for k in commodities) <= d[a])
		@constraint(m, nodecapacitycon_ip[n in nodes], sum(sum(q[k] * x_ip[k, a] for a in union(A_plus[n], A_minus[n])) for k in commodities) <= p[n])
		@constraint(m, flowbalance_ip[k in commodities, i in nodes], sum(x_ip[k, a] for a in A_plus[i]) - sum(x_ip[k, a] for a in A_minus[i]) == b[i,k])
			
		#Solve IP
		status_ip = optimize!(m)

		#Check status
		if (termination_status(m) == MOI.OPTIMAL) & (solve_time(m) < timegoal)
			total_time = solve_time(m)
			println("Too short = ", total_time)
			bestfeasiblegamma = currgamma
			currgamma = round((currgamma + bestinfeasiblegamma) / 2, digits=0)
		elseif (termination_status(m) == MOI.OPTIMAL) & (solve_time(m) >= timegoal)
			total_time = solve_time(m)
			bestfeasiblegamma = currgamma
			println("Good instance = ", total_time)
			goodinstance_flag = 1
		elseif termination_status(m) == MOI.TIME_LIMIT
			total_time = solve_time(m)
			bestfeasiblegamma = currgamma
			println("Good instance = ", total_time)
			goodinstance_flag = 1
		elseif (termination_status(m) == MOI.INFEASIBLE_OR_UNBOUNDED) || (termination_status(m) == MOI.INFEASIBLE)
			println("Infeasible")
			total_time = "Inf"
			bestinfeasiblegamma = currgamma
			currgamma = round((currgamma + bestfeasiblegamma) / 2, digits=0)
		end

		df = DataFrame(experiment_id = [runid],
			randomseed = [randomseedval],
			numcom = [numcom], 
			numnodes = [numnodes], 
			radius = [radius], 
			gamma_init = [gamma_node_init], 
			gamma_arc = [gamma_arc], 
			gamma_node = [currgamma],
			optgap = [optgap], 
			numarcs = [numarcs], 
			total_time = [total_time],
			goodinstance_flag = [goodinstance_flag]
           )

		if iter == 1
			CSV.write(outputfilename, df)
		else
			CSV.write(outputfilename, df, append=true)
		end

		iter += 1

	end

	return currgamma, goodinstance_flag

end