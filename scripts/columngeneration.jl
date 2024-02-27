
using Base.Iterators: partition

#--------------------------------------------------------------------------------------------#

function orderpathinitialization(c_mag, mcfinstance, numarcs_dummy)

	delta = Dict()
	pathcost = Dict()
	pathSet = Dict()
	for k in mcfinstance.commodities
		pathSet[k] = []
		pathcost[k] = []
	end

	#Generate initial dummy paths
	arcindex = mcfinstance.numarcs + 1
	for k in mcfinstance.commodities
		p = 1
		push!(pathSet[k], p)
		for a in 1:numarcs_dummy
			if a == arcindex
				delta[k,a,p] = 1
			else
				delta[k,a,p] = 0
			end
		end
		push!(pathcost[k], c_mag[arcindex, k])
	    arcindex += 1
	end

	return delta, pathcost, pathSet

end

#--------------------------------------------------------------------------------------------#

function columngeneration!(pathSet, pathcost, delta, mcfinstance, numarcs_dummy)

	#Initialize
	cg_iteration = 1
	rmp_time, pp_time, pp_time_par = 0, 0, 0
	parallel_time = 0
	listlength = convert(Int64, ceil(length(mcfinstance.commodities)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)
	M = preprocessreducedcostsets(mcfinstance)

	#Create the sparse master problem
	m = Model(Gurobi.Optimizer)
	set_optimizer_attribute(m, "OutputFlag", 0)
	@variable(m, 0 <= y[k in mcfinstance.commodities, p in pathSet[k]]) #Upper bound redundant
	@objective(m, Min, sum(sum(pathcost[k][p] * mcfinstance.q[k] * y[k,p] for p in pathSet[k]) for k in mcfinstance.commodities) )
	@constraint(m, assignment[k in 1:numcom], sum(y[k,p] for p in pathSet[k]) == 1)
	@constraint(m, arccapacity[a = 1:mcfinstance.numarcs], sum(sum(mcfinstance.q[k] * delta[k,a,p] * y[k,p] for p in pathSet[k]) for k in mcfinstance.commodities) <= mcfinstance.d[a])
	@constraint(m, nodecapacity[n in 1:numnodes], sum(sum(sum(mcfinstance.q[k] * delta[k,a,p] * y[k,p] for a in union(mcfinstance.A_plus[n], mcfinstance.A_minus[n])) for p in pathSet[k]) for k in mcfinstance.commodities) <= mcfinstance.p[n])

	fullalgstarttime = time()

	#Main loop
	while 1==1

		#Solve RMP
		optimize!(m)
		rmp_time += solve_time(m)
		rmpobj = objective_value(m)
		println("CG iteration ", cg_iteration, " RMP objective = ", rmpobj)

		#Snapshot of current arcs
		pathSet_snapshot = deepcopy(pathSet)

		#Get dual values
		alpha = dual.(arccapacity)
		beta = dual.(assignment)
		iota = dual.(nodecapacity)

		#Find arc reduced costs
		#arcredcosts = Dict()
		#for k in mcfinstance.commodities, a in 1:mcfinstance.numarcs
		#	i, j = mcfinstance.arcs[a]
		#	arcredcosts[k,a] = mcfinstance.c[a, k] * mcfinstance.q[k] - mcfinstance.q[k] * alpha[a] - mcfinstance.q[k] * iota[i] - mcfinstance.q[k] * iota[j]
		#end
        arcredcosts = transpose(mcfinstance.c[:, 1] * transpose(mcfinstance.q) - alpha * transpose(mcfinstance.q) + M.iota * iota * transpose(mcfinstance.q))

		#Solve pricing subproblem
		min_rc_list = []
		dptimelist = []
		pathsadded = 0
		for k in mcfinstance.commodities
			thisiterstarttime = time()
			minreducedcost_k, shortestpathnodes, shortestpatharcs = findshortestpath_mcf(0, k, arcredcosts, numnodes, mcfinstance.numarcs, mcfinstance.arcs, mcfinstance.arcLookup, mcfinstance.Origin, mcfinstance.Destination)
			push!(dptimelist, time() - thisiterstarttime)

			#Complete the reduced cost calculation 
			fullminreducedcost = minreducedcost_k - beta[k]
			push!(min_rc_list, fullminreducedcost)

			#Add new arcs to current arc sets
			if fullminreducedcost < -0.001

				#Define new path
				p = last(pathSet[k]) + 1
				push!(pathSet[k], p)
				for a in shortestpatharcs
					delta[k,a,p] = 1
				end
				for a in setdiff(1:numarcs_dummy, shortestpatharcs)
					delta[k,a,p] = 0
				end
				newpathcost = sum(mcfinstance.c[a,k] for a in shortestpatharcs)
				push!(pathcost[k], newpathcost)

				#Add new path variable
				global y[k,p] = @variable(m, lower_bound = 0) #, upper_bound = 1)
				set_name(y[k,p], string("y[",k,",",p,"]")) 
				#set_objective_function(m, objective_function(m) + pathcost[k][p] * mcfinstance.q[k] * y[k,p])
                set_objective_coefficient(m, y[k,p], pathcost[k][p] * mcfinstance.q[k])
				for a in shortestpatharcs
					set_normalized_coefficient(arccapacity[a], y[k,p], mcfinstance.q[k])
				end
				for n in mcfinstance.nodes
					totalcoeff = 0
					for a in intersect(shortestpatharcs, union(mcfinstance.A_plus[n], mcfinstance.A_minus[n]))
						totalcoeff += mcfinstance.q[k]
					end
					set_normalized_coefficient(nodecapacity[n], y[k,p], totalcoeff)
				end
				set_normalized_coefficient(assignment[k], y[k,p], 1.0)

				pathsadded += 1
				
			end
		end

		println("Paths added = ", pathsadded)

	 	min_rc = minimum(min_rc_list)

		#"Parallelize"
		shuffleddptimes = shuffle_partition(length(mcfinstance.commodities))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		pp_time_par += maximum(dptimelistsums)
		parallel_time += sum(dptimelistsums) -  maximum(dptimelistsums)

		#Termination criterion
		if min_rc >= -0.0001
			println("No negative reduced costs found!")
			break
		end

		cg_iteration += 1

	end

	fullalgtime = time() - fullalgstarttime - parallel_time

	return objective_value(m), cg_iteration, pathSet, pathcost, delta, rmp_time, pp_time_par, fullalgtime, y

end