
using Base.Iterators: partition

#--------------------------------------------------------------------------------------------#

function colgeninitialize(mcfinstance)
	
	c_mag, d_mag = copy(mcfinstance.c), copy(mcfinstance.d)
	commArcSet = [[] for k in mcfinstance.commodities]
	A_minus_k, A_plus_k = Dict(), Dict()
	for i in mcfinstance.nodes, k in mcfinstance.commodities
		A_minus_k[k, i] = []
		A_plus_k[k, i] = []
	end

	arcindex = mcfinstance.numarcs + 1
	dummydeletions = []
	#Dummy arcs to initialize column generation
	for k in mcfinstance.commodities
		push!(A_minus_k[k, mcfinstance.Destination[k]], arcindex)
		push!(A_plus_k[k, mcfinstance.Origin[k]], arcindex)
		#c_mag[arcindex, k] = 100000
		push!(commArcSet[k], arcindex)
		push!(dummydeletions, (k, arcindex, mcfinstance.Origin[k], mcfinstance.Destination[k]))
		arcindex += 1
	end
	numarcs_dummy = arcindex - 1

	for a in mcfinstance.numarcs+1:numarcs_dummy
		c_mag = vcat(c_mag, transpose([100000 for k in mcfinstance.commodities]))
	end

	#Adjust capacities
	for a in mcfinstance.numarcs+1:numarcs_dummy
		push!(d_mag, sum(mcfinstance.q))
	end

    magarcs = (A=commArcSet, A_minus=A_minus_k, A_plus=A_plus_k)

	return c_mag, d_mag, magarcs, numarcs_dummy, dummydeletions
	
end

#--------------------------------------------------------------------------------------------#

function findshortestpath_mcf(exactflag, k, arcredcosts, numnodes, numarcs, arcs, arcLookup, Origin, Destination)

	#Find origin and destination nodes 
	orig = Origin[k]
	dest = Destination[k]

	#Initialize shortest path algorithm (Bellman-Ford)
	currdistance = repeat([999999999.0],outer=[numnodes])
	currdistance[orig] = 0
	prevnode, prevarc = zeros(numnodes), zeros(numnodes)
	
	if exactflag == 0
		maxspiter = 10
	elseif exactflag == 1
		maxspiter = min(50, numnodes)
	end

	for iteration in 1:maxspiter #numnodes #Max path length
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

#--------------------------------------------------------------------------------------------#

function multiarcgeneration!(mcfinstance, magarcs, c_mag, d_mag, numarcs_dummy, dummydeletions)

	#Initialize
	mag_iteration = 1
	smp_time, magsp_time, magsp_time_par = 0, 0, 0
	parallel_time = 0
	listlength = convert(Int64, ceil(length(commodities)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)

	#Create the sparse master problem
	m = Model(Gurobi.Optimizer)
	set_optimizer_attribute(m, "OutputFlag", 0)
	@variable(m, 0 <= x[k in commodities, a in magarcs.A[k]]) #Upper bound redundant
	@objective(m, Min, sum(sum(c_mag[a,k] * mcfinstance.q[k] * x[k,a] for a in magarcs.A[k]) for k in mcfinstance.commodities) )
	@constraint(m, arccapacity[a = 1:numarcs_dummy], 0 <= d_mag[a])
	@constraint(m, flowbalance[k in mcfinstance.commodities, i in mcfinstance.nodes], sum(x[k, a] for a in magarcs.A_plus[k, i]) - sum(x[k, a] for a in magarcs.A_minus[k, i]) == mcfinstance.b[i,k])
	@constraint(m, nodecapacity[n in mcfinstance.nodes], 0 <= mcfinstance.p[n])

	#Main loop
    fullalgstarttime = time()

	time1, time2, time3, time4, time5 = 0,0,0,0,0

	while 1==1

		tempstart = time()
		#Solve SMP
		optimize!(m)
		smp_time += solve_time(m)
		rmpobj = objective_value(m)
		println("MAG iteration ", mag_iteration, " SMP objective = ", rmpobj)
		time1 += time() - tempstart

		#Get dual values
		tempstart = time()
		alpha1 = dual.(arccapacity)
		beta = dual.(flowbalance)
		iota = dual.(nodecapacity)
		time2 += time() - tempstart

		#Find arc reduced costs
		tempstart = time()
		arcredcosts = Dict()
		for k in mcfinstance.commodities, a in 1:mcfinstance.numarcs
			i, j = mcfinstance.arcs[a]
			arcredcosts[k,a] = c_mag[a, k] * mcfinstance.q[k] - mcfinstance.q[k] * alpha[a] - beta[k, i] + beta[k, j] - mcfinstance.q[k] * iota[i] - mcfinstance.q[k] * iota[j]
		end
		#=
		#Optimized MCF
		arcredcosts = zeros()
		for a in 1:mcfinstance.numarcs
			arcredcosts[:,a] += c_mag[:, 1] * transpose(mcfinstance.q)
		end
		arcredcosts += c_mag[:, 1] * transpose(mcfinstance.q) - alpha1 * transpose(mcfinstance.q)

		#Optimized relay
		arcredcosts = zeros(numorders, extendednumarcs)
		for i in orders
			arcredcosts[i,:] += c[1:extendednumarcs] + M.alpha[i] * alpha[i,:] + M.beta[i] * beta[i] + M.gamma[i] * gamma[i] + M.psi[i] * psi[i]
			arcredcosts[i,:] += M.theta * theta + M.nu * nu + M.mu * mu + M.xi * xi
		end=#
		time3 += time() - tempstart

		#Solve MAG subproblem
		min_rc_list = []
		addarcs = []
		dptimelist = []
		for k in mcfinstance.commodities
			thisiterstarttime = time()
			minreducedcost_k, shortestpathnodes, shortestpatharcs = findshortestpath_mcf(0, k, arcredcosts, numnodes, mcfinstance.numarcs, mcfinstance.arcs, mcfinstance.arcLookup, mcfinstance.Origin, mcfinstance.Destination)
			push!(min_rc_list, minreducedcost_k)

			#Add new arcs to current arc sets
			if minreducedcost_k < -0.001
				for a in shortestpatharcs
					if !(a in magarcs.A[k])
                        push!(addarcs, (k,a))
						push!(magarcs.A[k], a)
						i, j = mcfinstance.arcs[a]
						push!(magarcs.A_plus[k, i], a)
						push!(magarcs.A_minus[k, j], a)
					end
				end
			end

			push!(dptimelist, time() - thisiterstarttime)
		end

	 	min_rc = minimum(min_rc_list)

		#"Parallelize"
		shuffleddptimes = shuffle_partition(length(mcfinstance.commodities))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		magsp_time_par += maximum(dptimelistsums)
		parallel_time += sum(dptimelistsums) -  maximum(dptimelistsums)

		#Add new variables 
		tempstart = time()
		for (k,a) in addarcs
            #Create a new variable for the arc
            global x[k,a] = @variable(m, lower_bound = 0) #, upper_bound = 1)
            set_name(x[k,a], string("x[",k,",",a,"]")) 

            #Add to the objective
            set_objective_function(m, objective_function(m) + c_mag[a,k] * mcfinstance.q[k] * x[k,a])

            #Add new variable to order constraints
            i = mcfinstance.arcs[a][1]
            j = mcfinstance.arcs[a][2]
            set_normalized_coefficient(arccapacity[a], x[k,a], mcfinstance.q[k])
            set_normalized_coefficient(flowbalance[k,i], x[k,a], 1.0)
            set_normalized_coefficient(flowbalance[k,j], x[k,a], -1.0)
            set_normalized_coefficient(nodecapacity[i], x[k,a], mcfinstance.q[k])
            set_normalized_coefficient(nodecapacity[j], x[k,a], mcfinstance.q[k])
		end
		time4 += time() - tempstart

		#Termination criterion
		tempstart = time()
		if min_rc >= -0.001
			println("No negative reduced costs found!")
			break
		end
		
		mag_iteration += 1
		time5 += time() - tempstart

	end

    fullalgtime = time() - fullalgstarttime - parallel_time
	
	#Remove dummy arcs
	for (k, a, i, j) in dummydeletions
		remove!(magarcs.A[k], a)
        remove!(magarcs.A_plus[k, i], a)
        remove!(magarcs.A_minus[k, j], a)
	end

	println("Time 1 = ", time1)
	println("Time 2 = ", time2)
	println("Time 3 = ", time3)
	println("Time 4 = ", time4)
	println("Time 5 = ", time5)
	println("Full time = $fullalgtime vs. ", time1+time2+time3+time4+time5+magsp_time_par)

	return objective_value(m), mag_iteration, magarcs, smp_time, magsp_time_par, fullalgtime, x

end

#--------------------------------------------------------------------------------------------#
