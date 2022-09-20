
using Base.Iterators: partition

#----------------------------------------------------------------------------------------------------------#

function solvemcfmodel(lprelaxation_flag, commArcSet, A_plus_k, A_minus_k, c_ip)

	#lprelaxation_flag, commArcSet, A_plus_k, A_minus_k, c_ip = 1, commArcSet_converged, A_plus_k_converged, A_minus_k_converged, c_cg

	m = Model(Gurobi.Optimizer)
	set_optimizer_attribute(m, "TimeLimit", 60*60*95)
	set_optimizer_attribute(m, "OutputFlag", 0)
	set_optimizer_attribute(m, "MIPGap", 0.01)

	if lprelaxation_flag == 1
		@variable(m, 0 <= x[k in commodities, a in commArcSet[k]] <= 1)
	elseif lprelaxation_flag == 0
		@variable(m, x[k in commodities, a in commArcSet[k]], Bin)
	end

	@objective(m, Min, sum(sum(c_cg[a,k] * q[k] * x[k,a] for a in commArcSet[k]) for k in commodities) )

	@constraint(m, arccapacity[a in 1:numarcs], sum(q[k] * x[k, a] for k in commodities if a in commArcSet[k]) <= d[a])
	@constraint(m, flowbalance[k in commodities, n in nodes], sum(x[k, a] for a in A_plus_k[k, n]) - sum(x[k, a] for a in A_minus_k[k, n]) == b[n,k])
	@constraint(m, nodecapacity[n in nodes], sum(sum(q[k] * x[k, a] for a in setdiff(union(A_plus_k[k, n], A_minus_k[k, n]), numarcs+1:numarcs_dummy)) for k in commodities) <= qnode[n])

	status_ip = optimize!(m)

	println("Objective = ", getobjectivevalue(m), " in ", solve_time(m), " seconds")

	return getobjectivevalue(m), solve_time(m)

end

#----------------------------------------------------------------------------------------------------------#

function colgeninitialize(c, d, numarcs)
	
	c_cg, d_cg = copy(c), copy(d)
	commArcSet = [[] for k in commodities]
	A_minus_k, A_plus_k = Dict(), Dict()
	for i in nodes, k in commodities
		A_minus_k[k, i] = []
		A_plus_k[k, i] = []
	end

	arcindex = numarcs + 1
	#Dummy arcs to initialize column generation
	for k in commodities
		push!(A_minus_k[k, Destination[k]], arcindex)
		push!(A_plus_k[k, Origin[k]], arcindex)
		c_cg[arcindex, k] = 100000
		push!(commArcSet[k], arcindex)
		arcindex += 1
	end
	numarcs_dummy = arcindex - 1

	#Adjust capacities
	for a in numarcs+1:numarcs_dummy
		push!(d_cg, sum(q))
	end

	return c_cg, d_cg, commArcSet, A_minus_k, A_plus_k, numarcs_dummy
	
end

#----------------------------------------------------------------------------------------------------------#

function sparsemasterproblem(commArcSet, A_plus_k, A_minus_k, c_cg)

	m = Model(Gurobi.Optimizer)
	set_optimizer_attribute(m, "OutputFlag", 0)

	@variable(m, 0 <= x[k in commodities, a in commArcSet[k]]) #Upper bound redundant

	@objective(m, Min, sum(sum(c_cg[a,k] * q[k] * x[k,a] for a in commArcSet[k]) for k in commodities) )

	@constraint(m, arccapacity[a = 1:numarcs_dummy], 0 <= d_cg[a])
	@constraint(m, flowbalance[k in commodities, i in nodes], sum(x[k, a] for a in A_plus_k[k, i]) - sum(x[k, a] for a in A_minus_k[k, i]) == b[i,k])
	@constraint(m, nodecapacity[n in nodes], 0 <= qnode[n])

	smpconstraints = (con_arccapacity = arccapacity,
		con_flowbalance = flowbalance,
		con_nodecapacity = nodecapacity
		)

	return m, x, smpconstraints

end

#----------------------------------------------------------------------------------------------------------#

function multivariablegeneration!(commArcSet, A_plus_k, A_minus_k, onearcatatime_flag, c_cg)

	#Initialize
	mvg_iteration = 1
	smp_time, mvgsp_time, mvgsp_time_par = 0, 0, 0
	parallel_time = 0
	listlength = convert(Int64, ceil(length(commodities)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)

	#Create the sparse master problem
	m, x, smpconstraints = sparsemasterproblem(commArcSet, A_plus_k, A_minus_k, c_cg) 

	fullalgstarttime = time()

	#Main loop
	while 1==1

		#Solve RMP
		optimize!(m)
		smp_time += solve_time(m)
		rmpobj = getobjectivevalue(m)
		println("MVG iteration ", mvg_iteration, " SMP objective = ", rmpobj)

		#Snapshot of current arcs
		commArcSet_snapshot = deepcopy(commArcSet)

		#Get dual values
		alpha = getdual.(smpconstraints.con_arccapacity)
		beta = getdual.(smpconstraints.con_flowbalance)
		iota = getdual.(smpconstraints.con_nodecapacity)

		#Find arc reduced costs
		arcredcosts = Dict()
		for k in commodities, a in 1:numarcs
			i, j = arcs[a]
			arcredcosts[k,a] = c_cg[a, k] * q[k] - q[k] * alpha[a] - beta[k, i] + beta[k, j] - q[k] * iota[i] - q[k] * iota[j]
		end

		#Solve MVG subproblem
		min_rc_list = []
		addarcs = []
		dptimelist = []
		for k in commodities
			if onearcatatime_flag == 0
				thisiterstarttime = time()
				minreducedcost_k, shortestpathnodes, shortestpatharcs = findshortestpath_mcf(k, arcredcosts, numnodes, numarcs, arcs, arcLookup, Origin, Destination)
				push!(dptimelist, time() - thisiterstarttime)
			elseif onearcatatime_flag == 1
				thisiterstarttime = time()
				mostnegativearc = [a for a in 1:numarcs][argmin([arcredcosts[k,a] for a in 1:numarcs])]
				minreducedcost_k, shortestpatharcs, shortestpathnodes = arcredcosts[k, mostnegativearc], [mostnegativearc], []
				push!(dptimelist, time() - thisiterstarttime)
			end
			push!(min_rc_list, minreducedcost_k)

			#Add new arcs to current arc sets
			if minreducedcost_k < -0.001
				for a in shortestpatharcs
					if !(a in commArcSet[k])
						push!(commArcSet[k], a)
						i, j = arcs[a]
						push!(A_plus_k[k, i], a)
						push!(A_minus_k[k, j], a)
					end
				end
			end
		end

	 	min_rc = minimum(min_rc_list)
		println("Min reduced cost = ", min_rc) 

		#"Parallelize"
		shuffleddptimes = shuffle_partition(length(commodities))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		mvgsp_time_par += maximum(dptimelistsums)
		parallel_time += sum(dptimelistsums) -  maximum(dptimelistsums)

		#Add new variables 
		for k in commodities
			for a in setdiff(commArcSet[k], commArcSet_snapshot[k])

				#Create a new variable for the arc
				global x[k,a] = @variable(m, lower_bound = 0) #, upper_bound = 1)
				set_name(x[k,a], string("x[",k,",",a,"]")) 

				#Add to the objective
				set_objective_function(m, objective_function(m) + c_cg[a,k] * q[k] * x[k,a])

				#Add new variable to order constraints
				i = arcs[a][1]
				j = arcs[a][2]
				set_coefficient(smpconstraints.con_arccapacity[a], x[k,a], q[k])
				set_coefficient(smpconstraints.con_flowbalance[k,i], x[k,a], 1.0)
				set_coefficient(smpconstraints.con_flowbalance[k,j], x[k,a], -1.0)
				set_coefficient(smpconstraints.con_nodecapacity[i], x[k,a], q[k])
				set_coefficient(smpconstraints.con_nodecapacity[j], x[k,a], q[k])
			end
		end

		#Termination criterion
		if min_rc >= -0.001
			println("No negative reduced costs found!")
			break
		end

		mvg_iteration += 1

	end

	fullalgtime = time() - fullalgstarttime - parallel_time

	return getobjectivevalue(m), mvg_iteration, commArcSet, A_plus_k, A_minus_k, smp_time, mvgsp_time_par, fullalgtime, x

end

#----------------------------------------------------------------------------------------------------------#

function orderpathinitialization(c_cg)

	delta = Dict()
	pathcost = Dict()
	pathSet = Dict()
	for k in commodities
		pathSet[k] = []
		pathcost[k] = []
	end

	#Generate initial dummy paths
	arcindex = numarcs + 1
	for k in commodities
		p = 1
		push!(pathSet[k], p)
		for a in 1:numarcs_dummy
			if a == arcindex
				delta[k,a,p] = 1
			else
				delta[k,a,p] = 0
			end
		end
		push!(pathcost[k], c_cg[arcindex, k])
	    arcindex += 1
	end

	return delta, pathcost, pathSet

end

#----------------------------------------------------------------------------------------------------------#

function restrictedmasterproblem(pathSet, pathcost, delta)

	m = Model(Gurobi.Optimizer)
	set_optimizer_attribute(m, "OutputFlag", 0)

	@variable(m, 0 <= y[k in commodities, p in pathSet[k]]) #Upper bound redundant

	@objective(m, Min, sum(sum(pathcost[k][p] * q[k] * y[k,p] for p in pathSet[k]) for k in commodities) )

	@constraint(m, assignment[k in commodities], sum(y[k,p] for p in pathSet[k]) == 1)
	@constraint(m, arccapacity[a = 1:numarcs], sum(sum(q[k] * delta[k,a,p] * y[k,p] for p in pathSet[k]) for k in commodities) <= d[a])
	@constraint(m, nodecapacity[n in nodes], sum(sum(sum(q[k] * delta[k,a,p] * y[k,p] for a in union(A_plus[n], A_minus[n])) for p in pathSet[k]) for k in commodities) <= qnode[n])

	rmpconstraints = (con_arccapacity = arccapacity,
		con_assignment = assignment,
		con_nodecapacity = nodecapacity
		)

	return m, y, rmpconstraints

end

#----------------------------------------------------------------------------------------------------------#

function pathbasedcolumngeneration!(pathSet, pathcost, delta)

	#Initialize
	pbcg_iteration = 1
	rmp_time, pp_time, pp_time_par = 0, 0, 0
	parallel_time = 0
	listlength = convert(Int64, ceil(length(commodities)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)

	#Create the sparse master problem
	m, y, rmpconstraints = restrictedmasterproblem(pathSet, pathcost, delta)

	fullalgstarttime = time()

	#Main loop
	while 1==1

		#Solve RMP
		optimize!(m)
		rmp_time += solve_time(m)
		rmpobj = getobjectivevalue(m)
		println("PBCG iteration ", pbcg_iteration, " RMP objective = ", rmpobj)

		#Snapshot of current arcs
		pathSet_snapshot = deepcopy(pathSet)

		#Get dual values
		alpha = getdual.(rmpconstraints.con_arccapacity)
		beta = getdual.(rmpconstraints.con_assignment)
		iota = getdual.(rmpconstraints.con_nodecapacity)

		#Find arc reduced costs
		arcredcosts = Dict()
		for k in commodities, a in 1:numarcs
			i, j = arcs[a]
			arcredcosts[k,a] = c[a, k] * q[k] - q[k] * alpha[a] - q[k] * iota[i] - q[k] * iota[j]
		end

		#Solve pricing subproblem
		min_rc_list = []
		addarcs = []
		dptimelist = []
		for k in commodities
			thisiterstarttime = time()
			minreducedcost_k, shortestpathnodes, shortestpatharcs = findshortestpath_mcf(k, arcredcosts, numnodes, numarcs, arcs, arcLookup, Origin, Destination)
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
				newpathcost = sum(c[a,k] for a in shortestpatharcs)
				push!(pathcost[k], newpathcost)

				#Add new path variable
				global y[k,p] = @variable(m, lower_bound = 0) #, upper_bound = 1)
				set_name(y[k,p], string("y[",k,",",p,"]")) 
				set_objective_function(m, objective_function(m) + pathcost[k][p] * q[k] * y[k,p])
				for a in shortestpatharcs
					set_coefficient(rmpconstraints.con_arccapacity[a], y[k,p], q[k])
				end
				for n in nodes
					totalcoeff = 0
					for a in intersect(shortestpatharcs, union(A_plus[n], A_minus[n]))
						totalcoeff += q[k]
					end
					set_coefficient(rmpconstraints.con_nodecapacity[n], y[k,p], totalcoeff)
				end
				set_coefficient(rmpconstraints.con_assignment[k], y[k,p], 1.0)
				
			end
		end

	 	min_rc = minimum(min_rc_list)
		println("Min reduced cost = ", min_rc) 

		#"Parallelize"
		shuffleddptimes = shuffle_partition(length(commodities))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		pp_time_par += maximum(dptimelistsums)
		parallel_time += sum(dptimelistsums) -  maximum(dptimelistsums)

		#Termination criterion
		if min_rc >= -0.001
			println("No negative reduced costs found!")
			break
		end

		pbcg_iteration += 1

	end

	fullalgtime = time() - fullalgstarttime - parallel_time

	return getobjectivevalue(m), pbcg_iteration, pathSet, pathcost, delta, rmp_time, pp_time_par, fullalgtime, y

end

#----------------------------------------------------------------------------------------------------------#

function solvepathbasedmcfmodel(lprelaxation_flag, pathSet, pathcost, delta)

	m = Model(Gurobi.Optimizer)
	set_optimizer_attribute(m, "TimeLimit", 60*60*95)
	set_optimizer_attribute(m, "OutputFlag", 0)
	set_optimizer_attribute(m, "MIPGap", 0.01)

	if lprelaxation_flag == 1
		@variable(m, 0 <= y[k in commodities, p in pathSet[k]] <= 1)
	elseif lprelaxation_flag == 0
		@variable(m, y[k in commodities, p in pathSet[k]], Bin)
	end

	@objective(m, Min, sum(sum(pathcost[k][p] * q[k] * y[k,p] for p in pathSet[k]) for k in commodities) )

	@constraint(m, assignment[k in commodities], sum(y[k,p] for p in pathSet[k]) == 1)
	@constraint(m, arccapacity[a = 1:numarcs], sum(sum(q[k] * delta[k,a,p] * y[k,p] for p in pathSet[k]) for k in commodities) <= d[a])
	@constraint(m, nodecapacity[n in nodes], sum(sum(sum(q[k] * delta[k,a,p] * y[k,p] for a in union(A_plus[n], A_minus[n])) for p in pathSet[k]) for k in commodities) <= qnode[n])

	status_ip = optimize!(m)

	println("Objective = ", getobjectivevalue(m), " in ", solve_time(m), " seconds")

	return getobjectivevalue(m), solve_time(m)

end

#----------------------------------------------------------------------------------------------------------#

function generatefullarcssets()

	commArcSet_full = [union(1:numarcs, numarcs+k) for k in commodities]
	A_minus_k_full, A_plus_k_full = Dict(), Dict()
	for i in nodes, k in commodities
		A_minus_k_full[k, i] = copy(A_minus[i])
		A_plus_k_full[k, i] = copy(A_plus[i])
	end

	return commArcSet_full, A_plus_k_full, A_minus_k_full 

end

#----------------------------------------------------------------------------------------------------------#

function writeresults(solutionmethod, lp_obj, ip_obj, mp_time, sp_time, ip_time, full_time, numiter, arccount, pathcount, first_flag)

	#solutionmethod, lp_obj, ip_obj, mp_time, sp_time, ip_time, numiter, arccount, pathcount, first_flag = "mvg", mvg_lp, mvg_obj, smp_time, mvgsp_time_par, mvgip_time, mvg_iterations, sum(length(commArcSet_converged[k]) for k in commodities), 0, 1

	df = DataFrame(ID = [runid],
				numcommodities = [numcom],
				numnodes = [numnodes],
				radius = [radius],
				numarcs = [numarcs],
				method = [solutionmethod],
				lp_obj = [lp_obj],
				ip_obj = [ip_obj],
				mp_time = [mp_time],
				sp_time = [sp_time],
				ip_time = [ip_time],
				full_time = [full_time],
				iterations = [numiter],
				finalarccount = [arccount],
				finalpathcount = [pathcount]		
	           )

	if first_flag == 1
		CSV.write(string("outputs/algorithmcomparison/mcf_exp", runid, ".csv"), df) 
	else 
		CSV.write(string("outputs/algorithmcomparison/mcf_exp", runid, ".csv"), df, append=true)
	end

end