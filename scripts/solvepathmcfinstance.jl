
function solvepathmcfinstance(lprelaxation_flag, pathSet, pathcost, delta, mcfinstance)

	m = Model(Gurobi.Optimizer)
	set_optimizer_attribute(m, "TimeLimit", 60*60)
	set_optimizer_attribute(m, "OutputFlag", 0)
	set_optimizer_attribute(m, "MIPGap", 0.01)

	if lprelaxation_flag == 1
		@variable(m, 0 <= y[k in mcfinstance.commodities, p in pathSet[k]] <= 1)
	elseif lprelaxation_flag == 0
		@variable(m, y[k in mcfinstance.commodities, p in pathSet[k]], Bin)
	end

	@objective(m, Min, sum(sum(pathcost[k][p] * mcfinstance.q[k] * y[k,p] for p in pathSet[k]) for k in mcfinstance.commodities) )

	@constraint(m, assignment[k in mcfinstance.commodities], sum(y[k,p] for p in pathSet[k]) == 1)
	@constraint(m, arccapacity[a = 1:mcfinstance.numarcs], sum(sum(mcfinstance.q[k] * delta[k,a,p] * y[k,p] for p in pathSet[k]) for k in mcfinstance.commodities) <= mcfinstance.d[a])
	@constraint(m, nodecapacity[n in mcfinstance.nodes], sum(sum(sum(mcfinstance.q[k] * delta[k,a,p] * y[k,p] for a in union(mcfinstance.A_plus[n], mcfinstance.A_minus[n])) for p in pathSet[k]) for k in mcfinstance.commodities) <= mcfinstance.p[n])

	status_ip = optimize!(m)

	println("Objective = ", objective_value(m), " in ", solve_time(m), " seconds")

	return objective_value(m), solve_time(m)

end
