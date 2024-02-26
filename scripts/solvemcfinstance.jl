
function solvemcfinstance(mcfinstance, lprelaxflag, timegoal, nodecapacity_flag, usage, reducedarcs)

    #Write MCF IP
    m = Model(Gurobi.Optimizer)

    set_optimizer_attribute(m, "TimeLimit", timegoal+1)
    set_optimizer_attribute(m, "OutputFlag", 0)
    set_optimizer_attribute(m, "LogToConsole", 0)
    set_optimizer_attribute(m, "MIPGap", opt_gap)
    
    if usage == "reducedsolve"
        @variable(m, x[k in mcfinstance.commodities, a in reducedarcs.A[k]], Bin)
    elseif lprelaxflag == 1
        @variable(m, x[k in mcfinstance.commodities, a in 1:mcfinstance.numarcs] >= 0)
    elseif lprelaxflag == 0
        @variable(m, x[k in mcfinstance.commodities, a in 1:mcfinstance.numarcs], Bin)
    end
    
    if usage == "reducedsolve"
        @objective(m, Min, sum(sum(mcfinstance.c[a,k] * mcfinstance.q[k] * x[k,a] for a in reducedarcs.A[k]) for k in mcfinstance.commodities) )
        @constraint(m, flowbalance[k in mcfinstance.commodities, n in mcfinstance.nodes], sum(x[k, a] for a in reducedarcs.A_plus[k,n]) - sum(x[k, a] for a in reducedarcs.A_minus[k,n]) == mcfinstance.b[n,k])
        @constraint(m, arccapacitycon[a = 1:mcfinstance.numarcs], 0 <= mcfinstance.d[a])
        for a in 1:mcfinstance.numarcs, k in mcfinstance.commodities
            if a in reducedarcs.A[k]
                set_normalized_coefficient(arccapacitycon[a], x[k,a], mcfinstance.q[k])
            end
        end   
        if nodecapacity_flag == 1
            @constraint(m, nodecapacitycon[n in mcfinstance.nodes], sum(sum(mcfinstance.q[k] * x[k, a] for a in union(reducedarcs.A_plus[k,n], reducedarcs.A_minus[k,n])) for k in mcfinstance.commodities) <= mcfinstance.p[n])
        end
    else
        @objective(m, Min, sum(sum(mcfinstance.c[a,k] * mcfinstance.q[k] * x[k,a] for a in 1:mcfinstance.numarcs) for k in mcfinstance.commodities) )
        @constraint(m, flowbalance[k in mcfinstance.commodities, i in mcfinstance.nodes], sum(x[k, a] for a in mcfinstance.A_plus[i]) - sum(x[k, a] for a in mcfinstance.A_minus[i]) == mcfinstance.b[i,k])
        @constraint(m, arccapacitycon[a = 1:mcfinstance.numarcs], sum(mcfinstance.q[k] * x[k, a] for k in mcfinstance.commodities) <= mcfinstance.d[a])
        if nodecapacity_flag == 1
            @constraint(m, nodecapacitycon[n in mcfinstance.nodes], sum(sum(mcfinstance.q[k] * x[k, a] for a in union(mcfinstance.A_plus[n], mcfinstance.A_minus[n])) for k in mcfinstance.commodities) <= mcfinstance.p[n])
        end
    end
    
    #Solve IP
    status_ip = optimize!(m)

    #Check if needs more time
    if (usage == "tuning") & (termination_status(m) == MOI.TIME_LIMIT) && !(has_values(m))
        println("Needs more time...")
        set_optimizer_attribute(m, "TimeLimit", timegoal*3)
        set_optimizer_attribute(m, "SolutionLimit", 1)
        optimize!(m)
    end

    if usage == "tuning"

        return m, termination_status(m), solve_time(m), has_values(m)
   
    elseif usage in ["fullsolve", "reducedsolve"]
        
        if (termination_status(m) == MOI.OPTIMAL) || ((termination_status(m) == MOI.TIME_LIMIT) && (hasvalues(m))  )
            if lprelaxflag == 1
                println("LP objective = ", objective_value(m))
            elseif lprelaxflag == 0
                println("IP objective = ", objective_value(m))
            end
            return objective_value(m), value.(x), termination_status(m), solve_time(m), has_values(m)
        elseif (termination_status(m) == MOI.INFEASIBLE_OR_UNBOUNDED) || (termination_status(m) == MOI.INFEASIBLE)
            if lprelaxflag == 1
                println("LP infeasible")
            elseif lprelaxflag == 0
                println("IP infeasible")
            end
            return 1e10, [], termination_status(m), solve_time(m), has_values(m)
        elseif (termination_status(m) == MOI.TIME_LIMIT) && !(hasvalues(m))
            if lprelaxflag == 1
                println("LP infeasible")
            elseif lprelaxflag == 0
                println("IP infeasible")
            end
            return 1e10, [], termination_status(m), solve_time(m), has_values(m)
        end

    end

end

