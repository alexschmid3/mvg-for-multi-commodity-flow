
function checkoptimizationoutput(terminationstatus, solvetime, hasvalues, currgamma, bestinfeasiblegamma, bestfeasiblegamma, timegoal)

    if (terminationstatus == MOI.OPTIMAL) & (solvetime < timegoal)
        println("Too short = ", solvetime)
        return 0, solvetime, round((currgamma + bestinfeasiblegamma) / 2, digits=0), bestinfeasiblegamma, currgamma

    elseif (terminationstatus == MOI.OPTIMAL) & (solvetime >= timegoal)
        println("Good instance = ", solvetime)
        return 1, solvetime, currgamma, bestinfeasiblegamma, currgamma

    elseif (terminationstatus == MOI.TIME_LIMIT) && (hasvalues)        
        println("Good instance = ", solvetime)
        return 1, solvetime, currgamma, bestinfeasiblegamma, currgamma

    elseif (terminationstatus == MOI.INFEASIBLE_OR_UNBOUNDED) || (terminationstatus == MOI.INFEASIBLE)
        println("Infeasible")
        return 0, "Inf", round((currgamma + bestfeasiblegamma) / 2, digits=0), currgamma, bestfeasiblegamma

    elseif (terminationstatus == MOI.TIME_LIMIT) && !(hasvalues)
        println("Infeasible by time limit")
        return 0, "Inf", round((currgamma + bestfeasiblegamma) / 2, digits=0), currgamma, bestfeasiblegamma

    elseif (terminationstatus == MOI.SOLUTION_LIMIT) && (hasvalues) & (solvetime >= timegoal)
        println("Good instance = ", solvetime)
        return 1, solvetime, currgamma, bestinfeasiblegamma, currgamma

    elseif (terminationstatus == MOI.SOLUTION_LIMIT) && (hasvalues) & (solvetime < timegoal)
        println("Too short = ", solvetime)
        return 0, solvetime, round((currgamma + bestinfeasiblegamma) / 2, digits=0), bestinfeasiblegamma, currgamma

    end

end

#--------------------------------------------------------------------------------------------#

function writeinstancetofile(iter, gamma_arc, gamma_node, numarcs, total_time, goodinstance_flag)

    df = DataFrame(experiment_id = [runid],
			randomseed = [randomseedval],
			numcom = [numcom], 
			numnodes = [numnodes], 
			radius = [radius], 
			gamma_arc = [gamma_arc], 
			gamma_node = [gamma_node],
			opt_gap = [opt_gap], 
			numarcs = [numarcs], 
			total_time = [total_time],
			goodinstance_flag = [goodinstance_flag]
           )

    if iter == 1
        CSV.write(instancefilename, df)
    else
        CSV.write(instancefilename, df, append=true)
    end

end

#--------------------------------------------------------------------------------------------#

function arccapacitytuning(gamma_arc_init, gamma_node, timegoal, maxtuningiterations)

    coordinates, commodities, nodes, arcs, arcLookup, numarcs, A_minus, A_plus, c, b, q, Origin, Destination = randomizeinstance(randomseedval)
    arcperturbation, nodeperturbation = randomizeperturbations(numarcs, nodes, randomseedval)

    #Initialize capacity tuning
    currgamma = gamma_arc_init
    goodinstance_flag = 0
    bestfeasiblegamma = 0
    bestinfeasiblegamma = 50*gamma_arc_init
    iter = 1

    #Terminate if no tuning iterations
    if maxtuningiterations == 0
        d, p = setcapacities(currgamma, gamma_node, q, numarcs, nodes, arcperturbation, nodeperturbation)
        mcfinstance = (commodities=commodities, numarcs=numarcs, nodes=nodes, c=c, q=q, p=p, d=d, b=b, A_plus=A_plus, A_minus=A_minus, coordinates=coordinates, arcs=arcs, arcLookup=arcLookup, Origin=Origin, Destination=Destination)
        return currgamma, 1, mcfinstance, iter
    end

    #Capacity tuning loop
    while 1==1
        
        println("---------- Iteration $iter arc gamma = ", currgamma, " ----------")

        #Set the capacities based on the gammas
        d, p = setcapacities(currgamma, gamma_node, q, numarcs, nodes, arcperturbation, nodeperturbation)
        mcfinstance = (commodities=commodities, numarcs=numarcs, nodes=nodes, c=c, q=q, p=p, d=d, b=b, A_plus=A_plus, A_minus=A_minus, coordinates=coordinates, arcs=arcs, arcLookup=arcLookup, Origin=Origin, Destination=Destination)

        #Solve the instance 
        model, terminationstatus, solvetime, hasvalues = solvemcfinstance(mcfinstance, 0, timegoal, 0, "tuning", [])
        
        #Check Gurobi output and update gammas accordingly
        println("$terminationstatus, $solvetime, $hasvalues, $currgamma, $bestinfeasiblegamma, $bestfeasiblegamma, $timegoal")
        goodinstance_flag, total_time, currgamma, bestinfeasiblegamma, bestfeasiblegamma = checkoptimizationoutput(terminationstatus, solvetime, hasvalues, currgamma, bestinfeasiblegamma, bestfeasiblegamma, timegoal)

        #Write instance to file
        writeinstancetofile(iter, currgamma, gamma_node, numarcs, total_time, 0.5*goodinstance_flag)

        #Iterate
        iter += 1

        #Check termination criterion
        if (goodinstance_flag == 1) || (iter > maxtuningiterations)
            return bestfeasiblegamma, goodinstance_flag, mcfinstance, iter
        end

    end

    return bestfeasiblegamma, goodinstance_flag, (), iter

end

#--------------------------------------------------------------------------------------------#

function nodecapacitytuning(gamma_arc, gamma_node_init, timegoal, maxtuningiterations, startiter)

    coordinates, commodities, nodes, arcs, arcLookup, numarcs, A_minus, A_plus, c, b, q, Origin, Destination = randomizeinstance(randomseedval)
    arcperturbation, nodeperturbation = randomizeperturbations(numarcs, nodes, randomseedval)

    #Initialize capacity tuning
    currgamma = gamma_node_init
    goodinstance_flag = 0
    bestfeasiblegamma = 0
    bestinfeasiblegamma = 50*gamma_node_init
    iter = startiter + 1

    #Terminate if no tuning iterations
    if maxtuningiterations == 0
        d, p = setcapacities(gamma_arc, currgamma, q, numarcs, nodes, arcperturbation, nodeperturbation)
        mcfinstance = (commodities=commodities, numarcs=numarcs, nodes=nodes, c=c, q=q, p=p, d=d, b=b, A_plus=A_plus, A_minus=A_minus, coordinates=coordinates, arcs=arcs, arcLookup=arcLookup, Origin=Origin, Destination=Destination)
        return currgamma, 1, mcfinstance
    end

    #Capacity tuning loop
    while 1==1
        
        println("---------- Iteration $iter node gamma = ", currgamma, " ---------")

        #Set the capacities based on the gammas
        d, p = setcapacities(gamma_arc, currgamma, q, numarcs, nodes, arcperturbation, nodeperturbation)
        mcfinstance = (commodities=commodities, numarcs=numarcs, nodes=nodes, c=c, q=q, p=p, d=d, b=b, A_plus=A_plus, A_minus=A_minus, coordinates=coordinates, arcs=arcs, arcLookup=arcLookup, Origin=Origin, Destination=Destination)

        #Solve the instance 
        model, terminationstatus, solvetime, hasvalues = solvemcfinstance(mcfinstance, 0, timegoal, 1, "tuning", [])
        
        #Check Gurobi output and update gammas accordingly
        goodinstance_flag, total_time, currgamma, bestinfeasiblegamma, bestfeasiblegamma = checkoptimizationoutput(terminationstatus, solvetime, hasvalues, currgamma, bestinfeasiblegamma, bestfeasiblegamma, timegoal)

        #Write instance to file
        writeinstancetofile(iter, gamma_arc, currgamma, numarcs, total_time, goodinstance_flag)

        iter += 1

        #Check termination criterion
        if (goodinstance_flag == 1) || (iter > maxtuningiterations)
            return bestfeasiblegamma, goodinstance_flag, mcfinstance
        end

    end

    return bestfeasiblegamma, goodinstance_flag, ()

end

#--------------------------------------------------------------------------------------------#

function capacitytuning(gamma_arc_init, gamma_node_init, timegoal_arc, timegoal_node, maxtuningiterations)

    gamma_arc, goodinstance_flag_arc, arcmcfinstance, startiter = arccapacitytuning(gamma_arc_init, gamma_node_init, timegoal_arc, maxtuningiterations)
    if goodinstance_flag_arc == 1
        gamma_node, goodinstance_flag, finalmcfinstance = nodecapacitytuning(gamma_arc, gamma_node_init, timegoal_node, maxtuningiterations, startiter)
    elseif gamma_arc > 1e-4
        println("No good instance found during arc tuning, but proceeding to node tuning")
        gamma_node, goodinstance_flag, finalmcfinstance = nodecapacitytuning(gamma_arc, gamma_node_init, timegoal_node, maxtuningiterations, startiter) 
    else
        println("Sorry, no feasible instance found during arc tuning")
        return gamma_arc, gamma_node_init, arcmcfinstance, 0
    end

    if goodinstance_flag == 1
        return gamma_arc, gamma_node, finalmcfinstance, goodinstance_flag
    else
        println("Sorry, no good instance found during node tuning")
        return gamma_arc, gamma_node, finalmcfinstance, goodinstance_flag
    end

end
