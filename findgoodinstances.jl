
using JuMP, Gurobi, Plots, Random, CSV, DataFrames, Statistics, LinearAlgebra

include("scripts/randomizeinstance.jl")
include("scripts/randomizecapacities.jl")
include("scripts/capacitytuning.jl")
include("scripts/solvemcfinstance.jl")
include("scripts/writerunresults.jl")
include("scripts/multiarcgeneration.jl")
include("scripts/columngeneration.jl")
include("scripts/solvepathmcfinstance.jl")

#-------------------------------PARAMETERS--------------------------------#

runid = 90 #ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
println("runid = $runid")
paramsfilename = "data/findgoodinstances_2024.csv"
expparms = CSV.read(paramsfilename, DataFrame)
randomseedval = expparms[runid, 2]		
numcom = expparms[runid, 3]
numnodes = expparms[runid, 4]		
radius = expparms[runid, 5]	
gamma_arc_init = 150000	#expparms[runid, 6]	
gamma_node_init = 938 #expparms[runid, 7]
opt_gap = expparms[runid, 8]

println("Seed = ", randomseedval)	

maxorder = 10
minorder = 1
destdistpercentile = 0.2
mindistanceperturb = 0.9
maxdistanceperturb = 1.1
mincapacityperturb = 0.5
maxcapacityperturb = 1.5
commodities = 1:numcom

#timegoal_arcdict = Dict(100=>15, 200=>45, 400=>60)
#timegoal_nodedict = Dict(100=>40, 200=>90, 400=>120)
timegoal_arcdict = Dict(100=>30, 200=>60, 400=>75)
timegoal_nodedict = Dict(100=>120, 200=>120, 400=>180)

#Time goals for solving MCF instances
timegoal_arc = timegoal_arcdict[numcom] #numnodes/3 #0*60
timegoal_node = timegoal_nodedict[numcom] #numnodes/2 #30*60
maxtuningiterations = 0 #20

#Algorithm parameters
iptimelimit = 60*2

#Output file name
instancefilename = string("outputs/instancetuning_exp", runid, ".csv")
outputfilename = string("outputs/mcfalgorithms_exp", runid, ".csv")

#---------------------------GENERATE INSTANCE-----------------------------#

gamma_arc, gamma_node, mcfinstance, foundgoodinstance_flag = capacitytuning(gamma_arc_init, gamma_node_init, timegoal_arc, timegoal_node, maxtuningiterations)
# gamma_arc, gamma_node, mcfinstance, foundgoodinstance_flag = capacitytuning(30000, 150, 1, 1, 0)

#----------------------------SOLVE INSTANCE-------------------------------#

#include("scripts/drawmap.jl")
#drawmap("networkmap.png", mcfinstance, 2000, 2000)

fullarccount = sum(mcfinstance.numarcs for k in mcfinstance.commodities)
fullpathcount = 0

#LP
println("---------- LP ----------")
obj_lp, x_lp, termination_lp, solvetime_lp, hasvalues_lp = solvemcfinstance(mcfinstance, 1, iptimelimit, 1, "fullsolve", [])
writerunresults(outputfilename, "LP", mcfinstance, obj_lp, solvetime_lp, 0, 0, solvetime_lp, 0, fullarccount, fullpathcount, 1)

#IP
println("---------- IP ----------")
obj_ip, x_ip, termination_ip, solvetime_ip, hasvalues_ip = solvemcfinstance(mcfinstance, 0, iptimelimit, 1, "fullsolve", [])
writerunresults(outputfilename, "IP", mcfinstance, obj_ip, 0, 0, solvetime_ip, solvetime_ip, 0, fullarccount, fullpathcount, 0)

#MAG
println("---------- MAG ---------")
c_mag, d_mag, magarcs, numarcs_dummy, dummydeletions = colgeninitialize(mcfinstance)
obj_mag, mag_iterations, magarcs, smp_time, magsp_time_par, fullalgtime = multiarcgeneration!(mcfinstance, magarcs, c_mag, d_mag, numarcs_dummy, dummydeletions)
writerunresults(outputfilename, "MAG", mcfinstance, obj_mag, smp_time, magsp_time_par, 0, fullalgtime, mag_iterations, sum(length(magarcs.A[k]) for k in mcfinstance.commodities), 0, 0)

println("-------- MAG IP --------")
obj_magip, x_ip, termination_ip, solvetime_magip, hasvalues_ip = solvemcfinstance(mcfinstance, 0, iptimelimit, 1, "reducedsolve", magarcs)
writerunresults(outputfilename, "MAGIP", mcfinstance, obj_magip, 0, 0, solvetime_magip, solvetime_magip, 0, sum(length(magarcs.A[k]) for k in mcfinstance.commodities), 0, 0)

#CG
println("---------- CG ----------")
delta, pathcost, pathSet = orderpathinitialization(c_mag, mcfinstance, numarcs_dummy)
obj_cg, cg_iterations, pathSet_converged, pathcost_converged, delta_converged, rmp_time, cgsp_time_par, cg_fulltime, y_pb = columngeneration!(pathSet, pathcost, delta, mcfinstance, numarcs_dummy)
writerunresults(outputfilename, "CG", mcfinstance, obj_cg, rmp_time, cgsp_time_par, 0, cg_fulltime, cg_iterations, 0, sum(length(pathSet_converged[k]) for k in mcfinstance.commodities), 0)

println("-------- CG IP ---------")
obj_cgip, solvetime_cgip = solvepathmcfinstance(0, pathSet_converged, pathcost_converged, delta_converged, mcfinstance)
writerunresults(outputfilename, "CGIP", mcfinstance, obj_cgip, 0, 0, solvetime_cgip, solvetime_cgip, 0, 0, sum(length(pathSet_converged[k]) for k in mcfinstance.commodities), 0)
