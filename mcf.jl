
using JuMP, Gurobi, Plots, Random, CSV, DataFrames, Statistics, LinearAlgebra

include("scripts/mcfinstancegeneration.jl")

#-------------------------------GENERATE RANDOM INSTANCE--------------------------------#

runid = ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
println("runid = $runid")
paramsfilename = "data/findgoodinstances.csv"
expparms = CSV.read(paramsfilename, DataFrame)
randomseedval = expparms[runid, 2]		
numcom = expparms[runid, 3]
numnodes = expparms[runid, 4]		
radius = expparms[runid, 5]	
gamma_arc_init = expparms[runid, 6]	
gamma_node_init = expparms[runid, 7]
opt_gap = expparms[runid, 8]

println("Seed = ", randomseedval)	

Random.seed!(randomseedval)
maxorder = 10
minorder = 1
destdistpercentile = 0.2
mindistanceperturb = 0.9
maxdistanceperturb = 1.1
mincapacityperturb = 0.5
maxcapacityperturb = 1.5
commodities = 1:numcom

outputfilename = string("outputs/ipruntimes_exp", runid, ".csv")

#---------------------------TUNE GAMMA AND COMPLETE INSTANCE----------------------------#

maxiter = 20
timegoal1 = 5
 
Random.seed!(randomseedval)
gamma_arc, goodinstance_flag = findgoodinstance_arctuning(gamma_arc_init, gamma_node_init, opt_gap, maxiter, timegoal1, numcom, numnodes, maxorder, minorder, radius, destdistpercentile, maxdistanceperturb, mindistanceperturb, maxcapacityperturb, mincapacityperturb)

timegoal2 = 5

Random.seed!(randomseedval)
gamma_node, goodinstance_flag = findgoodinstance_nodetuning(gamma_arc, gamma_node_init, opt_gap, maxiter, timegoal2, numcom, numnodes, maxorder, minorder, radius, destdistpercentile, maxdistanceperturb, mindistanceperturb, maxcapacityperturb, mincapacityperturb)

#---------------------------------------------------------------------------------------#

println("Done!")