
using JuMP, Gurobi, Plots, Random, CSV, DataFrames, Statistics, LinearAlgebra

include("scripts/mcfinstancegeneration.jl")
include("scripts/algorithms.jl")

#-------------------------------------- PARAMETERS -------------------------------------#

rowid = 1 #ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
paramsfilename = "data/goodinstances.csv"
expparms = CSV.read(paramsfilename, DataFrame)
runid = expparms[rowid, 1] 	
randomseedval = expparms[rowid, 2] 	
numcom = expparms[rowid, 3]
numnodes = expparms[rowid, 4]		
radius = expparms[rowid, 5]	
gamma_arc = expparms[rowid, 6]	
gamma_node = expparms[rowid, 7]	
opt_gap = expparms[rowid, 8] 

maxorder = 10
minorder = 1
destdistpercentile = 0.2
mindistanceperturb = 0.9
maxdistanceperturb = 1.1
mincapacityperturb = 0.5
maxcapacityperturb = 1.5

outputfilename = string("outputs/allruns_exp", runid, ".csv")

time()

#-------------------------------GENERATE RANDOM INSTANCE--------------------------------#

Random.seed!(randomseedval)
coordinates, commodities, nodes, arcs, arcLookup, numarcs, A_minus, A_plus, c, b, q, Origin, Destination = createinstance_mcf(radius, destdistpercentile, maxdistanceperturb, mindistanceperturb)
arcperturbation, nodeperturbation = randomizecapacities_mcf(numarcs, nodes, maxcapacityperturb, mincapacityperturb)
d, qnode = setcapacities_mcf(gamma_arc, gamma_node, q, numarcs, nodes, arcperturbation, nodeperturbation)

#----------------------------------------- MVG -----------------------------------------#

onearcatatime_flag = 0
c_cg, d_cg, commArcSet, A_minus_k, A_plus_k, numarcs_dummy = colgeninitialize(c, d, numarcs)
mvg_lp, mvg_iterations, commArcSet_mvg, A_plus_k_mvg, A_minus_k_mvg, smp_time, mvgsp_time_par, mvg_fulltime = multivariablegeneration!(commArcSet, A_plus_k, A_minus_k, onearcatatime_flag, c_cg)
mvg_obj, mvgip_time = solvemcfmodel(0, commArcSet_mvg, A_plus_k_mvg, A_minus_k_mvg, c_cg)
writeresults("mvg", mvg_lp, mvg_obj, smp_time, mvgsp_time_par, mvgip_time, mvg_fulltime, mvg_iterations, sum(length(commArcSet_mvg[k]) for k in commodities), 0, 1)

#----------------------------------------- SVG -----------------------------------------#

onearcatatime_flag = 1
c_cg, d_cg, commArcSet, A_minus_k, A_plus_k, numarcs_dummy = colgeninitialize(c, d, numarcs)
svg_lp, svg_iterations, commArcSet_svg, A_plus_k_svg, A_minus_k_svg, smp_time, svgsp_time_par, svg_fulltime = multivariablegeneration!(commArcSet, A_plus_k, A_minus_k, onearcatatime_flag, c_cg)
svg_obj, svgip_time = solvemcfmodel(0, commArcSet_svg, A_plus_k_svg, A_minus_k_svg, c_cg)
writeresults("svg", svg_lp, svg_obj, smp_time, svgsp_time_par, svgip_time, svg_fulltime, svg_iterations, sum(length(commArcSet_svg[k]) for k in commodities), 0, 0)

#----------------------------------------- PBCG ----------------------------------------#

delta, pathcost, pathSet = orderpathinitialization(c_cg)
pbcg_lp, pbcg_iterations, pathSet_converged, pathcost_converged, delta_converged, rmp_time, pp_time_par, pbcg_fulltime = pathbasedcolumngeneration!(pathSet, pathcost, delta)
pbcg_obj, pbcgip_time = solvepathbasedmcfmodel(0, pathSet_converged, pathcost_converged, delta_converged)
writeresults("pbcg", pbcg_lp, pbcg_obj, rmp_time, pp_time_par, pbcgip_time, pbcg_fulltime, pbcg_iterations, 0, sum(length(pathSet[k]) for k in commodities), 0)

#----------------------------------------- IP ------------------------------------------#

commArcSet_full, A_plus_k_full, A_minus_k_full = generatefullarcssets()
ip_obj, ip_time = solvemcfmodel(0, commArcSet_full, A_plus_k_full, A_minus_k_full, c)
writeresults("ip", 0, ip_obj, 0, 0, ip_time, ip_time, 0, sum(length(commArcSet_full[k]) for k in commodities), 0, 0)

#---------------------------------------------------------------------------------------#

println("Done!")