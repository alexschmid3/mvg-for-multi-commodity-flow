
function writerunresults(filename, solutionmethod, mcfinstance, obj, mp_time, sp_time, ip_time, full_time, numiter, arccount, pathcount, first_flag)

    df = DataFrame(ID = [runid],
				numcommodities = [numcom],
				numnodes = [numnodes],
				radius = [radius],
				numarcs = [mcfinstance.numarcs],
				method = [solutionmethod],
				obj = [obj],
				mp_time = [mp_time],
				sp_time = [sp_time],
				ip_time = [ip_time],
				full_time = [full_time], 
				iterations = [numiter],
				finalarccount = [arccount],
				finalpathcount = [pathcount]		
	           )

	if first_flag == 1
		CSV.write(filename, df) 
	else 
		CSV.write(filename, df, append=true)
	end

end