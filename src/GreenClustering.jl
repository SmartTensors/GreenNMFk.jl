
"""
Purpose of function

$(DocumentFunction.documentfunction(clustering_the_solutions;
argtext=Dict("number_of_sources"=>"",
			""=>"",
			""=>"",
            ""=>"")))
Returns:
- 
"""

function clustering_the_solutions(number_of_sources::Number, nd::Number, sol::Matrix, normF::Vector, Qyes::Number)

	# Define initial variables
	local number_of_clust_sim, vect_index, sol1, cent, ind
	local mean_savg = 0
	minsil_old = -2

	# Determine the number of quantiles
	if (Qyes == 1)
		steps = 10
		quants = Base.quantile(normF, linspace(0.2, 0.01, steps))
	else
		steps = 1
		quants = maximum(normF)
	end
	
	# Calculation of silhouettes in different quantiles
	for p = 1:steps	
		if GreenNMFk.io_level > 1 println("  Quantile $(p)/$(steps)...") end
		
		number_of_clust_sim = p
		ind = find(normF[normF .<= quants[p]])
		sol1 = sol[ind,:]
		
		# Solution space for:
	    #                   1  2  3  4  5  6  7  8  9 10 11 12
	    # three sources = [Ux Dx Dy A1 X1 Y1 A2 X2 Y2 A3 X3 Y3]
	    # two sources   = [Ux Dx Dy A1 X1 Y1 A2 X2 Y2]
	    # one source    = [Ux Dx Dy A1 X1 Y1]
		
		just_source = sol1[:,4:end]
		sources_3D = Array{Float64}(3, Int(size(just_source, 2)/3), size(just_source, 1))
		
		local col, col_sources # If not defined, they go out of scope in for-loop
		
		for J = 1:size(sol1, 1)
			hold = just_source[J,:]
			
			for kJ = 1:Int(size(just_source, 2)/3)
				if kJ == 1
					col = hold[kJ*3-2:kJ*3]
				else
					col = [col; hold[kJ*3-2:kJ*3]]
				end
			end
			
			sources_3D[:,:,J] = col'
			
			if (J == 1)
				col_sources = col
			else
				col_sources = [col_sources; col]
			end
		end
		
		# Ludmil clustering
		if GreenNMFk.io_level > 1 println("    - Ludmil clustering") end
		_, vect_index, cent = ludmil_cluster(sources_3D, col_sources, 100)
		vect_index = Array{Int64,1}(vec(vect_index))
		
		# Silhouette clustering
		if GreenNMFk.io_level > 1 println("    - Silhouette clustering") end
		clusterassignments, centroids = NMFk.clustersolutions_old([vec(col_sources)'], true)
		W, H, clustersilhouettes = NMFk.finalize(col_sources, [vect_index'], clusterassignments, true)

		savg = grpstats(clustersilhouettes, vect_index)
		minsil = savg
		
		if ((mean(minsil) < (minsil_old / 2)) || (size(ind, 1) < 5) || (min(minsil) > 0.95))
			error("Broken at quantile $(p)")
			break
		end
		
		minsil_old = mean(minsil)
	end
	
	idx1 = (vect_index .== 1) # Boolean array, true where vect_index[i] == 1
	sol_matrix = sol1[vect_index[idx1],:]
	avg_sol = [mean(sol_matrix[:,i]) for i=1:size(sol_matrix,2)]
	
	solution = Array{Float64}(number_of_sources,6)

	# Condense into single line?
	for i = 1:number_of_sources
		solution[i,:] = [cent[i,:] avg_sol[1:3]]
	end
	
	reconstr = mean(normF[ind])
	mean_savg = mean(savg)
	
	# Save a JLD file with function results
	if save_output
		outfile = "Solution_$(nd)det_$(number_of_sources)sources.jld"
		
		info("Saving results to $(working_dir)/$(outfile)")
		JLD.save(joinpath(working_dir, outfile), "Solution", solution, "VectIndex", vect_index, "Cent", cent, "ss", ss, "savg", savg, "reconstr", reconstr, "mean_savg", mean_savg, "number_of_clust_sim", number_of_clust_sim)
	end
	
	return solution, vect_index, cent, reconstr, mean_savg, number_of_clust_sim
end
		
		