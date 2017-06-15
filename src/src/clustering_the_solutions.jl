function clustering_the_solutions(number_of_sources, nd, sol, normF, Qyes)
	minsil_old = -2
	
	# Determine the number of quantiles
	if (Qyes == 1)
		steps = 10
		quants = quantile(normF, linspace(0.2, 0.01, steps))
	else
		steps = 1
		quants = max(normF)
	end
	
	# Calculation of silhouettes in different quantiles
	for p = 1:steps
		println("Run $(p) of $(steps)")
		ind = find(normF <= quants(p))
		sol1 = sol(ind,:)
		
	    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	    #            1  2  3  4  5  6   7  8  9  10 11 12
	    #Each sol for three sources = [Ux Dx Dy A1 X1 Y1 A2 X2 Y2 A3 X3 Y3]
	    #Each sol for two sources   = [Ux Dx Dy A1 X1 Y1 A2 X2 Y2]
	    #Each sol for one source    = [Ux Dx Dy A1 X1 Y1]
		
		just_source = sol1(:,4:end)
		sources_3D = zeros(3, size(just_source, 2)/3, size(just_source, 1))
		
		for J = 1:size(sol1, 1)
			hold = just_source[J,:]
			
			for kJ = 1:size(just_source, 2) / 3
				if (kJ == 1)
					col = hold(kJ*3-2:kJ*3)
				else
					col = [col; hold(kJ*3-2:kJ*3)]
				end
			end
			
			sources_3D[:,:,J] = col'
			
			if (J == 1)
				col_sources = col
			else
				col_sources = [col_sources; col]
			end
		end
		
		[~, vect_index, cent] = ludmil_cluster(sources_3D, col_sources, 100)
		[ss, h] = silhouette[col_sources, vect_index]
		savg = grpstats[ss, vect_index]
		minsil = savg
		
		if ((mean(minsil) < (minsil_old / 2)) || (size(ind, 1) < 5) || (min(minsil) > 0.95))
			break
		end
		
		minsil_old = mean(minsil)
	end
	
	number_of_clust_sim = p
	
	idx1 = vect_index==1
	avg_sol = mean(sol1[vect_index[idx1],:])
	solution = zeros(number_of_sources,6)
	
	# Condense into single line?
	for i = 1:number_of_sources
		solution[i,:] = [cent[i,:] avg_sol[1:3]]
	end
	
	reconstr = mean(normF[ind])
	mean_savg = mean(savg)
	
	file_name1 = sprintf('./Results/Solution_%ddet_%dsources.mat',nd, number_of_sources);
	save(file_name1, 'Solution', 'VectIndex', 'Cent', 'ss', 'savg', 'reconstr', 'mean_savg','number_of_clust_sim');
end
		
		