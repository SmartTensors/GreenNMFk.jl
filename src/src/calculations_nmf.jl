# Purpose: 	Calculation of the solution with j sources
#			Minimizes the function: ∑ᵢ(MixFn[i]- ∑ⱼ(Sources[i,j]))² 
# 
# Input:
#	number_of_sources - Number of sources
#	nd 	 - Number of detectors
#	Nsim - Number of simulations
#	aa	 - Multiple for initial random conditions
#	xD	 - Detector positions
#	t0	 - Initial time of sources
#	time - Time vector from t_initial to t_final
#	S	 - Observation matrix (returned from initial_conditions)
#	numT - Number of time points

# Returns:
#	
function calculations_nmf_v02(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT)
	
	GreenNMFk.log("\nRunning NMF calculation...")
	GreenNMFk.log("-----------------------------------------")
	
	sol = zeros(Nsim, 3 * number_of_sources + 3) # Solution matrix
	normF = zeros(Nsim, 1) # 
	normCut = 0
	Qyes = 0
	
	# Define the function that will be minimized
	# funF = ∑ᵢ(MixFn[i]- ∑ⱼ(Sources[i,j]))² 
	GreenNMFk.log("  Finding function to be minimized")
	
	# Create empty array for which each new function iteration will be appended
	# The idea is to avoid an infinite recursion, i.e.:
	# a = x->2x, b = x->x^2
	# a = x->a(x)+b(x) != x->2x+x^2, instead it tries to call itself
	fun_sum = Any[]
	
	for i=1:nd
		if number_of_sources == 1
			Mixfn = x -> source(time, x[4:6], xD[i,:], x[1:2], t0, x[3])
		else
			for d=1:number_of_sources
				if (d == 1)
					Mixfn   = x -> source(time, x[4:6], xD[i,:], x[1:2], t0, x[3])
				else
					mixfun2 = x -> source(time, x[d*3+1:d*3+3], xD[i,:], x[1:2], t0, x[3])
					# This may have to be changed when testing n_sources > 1
					# mimic the push! functionality
					Mixfn   = let Mixfn = Mixfn; x -> Mixfn(x) + mixfun2(x); end
				end
			end
		end
		
		if (i == 1)
			push!(fun_sum, x -> ([Mixfn(x) zeros(1, (nd-1)*numT)] - S[i,:]))
			#funF = x -> ([Mixfn(x) zeros(1, (nd-1)*numT)] - S[i,:]')
		else
			fun2 = x -> ([zeros(1, (i-1)*numT) Mixfn(x) zeros(1, (nd-i)*numT)] - S[i,:])
			push!(fun_sum, x->fun_sum[i](x)+fun2(x))
			#funF = let funF = funF; x -> funF(x) + fun2(x); end
		end
	end
	
	# Define funF as the final element in the function vector
	funF = fun_sum[end]

	GreenNMFk.log("  Calculating simulation parameters")
	
	# Defining the lower and upper boundary for the minimization
	lb = [0 0 0] # General lower boundary - [Ux Dx Dy A X Y]
	ub = [1 1 1] # General upper boundary - [Ux Dx Dy A X Y]
	
	# This loop is on the number of sources we investigate
	# We need limits for all sources (ampl and coord)
	for jj = 1:number_of_sources
		lb = [lb 0 -aa -aa] # General lower boundary [Ux Dx Dy A X Y]
		ub = [ub 1.5 aa aa] # General upper boundary [Ux Dx Dy A X Y]
	end
	
	# The norm of the observational matrix / vector
	AA = 0
	SS = 0
	
	for i = 1:nd
		SS = S[i,:].^2
		AA = AA + sum(SS)
	end
	
	real_num = 0
	sol_real = []
	normF_real = []
	normF1 = []
	sol_all = []
	j_all = 0
	DidItGoBack = 0
	
	cutNum = 5
	
	GreenNMFk.log("  Running NMF calculations")
	while ((real_num < cutNum) && (j_all < 10 * Nsim))
		#options = optimset("MaxFunEvals", 3000) # Set maximum times to evaluate function at 3000
		initCON = zeros(Nsim, 3 * number_of_sources + 3)
		
		GreenNMFk.log("  -> Calculating initial conditions")
		for k = 1:Nsim
			x_init = [rand(1) 2*aa*(0.5 - rand(1,2))]
			
			for d = 1:number_of_sources
				# The size is 3*number_of_sources+3 for the IC
				x_init = [x_init rand() 2*aa*(0.5 - rand(1,2))]
			end
			
			initCON[k,:] = x_init
		end
		
		# Iterates the NMF runs
		#TODO: implement parallel
		#parfor k = 1:Nsim
		GreenNMFk.log("  -> Calculating the non-linear least squares fit")
		
		model_NMF = JuMP.Model(solver=Ipopt.IpoptSolver())
		@JuMP.variable(model_NMF, lb[i] <= x[i=1:6] <= ub[i])
		JuMP.register(model_NMF, :funF, 1, funF, autodiff=true)
		
		@JuMP.NLobjective(model_NMF, Min, funF(x[1]))
		status = JuMP.solve(model_NMF)
		println("x = ", getValue(x))
		
		GreenNMFk.log(status)
		
		for k = 1:Nsim
			println(initCON[k,:])
#			sol[k,:], normF[k] = lsqnonlin(funF,initCON[k,:],lb,ub,options)
		end
		
		normF_abs = normF
		normF = sqrt(normF./AA).*100
		
		normCut = 0.1

		# Find the indices of normF where the element is less than normCut
		index_real = find(normF[normF .< normCut]) #find(normF < normCut)
		
		real_num = real_num + length(index_real)
		normF_real = [normF_real; normF[index_real]]
		sol_real = [sol_real; sol[index_real,:]]
		normF1 = [normF1; normF]
		sol_all = [sol_all; sol]
		
		j_all = j_all + Nsim
		
		if ((j_all == 10 * Nsim) && (real_num > 0) && (real_num < cutNum))
			DidItGoBack += 1
			j_all = 0
		end
	end
		
	if (real_num < cutNum)
		sol_real = sol_all
		normF_real = normF1
		Qyes = 1
	end
		
	# Save a JLD file with initial condition variables
	if save_output
		outfile = "Results_$(nd)det_$(number_of_sources)sources.jld"
		
		GreenNMFk.log("Saving results to $(working_dir)/$(outfile)")
		JLD.save(joinpath(working_dir, outfile), "sol", sol, "normF", normF, "S", S, "lb", lb, "ub", ub, "AA", AA, "sol_real", sol_real, "normF_real", normF_real, "normF_abs", normF_abs, "DidItGoBack", DidItGoBack, "Qyes", Qyes)
	end
	
	return sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all, normF_abs, Qyes
end
		