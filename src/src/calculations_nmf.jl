# TODO: Parallel solver
#		See this link: https://juliaeconomics.com/2014/06/18/parallel-processing-in-julia-bootstrapping-the-mle/

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

	# Number of independent variables to solve
	nvar = 3 + number_of_sources * 3

	sol = Array{Float64}(Nsim, 3 * number_of_sources + 3) # Solution matrix
	normF = Array{Float64}(Nsim, 1)
	normF_abs = Array{Float64}(Nsim, 1)
	
	normCut = 0
	Qyes = 0

	# Initialize arrays
	sol_real = []
	normF_real = []
	normF1 = []
	sol_all = []
	local real_num = j_all = DidItGoBack = 0

	cutNum = 5

	# Define the function that will be minimized
	# funF = ∑ᵢ(MixFn[i]- ∑ⱼ(Sources[i,j]))²
	GreenNMFk.log("  Define the function to be minimized")
	function nl_func(a...)
		x = collect(a)
		min_sum = 0
		fun_sum = 0
		for i=1:nd
			if number_of_sources == 1
				min_sum += source(time, x[4], x[5:6], xD[i,:], x[1], x[2], t0, x[3]) # replace with true variable names
			else
				for d=1:number_of_sources
					if (d == 1)
						min_sum += source(time, x[4], x[5:6], xD[i,:], x[1], x[2], t0, x[3])
					else
						min_sum += source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
					end
				end
			end

			if (i == 1)
				fun_sum += [min_sum; zeros((nd-1)*numT)] - S[i,:]
			else
				fun_sum += [zeros((i-1)*numT); min_sum; zeros((nd-i)*numT)] - S[i,:]
			end
		end
		return sum(fun_sum.^2)
	end

	GreenNMFk.log("  Calculating simulation parameters")
	# Defining the lower and upper boundary for the minimization
	lb = [0 0 0] # replace with true variable names
	ub = [1 1 1] # replace with true variable names

	# This loop is on the number of sources we investigate
	# We need limits for all sources (ampl and coord)
	for jj = 1:number_of_sources
		lb = [lb 0 -aa -aa] # replace with true variable names
		ub = [ub 1.5 aa aa] # replace with true variable names
	end

	# The norm of the observational matrix / vector
	AA = 0
	SS = 0

	for i = 1:nd
		SS = S[i,:].^2
		AA = AA + sum(SS)
	end

	GreenNMFk.log("  Running NMF calculations")
	while ((real_num < cutNum) && (j_all < 10 * Nsim))
		#options = optimset("MaxFunEvals", 3000) # Set maximum times to evaluate function at 3000
		initCON = Array{Float64}(Nsim, 3 * number_of_sources + 3)

		GreenNMFk.log("  -> Calculating initial conditions")
		for k = 1:Nsim
			x_init = [rand(1) 2*aa*(0.5 - rand(1,2))]

			for d = 1:number_of_sources
				# The size is 3*number_of_sources+3 for the IC
				x_init = [x_init rand() 2*aa*(0.5 - rand(1,2))]
			end

			initCON[k,:] = x_init
		end
		
		GreenNMFk.log("  -> Running solver")
		# Run solver once for every number of simulations

		for k = 1:Nsim
			
			local solutions
			result = 0 # For try/catch: 0 if failure, 1 for success
			max_iters = 3000 # Maximum # of iterations for solver - decreases on failure
			
			# Try/catch NL solver
			nl_solver(iters) = try
				model_NMF = JuMP.Model(solver=Ipopt.IpoptSolver(print_level=3, max_iter=iters))
		    	JuMP.register(model_NMF, :nl_func, nvar, nl_func, autodiff=true)
		    	@JuMP.variable(model_NMF, lb[i] <= x[i=1:nvar] <= ub[i], start=initCON[k,i])
		    	JuMP.setNLobjective(model_NMF, :Min, Expr(:call, :nl_func, [x[i] for i=1:nvar]...))
		
		    	JuMP.solve(model_NMF)
				
				result = 1
		    	return [JuMP.getvalue(x[i]) for i=1:nvar]
			catch y
				if isa(y, DomainError)
					result = 0
					return 0
				end
			end
		
			# While solver fails...
			while result == 0
				GreenNMFk.log("  --> Run $(k)/$(Nsim): Trying with $(max_iters) iterations")
				solutions = nl_solver(max_iters)
				max_iters = max_iters - 200
				
				if max_iters < 500
					GreenNMFk.log("Iterations too small: stopping")
					exit(1)
				end
			end
			
			# If solver was successful
			if result == 1
				sol[k,:] = solutions
				GreenNMFk.log("\n  --> Run $(k) solutions: $(sol[k,:])")
			end
		end

		normF_abs = normF
		normF = sqrt(normF./AA).*100

		normCut = 0.1

		# Find the indices of normF where the element is less than normCut
		index_real = find(normF .< normCut)

		real_num = real_num + length(index_real)
		normF_real = [normF_real; normF[index_real]]
		
		if sol_real == []
			sol_real = sol[index_real,:]
			sol_all = sol
		else
			sol_real = [sol_real; sol[index_real,:]]
			sol_all = [sol_all; sol]
		end
		
		normF1 = [normF1; normF]

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
