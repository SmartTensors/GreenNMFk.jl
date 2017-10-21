
"""
Calculation of the solution with j sources.
Minimizes the function sum_i((MixFn[i]) - sum_j(Sources[i,j]))**2

$(DocumentFunction.documentfunction(calculations_nmf;
argtext=Dict("number_of_sources"=>"Number of sources",
			"nd"=>"Number of detectors",
			"Nsim"=>"Number of simulations to run",
			"aa"=>"Boundary condition coefficient",
			"xD"=>"Detector positions",
			"t0"=>"Initial time of sources",
			"time"=>"Time vector",
			"S"=>"Observation matrix",
			"numT"=>"Number of time points",
			"x_true"=>"The true solutions vector"),
keytext=Dict("x_init"=>"Initial conditions for solver (of length 3*number_of_sources+3)",
			 "tol"=>"Solver tolerance")))

Returns:
- sol : solutions matrix
- normF : squared 2-norm of the residual at x
- lb : lower bounds
- ub : upper bounds
- AA
- sol_real
- normF_real
- normF1
- sol_all
- normF_abs
- Qyes
"""

function calculations_nmf(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT, x_true; tol::Float64=1e-3, epsilon=0.1, x_init=nothing)
	
	srand()

	#--- Function definitions ---------------------------
	# rg(Vector{Float64}) 
	# Populates a vector containing random values within [0,1]
	function rg(x::Vector{Float64}=Vector{Float64}(0))
		if length(x) == 0
			x_init = rand(3)'
			for d = 1:number_of_sources
				x_init = [x_init rand()*1.5 aa * (rand()-0.5) aa * (rand()-0.5)]
			end
		else
			x_init = x + (rand(length(x)) * 0.001)
		end
		return x_init
	end

	# nl_func(x...)
	# Defines a function to be solved with a solver
	function nl_func(a...)
		x = collect(a)

		fun_sum2 = 0
		min_sum = zeros(numberoftimes)

		for i=1:nd
			for d=1:number_of_sources
				min_sum += source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
			end
			
			fun_sum2 += sum((min_sum .- S[i,(i-1)*numT+1:i*numT]).^2)
			fun_sum2 += sum(S[i,1:i*numT:(i-1)*numT].^2) + sum(S[i,i*numT+1:nd*numT].^2)
		end

		return fun_sum2
	end

	# nl_func_mads(x)
	# Defines a function to be solved with MADS solver
	function nl_func_mads(x)
		fun_sum = zeros(nd*numT)
		min_sum = zeros(numberoftimes)

		for i=1:nd
			for d=1:number_of_sources
				if (d == 1)
					min_sum = source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
				else
					min_sum += source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
				end
			end
			
			fun_sum = [zeros((i-1)*numT); min_sum; zeros((nd-i)*numT)] .- S[i,:]
		end

		return fun_sum
	end

	#--- Variable definitions ----------------------------
	# Number of independent variables to solve
	nvar = 3 + number_of_sources * 3
	normCut = 0.1
	Qyes = 0
	cutNum = 5

	numberoftimes = length(time)

	sol = Array{Float64}(Nsim, 3 * number_of_sources + 3) # Solution matrix
	normF = Array{Float64}(Nsim, 1)
	normF_abs = Array{Float64}(Nsim, 1)

	# Initialize arrays
	sol_real::Array{Float64} = []
	normF_real::Array{Float64} = []
	normF1::Array{Float64} = []
	sol_all::Array{Float64} = []
	
	local real_num = j_all = DidItGoBack = 0
	
	# The norm of the observational matrix / vector
	AA = 0
	SS = 0

	for i = 1:nd
		SS = S[i,:].^2
		AA = AA + sum(SS)
	end

	# Defining the lower and upper boundary for the minimization
	lb = [1e-6 1e-6 1e-6] # replace with true variable names
	ub = [1 1 1] # replace with true variable names
	pl = [true true true]

	# This loop is on the number of sources we investigate
	# We need limits for all sources (ampl and coord)
	for jj = 1:number_of_sources
		lb = [lb 1e-6 -aa -aa] # replace with true variable names
		ub = [ub 1.5 aa aa] # replace with true variable names
		pl = [pl true false false]
	end

	ub = convert(Array{Float64}, ub)
	lb = convert(Array{Float64}, lb)

	if x_init == nothing
		x_init = similar(lb)#similar(x_true)
		for iii = 1:size(vec(ub))[1]
			#x_init[iii] = rand(Distributions.Uniform(lb[iii],ub[iii]))/10.0
			tmp_val = x_true[iii]
			x_init[iii] = rand(Distributions.Uniform(tmp_val-epsilon,tmp_val+epsilon))
		end
	end

	#--- Minimization solver loop -----------------------------
	for k = 1:Nsim
		local solutions
		of = Inf

		if GreenNMFk.io_level > 0
			print_with_color(:red,"\rSolving: $(k)/$(Nsim) runs..."," "^15,"\r")
		end

		#for iii = 1:size(vec(ub))[1]
		#	x_init[iii] = rand(Distributions.Uniform(lb[iii],ub[iii]))/10.0
		#end

		# Repeat the minimization until the minimum is greater than the tolerance
		while of > tol
			solutions, r = Mads.minimize(nl_func_mads, vec(x_init); upperbounds=vec(ub), lowerbounds=vec(lb), logtransformed=vec(pl), tolX=1e-6, tolG=1e-6, tolOF=1e-3)
			of = r.minimum
		end
		
		# Populate the solutions matrix with the solutions vector
		sol[k,:] = solutions
	end


	#--- Generate solution-dependant variables -----------------
	# normF is defined as the squared 2-norm of the residual at x: sum(fun(x).^2)
	for i=1:size(sol)[2]
		normF[i] = float(sum(sol[:,i].^2))
	end
	
	normF_abs = normF
	normF = sqrt(normF./AA).*100

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

	if (real_num < cutNum)
		sol_real = sol_all
		normF_real = normF1
		Qyes = 1
	end

	# Save a JLD file with initial condition variables
	if save_output
		outfile = "Results_$(nd)det_$(number_of_sources)sources.jld"
		
		if GreenNMFk.io_level > 0
			info("Saving results to $(working_dir)/$(outfile)")
		end
		
		JLD.save(joinpath(working_dir, outfile), "sol", sol, "normF", normF, "S", S, "lb", lb, "ub", ub, "AA", AA, "sol_real", sol_real, "normF_real", normF_real, "normF_abs", normF_abs, "DidItGoBack", DidItGoBack, "Qyes", Qyes)
	end

	return sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all, normF_abs, Qyes
end