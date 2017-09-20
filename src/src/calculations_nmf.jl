
"""
Calculation of the solution with j sources.
Minimizes the function sum_i((MixFn[i]) - sum_j(Sources[i,j]))**2

$(DocumentFunction.documentfunction(calculations_nmf;
argtext=Dict("number_of_sources"=>"",
			"nd"=>"number of detectors",
			"Nsim"=>"number of simulations to run",
			"aa"=>"boundary condition coefficient",
			"xD"=>"detector positions",
			"t0"=>"initial time of sources",
			"time"=>"time vector",
			"S"=>"observation matrix",
			"numT"=>"number of time points",
			"x_true"=>"",
            "tol"=>"(optional) solver tolerance")))
Returns:
- sol
- normF
- lb
- ub
- AA
- sol_real
- normF_real
- normF1
- sol_all
- normF_abs
- Qyes
"""

function calculations_nmf(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT, x_true; tol::Float64=1e-1)
	# Number of independent variables to solve
	nvar = 3 + number_of_sources * 3

	numberoftimes = length(time)

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

	function nl_func_mads(x)
		fun_sum = zeros(nd*numT)
		min_sum = zeros(numberoftimes)

		for i=1:nd
			for d=1:number_of_sources
				min_sum += source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
			end
			if i == 1
				fun_sum = [min_sum; zeros((nd-1)*numT)] .- S[i,:]
			else
				fun_sum += [zeros((i-1)*numT); min_sum; zeros((nd-i)*numT)] .- S[i,:]
			end
		end

		return fun_sum
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

	# The norm of the observational matrix / vector
	AA = 0
	SS = 0

	for i = 1:nd
		SS = S[i,:].^2
		AA = AA + sum(SS)
	end

	#TODO need to seperate MADS/Ipopt solvers
	#TODO need to remove result == whatever
	#TODO need to add clustering
	#TODO needs to be fixed to represent actual upper/lower bound arrays
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

	result = 0
	while result == 0

		for k = 1:Nsim
			local solutions
			of = Inf
			result = 0 # For try/catch: 0 if failure, 1 for success
			max_iters = 1000 # Maximum # of iterations for solver - decreases on failure

			# Try/catch NL solver
			# Depreciated for MADS solver
			function nl_solver(x_init)
				# model_NMF = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LD_MMA, maxeval=max_iters))
				model_NMF = JuMP.Model(solver=Ipopt.IpoptSolver(print_level=0, max_iter=max_iters))
				JuMP.register(model_NMF, :nl_func, nvar, nl_func, autodiff=true)
				@JuMP.variable(model_NMF, x[i=1:nvar], start=x_init[i])
				@JuMP.constraint(model_NMF, x[i=1:nvar] .<= ub[i=1:nvar]) # slows down if the initial guesses are
				@JuMP.constraint(model_NMF, x[i=1:nvar] .>= lb[i=1:nvar]) # the true values
				JuMP.setNLobjective(model_NMF, :Min, Expr(:call, :nl_func, [x[i] for i=1:nvar]...))

				JuMP.solve(model_NMF)
				of = JuMP.getobjectivevalue(model_NMF)
				result = 1
				return [JuMP.getvalue(x[i]) for i=1:nvar]
			end

			#tol = 1e-1
			x_init = rg()
			of = nl_func(x_init...)

			solutions = x_init

			while of > tol
				solutions, r = Mads.minimize(nl_func_mads, vec(x_true); upperbounds=vec(ub), lowerbounds=vec(lb), logtransformed=vec(pl), tolX=1e-12, tolG=1e-12, tolOF=1e-3)
				of = r.minimum
				result = 1
			end

			# If solver was successful
			if result == 1
				sol[k,:] = solutions
				#print("\n Run $(k)/$(Nsim) solutions: $(sol[k,:])")
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
