reload("GreenNMFk")

function compare(infile)
	prefix = split(infile, ".")[end-1]

	# Open corresponding Matlab and Julia files
	matlab_file = MAT.matopen(joinpath(GreenNMFk.matlab_dir, prefix * ".mat"))
	julia_file = JLD.load(joinpath(GreenNMFk.working_dir, infile))

	# Get variable names as dict
	mat_keys = MAT.names(matlab_file)

	# Iterate over keys and test matching ones
	for key in mat_keys
		mat_result = MAT.read(matlab_file, key)
		jld_result = julia_file[key]

		#println("key = $(key)")
		#println("$(key)_jld = $(jld_result)")
		#println("\n")
		#println("$(key)_mat = $(mat_result)")

		rtype = string(typeof(mat_result))
		if contains(rtype, "Array")

			# If dimensions don't match, throw an error
			# May be better just to let the test fail
			if size(mat_result) != size(jld_result)
				println("ERROR: Dimension mismatch for $(key) in $(prefix) ($(size(mat_result)) != $(size(jld_result)))")
			else
				println("Δ$(key) = $(mean(sum(mat_result .- jld_result)))")
			end

		elseif contains(rtype, "Float")
			println("Δ$(key) = $(mat_result - jld_result)")
		elseif contains(rtype, "Int")
			println("Δ$(key) = $(mat_result - jld_result)")
		end
	end
end

# Initialize simulation variables
max_number_of_sources = 9	# Number of sources
Nsim = 100 					# Number of simulations
numT = 80					# Number of time points
time = collect(linspace(0, 20, numT))

# Initialize equation parameters
u = 0.05 # (km/yr) - flow speed
D = [0.005 0.00125] # (km²/yr) - diffusion coefficient
t0 = -1 # Initial time of sources
noise = 0E-3 # Noise strength
As = [0.5; 0.5; 0.5; 0.5] # Amplitudes of real sources
ns = length(As) # 'Real' number of sources

aa = 1 # Multiple for initial random conditions
Xn = [[-0.3; -0.4] [0.4; -0.3] [ -0.1; 0.25] [ -0.3; 0.65]];
Xs = ones(length(As), 3)

# Ordering matrix of the sources: [A X Y]
for k = 1:ns
	Xs[k,:] = [As[k] Xn[1,k] Xn[2,k]]
end

xD = Array{Float64}(9, 2)
xD[1,:] = [0.0   0.0]
xD[2,:] = [-0.5 -0.5]
xD[3,:] = [-0.5  0.5]
xD[4,:] = [0.5   0.5]
xD[5,:] = [0.5  -0.5]
xD[6,:] = [0.0   0.5]
xD[7,:] = [0.0  -0.5]
xD[8,:] = [-0.5  0.0]
xD[9,:] = [0.5   0.0]

nd = size(xD, 1)

##########################################################

x_true = [D[1], D[2], u]
for k = 1:size(Xs,1)
	x_true = [x_true..., As[k], Xn[1,k], Xn[2,k]]
end

Sold, XF = GreenNMFk.initial_conditions(As,Xs,xD,D,t0,u,numT,noise,time)

S = GreenNMFk.create_problem(x_true, nd, numT, ns, xD, t0, time)

# Sold is wrong
# S is correct

number_of_sources = ns
Nsim = 100

srand(2015)

sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all, normF_abs, Qyes = GreenNMFk.calculations_nmf_v02(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT, x_true)

# Test the two functions:

println("Comparing initial conditions...")
# compare("xtrue_$(nd)det_$(length(As))sources.jld")

println("\n\nComparing NMFk results...")
compare("Results_$(nd)det_$(number_of_sources)sources.jld")

# GreenNMFk.test_results_mat("xtrue_$(nd)det_$(length(As))sources.jld")
GreenNMFk.test_results_mat("Results_$(nd)det_$(number_of_sources)sources.jld")