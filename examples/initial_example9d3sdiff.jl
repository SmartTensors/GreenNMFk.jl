import GreenNMFk

#Start stopwatch timer
tic()

# Initialize simulation variables
max_number_of_sources = 9	# Number of sources
Nsim = 100 					# Number of simulations
numT = 80					# Number of time points
time = collect(linspace(0, 20, numT))

# Initialize equation parameters
u = 0.05 # (km/yr) - flow speed
D = [0.005 0.00125] # (km²/yr) - diffusion coefficient
t0 = -10 # Initial time of sources
noise = 0E-3 # Noise strength
As = [0.5; 0.5; 0.5; 0.5] # Amplitudes of real sources
ns = length(As) # 'Real' number of sources

aa = 1 # Multiple for initial random conditions
Xn = Array{Float64}(2,length(As)) # Initialize Xn -> source coordinates
Xn = [[-0.3; -0.4] [0.4; -0.3] [ -0.1; 0.25] [ -0.3; 0.65]];
Xs = ones(length(As), 3)

# Ordering matrix of the sources: [A X Y]
for k = 1:size(Xs,1)
	Xs[k,:] = [As[k] Xn[1,k] Xn[2,k]]
end

xD = Array{Float64}(9, 2) # do not use zeros in Julia in not needed
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

S, XF = GreenNMFk.initial_conditions(As,Xs,xD,D,t0,u,numT,noise,time)

number_of_sources = 1
GreenNMFk.log("\nNumber of sources : $(number_of_sources)")
#sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all = GreenNMFk.calculations_nmf_v02(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT)

##########################################################################################
prefix = "Results_9det_1sources"
matlab_file = MAT.matopen(joinpath(GreenNMFk.matlab_dir, prefix * ".mat"))
sol = MAT.read(matlab_file, "sol")
normF = MAT.read(matlab_file, "normF")[:]
lb = MAT.read(matlab_file, "lb")
ub = MAT.read(matlab_file, "ub")
AA = MAT.read(matlab_file, "AA")
sol_real = MAT.read(matlab_file, "sol_real")
normF_real = MAT.read(matlab_file, "normF_real")
#normF1 = MAT.read(matlab_file, "normF1")
#sol_all = MAT.read(matlab_file, "sol_all")
##########################################################################################

GreenNMFk.log("Size of sol : $(size(sol))")

yy = Base.quantile(normF,0.25)
reconstr1 = mean(normF[normF .< yy])

ind = find(normF[normF .< yy])
sol1 = sol[ind,:]

avg_sol = mean(sol1)
solution = avg_sol
mean_savg = 1
number_of_clust_sim = 0

# Save a JLD file with simulation results
if GreenNMFk.save_output
	outfile = "Solution_$(nd)det_$(number_of_sources)sources.jld"

	GreenNMFk.log("-> Saving results to $(GreenNMFk.working_dir)/$(outfile)")
	JLD.save(joinpath(GreenNMFk.working_dir, outfile), "solution", solution, "reconstr1", reconstr1, "mean_savg", mean_savg, "number_of_clust_sum", number_of_clust_sim)
end

# Generate empty RECON and SILL_AVG vectors
RECON = zeros(max_number_of_sources,1)
SILL_AVG = zeros(max_number_of_sources,1)

RECON[1] = reconstr1
SILL_AVG[1] = 1

number_of_sources = 2

##########################################################################################
prefix = "Results_9det_2sources"
matlab_file = MAT.matopen(joinpath(GreenNMFk.matlab_dir, prefix * ".mat"))
sol = MAT.read(matlab_file, "sol")
normF = MAT.read(matlab_file, "normF")[:]
lb = MAT.read(matlab_file, "lb")
ub = MAT.read(matlab_file, "ub")
AA = MAT.read(matlab_file, "AA")
sol_real = MAT.read(matlab_file, "sol_real")
normF_real = MAT.read(matlab_file, "normF_real")[:]
Qyes = MAT.read(matlab_file, "Qyes")

#=
solution, vect_index, cent, reconstr, mean_savg, number_of_clust_sim = GreenNMFk.clustering_the_solutions(number_of_sources, nd, sol_real, normF_real, Qyes)
##########################################################################################

for jj = 2:max_number_of_sources

	sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all, normF_abs, Qyes = GreenNMFk.calculations_nmf_v02(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT)
	solution, vect_index, cent, reconstr, mean_savg, number_of_clust_sim = GreenNMFk.clustering_the_solutions(number_of_sources, nd, sol_real, normF_real, Qyes)

	RECON[number_of_sources] = reconstr
	SILL_AVG[number_of_sources] = mean_savg

	number_of_sources = number_of_sources + 1
	
end

#close all

RECON = RECON/Nsim

# Save a JLD file with simulation results
if save_output
	outfile = "All_$(nd)det_$(length(As))sources.jld"

	GreenNMFk.log("Saving results to $(working_dir)/$(outfile)")
	JLD.save(joinpath(working_dir, outfile), "RECON", RECON, "SILL_AVG", SILL_AVG)
end
=#
##########################################################################################
prefix = "Solution_$(nd)det_$(max_number_of_sources)sources"
matlab_file = MAT.matopen(joinpath(GreenNMFk.matlab_dir, prefix * ".mat"))
cent = MAT.read(matlab_file, "Cent")
solution = MAT.read(matlab_file, "Solution")

prefix = "All_$(nd)det_$(length(As))sources"
matlab_file = MAT.matopen(joinpath(GreenNMFk.matlab_dir, prefix * ".mat"))
RECON = MAT.read(matlab_file, "RECON")
SILL_AVG = MAT.read(matlab_file, "SILL_AVG")
##########################################################################################

if GreenNMFk.show_plots == true
	
	X = Array{Float64}(1:1:max_number_of_sources)
	plt = Gadfly.plot(Gadfly.layer(x=X, y=RECON[:], Gadfly.Geom.line), Gadfly.layer(x=X, y=SILL_AVG[:], Gadfly.Geom.line))
	
	display(plt)
	
end

Sf, Comp, Dr, Det, Wf = GreenNMFk.comp_res(cent, xD, solution, t0, numT, noise, S)

# End stopwatch timer

toc()