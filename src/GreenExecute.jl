"""
Execute Green-NMFk analysis for a set of sources and properties

$(DocumentFunction.documentfunction(execute;
argtext=Dict("Nsim"=>"number of simulations",
			"t0"=>"initial time",
			"As"=>"sources amplitude",
			"D"=>"diffusion coefficient",
			"u"=>"flow speed",
            "numT"=>"",
            "noise"=>"noise",
            "xD"=>"detector positions",
            "Xn"=>""),
keytext=Dict("aa"=>"Boundary conditions coefficient",
			 "ns"=>"*Real* number of sources",
			 "number_of_sources"=>"number of sources",
			 "x_init"=>"Initial conditions for solver (of length 3*number_of_sources+3)"
			)))
			
Returns:
- S
- sol
- normF
- sol_real
- normF_real
- lb
- ub
- Qyes
"""

function execute(Nsim::Integer, t0::Number, As::Vector, D::Array, u::Number, numT::Number, noise::Number, xD::Matrix, Xn::Matrix; aa=1, ns=nothing, number_of_sources=4, x_init=nothing)
    srand(2017)

    time = collect(linspace(0, 20, numT))
	solver_tol = 1e-3
	
    nd = size(xD,1)
    ns = length(As)
    if ns==nothing ns = length(As) end # 'Real' number of sources

    Xs = Array{Float64}(length(As),3)
    for k = 1:size(Xs,1)
        Xs[k,:] = [As[k] Xn[1,k] Xn[2,k]]
    end
    
    x_true = [D[1], D[2], u]
    for k = 1:size(Xs,1)
        x_true = [x_true..., As[k], Xn[1,k], Xn[2,k]]
	end
    
	if GreenNMFk.io_level > 0
		print_with_color(:green,"\nParameter space:\n")
		print("  Sources amplitude (As)     = $(As)\n")
		print("  Sources matrix length (Xs) = $(length(Xs))\n")
		print("  Detector positions (xD)    = $(xD)\n")
		print("  Diffusion coeff. (D)       = $(D) km²/year\n")
		print("  Initial time (t0)          = $(t0)\n")
		print("  Flow speed (u)             = $(u) km/year\n")
		print("  Noise (noise)              = $(noise)\n")
		print("  Number of detectors (nd)   = $(nd)\n\n")
	end
 
    S, XF, W = GreenNMFk.initialize(x_true, nd, numT, ns, xD, t0, time, As)

	if GreenNMFk.io_level > 0
		print_with_color(:green,"\nRunning Levenberg–Marquardt solver:\n")
		print("  Iterations                 = $(Nsim)\n")
		print("  Sources                    = $(number_of_sources)\n")
		print("  Detectors                  = $(nd)\n")
		print("  Tolerance                  = $(solver_tol)\n")
		print("\nInitializing...")
	end
	
    sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all, normF_abs, Qyes = GreenNMFk.calculations_nmf(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT, x_true, x_init=x_init, tol=solver_tol)
   
	print_with_color(:green,"\nClustering solutions:\n")
    solution, vect_index, cent, reconstr, mean_savg, number_of_clust_sim = GreenNMFk.clustering_the_solutions(number_of_sources+1, nd, sol_real, vec(normF_real), Qyes)

	print_with_color(:green,"\nCompRes:\n")
	Sf, Comp, Dr, Det, Wf = GreenNMFk.comp_res(cent, xD, solution, t0, numT, noise, S)

    return S, sol, normF, sol_real, normF_real, lb, ub, Qyes
end


"""
Performs unit testing on Green-NMFk functions and solver integrity
	
	$(DocumentFunction.documentfunction(test))
"""

function test()
    include(joinpath(GreenNMFk.nmfkdir, "test", "runtests.jl"))
	return nothing
end