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
            "Xn"=>"",
            "aa"=>"(optional) boundary conditions coefficient",
            "ns"=>"(optional) number of sources")))
Returns:
- 
"""

function execute(Nsim::Integer, t0::Number, As::Vector, D::Array, u::Number, numT::Number, noise::Number, xD::Matrix, Xn::Matrix, aa=1, ns=nothing)
    srand(2017)

    time = collect(linspace(0, 20, numT))
    
    nd = size(xD,1)
    ns = length(As)
    if ns==nothing ns = length(As) end # 'Real' number of sources
    number_of_sources = ns

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
		print("\nInitializing...")
	end
	
    sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all, normF_abs, Qyes = GreenNMFk.calculations_nmf(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT, x_true)
   
   	if GreenNMFk.io_level > 0
		print("\nMean solution = \n  [ ")
		for i=1:size(sol)[2]
        	print("$(round(mean(sol[:,i]),5)) ")
		end
		print("]\n")
	end

	#print_with_color(:green,"\nClustering solutions:\n")
    #solution, vect_index, cent, reconstr, mean_savg, number_of_clust_sim = GreenNMFk.clustering_the_solutions(number_of_sources+1, nd, sol_real, vec(normF_real), Qyes)

    return S, sol, normF, sol_real, normF_real, lb, ub, Qyes
end

"""
Test Green-NMFk algorithm against verified solutions

$(DocumentFunction.documentfunction(test))
"""

function test()
    nmfkdir = pwd()
    include(joinpath(nmfkdir, "test", "runtests.jl"))
	return nothing
end