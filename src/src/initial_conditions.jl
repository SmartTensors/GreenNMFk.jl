#---- Construction of the mixes at the detectors ----#

# With this function, we are ordering the signals
#  one after another in one long vector

#                        ---80---    ---80---  ---80---   --80--     --80---
#                           1           2         3         4          5
# The structure of S is [first mix   0,0...0,   0,0...0,  0,0...0    0,0...0]
#                       [0, 0...0   second mix  0,0...0,  0,0...0    0,0...0]
#                       [0, 0...0    0,0...0   third mix  0,0...0    0,0...0]
#                       [0, 0...0    0,0...0    0,0...0, fourth mix  0,0...0]
#                       [0, 0...0    0,0...0    0,0...0,  0,0...0  fifth mix]

# Initial conditions - calculation of observation matrix
# Input:
#	As    - Amplitudes of real sources
#	Xs    - Sources matrix
#	xD    - Detector positions
#	D     - Diffusion coefficient [km^2/year]
#	t0    - Initial time of sources
#	u     - Flow speed [km/year]
#	numT  - Number of time points
#	noise - Noise strength
#	time  - Time vector from t_initial to t_final
#
# Returns:
#	S	  - Observation matrix
#	XF	  - Observation matrix condensed to long vector

function initial_conditions(As,Xs,xD,D,t0,u,numT,noise,time)

	nd = size(xD,1) # The number of the detectors
	S  = Array{Float64}(nd, nd * numT) # The observation matrix S with unknown # of sources

	GreenNMFk.log("\nInitial conditions set:")
	GreenNMFk.log("-----------------------------------------")
	GreenNMFk.log("  Sources amplitude (As)     = $(As)")
	GreenNMFk.log("  Sources matrix length (Xs) = $(length(Xs))")
	GreenNMFk.log("  Detector positions (xD)    = $(xD)")
	GreenNMFk.log("  Diffusion coeff. (D)       = $(D) km²/year")
	GreenNMFk.log("  Initial time (t0)          = $(t0)")
	GreenNMFk.log("  Flow speed (u)             = $(u) km/year")
	GreenNMFk.log("  Noise (noise)              = $(noise)")
	GreenNMFk.log("  Number of detectors (nd)   = $(nd)")

	GreenNMFk.log("\nCalculating initial conditions...")
	GreenNMFk.log("-----------------------------------------")

	# Iterate through number of detectors
	for d=1:nd
		if (length(As) == 1)
			Mix = GreenNMFk.source(time, Xs[1,:], xD[d,:], D, t0, u) + noise * randn(size(time))
			S[d, :] = [zeros(1, (d - 1) * numT) Mix zeros(1, (nd - d) * numT)]
		else
			for i=1:length(As)
				if (i == 1)
					Mix = GreenNMFk.source(time, Xs[i,:], xD[d,:], D, t0, u) + noise * randn(size(time))
				else
					Mix = Mix + GreenNMFk.source(time, Xs[i,:], xD[d,:], D, t0, u)
				end
			end

			# Generate row in S
			# [0, 0...0    0,0...0   Mix'  0,0...0    0,0...0]
			S[d, :] = [zeros(1, (d-1)*numT)  Mix'  zeros(1, (nd-d)*numT)]
			GreenNMFk.log("  Filling column $(d) of $(nd) in observation matrix")
		end
	end

	GreenNMFk.log("\n  S dimensions: $(size(S))")

	# Condense S' and Xs' matrices into single vector
	XF = reshape(S', 1, size(S,2) * nd) # Long vector of length nd * numT
	xtrue = reshape(Xs', 1, size(Xs, 2) * size(Xs, 1))

	# Save a JLD file with initial condition variables
	if save_output
		outfile = "x_true_$(nd)det_$(length(As))sources.jld"

		GreenNMFk.log("\n  -> Saving results to $(working_dir)/$(outfile)")
		JLD.save(joinpath(working_dir, outfile), "XF", XF, "xtrue", xtrue, "S", S, "xD", xD, "D", D, "u", u)
	end

	return S, XF
end