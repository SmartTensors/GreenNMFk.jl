"""
Calculates observation matrix S. The observation matrix is a sparse matrix
with a diagonal populated by mixes.
$(DocumentFunction.documentfunction(initialize;
argtext=Dict("x"=>"the original source constants, [Ai Xi Yi]",
			"nd"=>"number of detectors",
			"numT"=>"number of time points",
			"ns"=>"number of sources",
			"xD"=>"detector positions",
			"t0"=>"init. time of sources",
            "time"=>"time vector")))
Returns:
- observation matrix
- vectorized observation matrix
- normalized observation matrix
"""

function initialize(x, nd, numT, ns, xD, t0, time, As)
	S = zeros(nd, nd * numT)
	min_sum = zeros(numT)
	W = Array{Float64}(nd, length(As))

	# Populate matrix over detector count
	#=
	for i=1:nd
		for d=1:ns
			min_sum += source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
			W[i,d] = sum(min_sum)
		end
		if i == 1
			S[i,:] += [min_sum; zeros((nd-1)*numT)]
		else
			S[i,:] += [zeros((i-1)*numT); min_sum; zeros((nd-i)*numT)]
		end
	end
=#

	# Iterate through number of detectors
	for i=1:nd
		for d=1:ns
			if (d == 1)
				min_sum = source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3]) + 0 * randn(size(time))
				#Mix = GreenNMFk.source(time, Xs[i], Xs[i,2:3], xD[d,:], D[1], D[2], t0, u)
			else
				min_sum += source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
				#Mix = Mix + GreenNMFk.source(time, Xs[i], Xs[i,2:3], xD[d,:], D[1], D[2], t0, u)
			end
			W[i,d] = sum(min_sum)
		end
		S[i, :] = [zeros(1, (i-1)*numT) min_sum' zeros(1, (nd-i)*numT)]
		#S[d, :] = [zeros(1, (d-1)*numT)  Mix'  zeros(1, (nd-d)*numT)]
	end


	# Build normalization matrix
	Wo = W
	for j = 2:size(W,2)
		Wo[:,j] = W[:,j] - W[:,j-1]
	end
	
	W = Wo
	for i = 1:size(W,1)
		W[i,:] = W[i,:] ./ sum(W[i,:])
	end

	# Condense S' and Xs' matrices into single vector
	XF = reshape(S', 1, size(S,2) * nd) # Long vector of length nd * numT
	
	# Save a JLD file with initial condition variables
	if save_output
		outfile = "init_conditions_$(nd)_detectors.jld"
		JLD.save(joinpath(working_dir, outfile), "XF", XF, "S", S)
	end
	
	return S, XF, W

end



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

function initial_conditions(As::Vector, Xs::Matrix, xD::Matrix, D::Matrix, t0::Number, u::Number, numT::Number, noise::Number, time::Vector)

	nd = size(xD,1) # The number of the detectors
	S  = Array{Float64}(nd, nd * numT) # The observation matrix S with unknown # of sources

	# Iterate through number of detectors
	for d=1:nd
		if (length(As) == 1)
			Mix = GreenNMFk.source(time, Xs[1], Xs[1,2:3], xD[d,:], D[1], D[2], t0, u) + noise * randn(size(time))
			S[d, :] = [zeros(1, (d - 1) * numT) Mix zeros(1, (nd - d) * numT)]
		else
			for i=1:length(As)
				if (i == 1)
					Mix = GreenNMFk.source(time, Xs[i], Xs[i,2:3], xD[d,:], D[1], D[2], t0, u) + noise * randn(size(time))
				else
					Mix = Mix + GreenNMFk.source(time, Xs[i], Xs[i,2:3], xD[d,:], D[1], D[2], t0, u)
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


#=

function initial_conditions_2(As::Vector, Xs::Matrix, xD::Matrix, D, t0::Number, u::Number, numT::Number, noise::Number, time::Vector)
	# Calculation of the observation matrix
	nd = size(xD, 1)
	S = Array{Float64}(nd, nd*numT)
	W = Array{Float64}(nd, length(As))
	
	local Mix
	
	# Below, we are ordering the signals consecutively in a single vector

	for i = 1:nd # Loop over number of detectors (nd)
		
		for d = 1:length(As) # Loop over the count of sources
			
			if (d == 1) # If first source...
				Mix = GreenNMFk.source(time, Xs[d], Xs[d,2:3], xD[i,:], D[1], D[2], t0, u)
			else # If second or higher source...
				Mix = Mix + GreenNMFk.source(time, Xs[d], Xs[d,2:3], xD[i,:], D[1], D[2], t0, u)
			end
			
			W[i,d] = sum(Mix)
		end
		
		S[i,:] = [zeros(1, (i-1)*numT) Mix' zeros(1, (nd-i)*numT)]

	end
	
	Wo = W
	for j = 2:size(W,2)
		Wo[:,j] = W[:,j] - W[:,j-1]
	end
	
	W = Wo
	for i = 1:size(W,1)
		W[i,:] = W[i,:] ./ sum(W[i,:])
	end
	
	# This is the single 'long vector' (of length # of detectors * # of time points)
	XF = reshape(S', 1, (size(S,2) * nd))
	
	# The original source constants, [Ai Xi Yi]
	xtrue = reshape(Xs', 1, (size(Xs, 2) * size(Xs, 1)))
	
	return S, XF, W
end

=#