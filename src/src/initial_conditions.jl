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

	# Iterate through number of detectors and build observation matrix
	for i=1:nd
		for d=1:ns
			if (d == 1)
				min_sum = source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3]) + 0 * randn(size(time))
			else
				min_sum += source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
			end
			W[i,d] = sum(min_sum)
		end
		S[i, :] = [zeros(1, (i-1)*numT) min_sum' zeros(1, (nd-i)*numT)]
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
		
		if GreenNMFk.io_level > 0
			println("Saving results to $(working_dir)/$(outfile)")
		end
		
		JLD.save(joinpath(working_dir, outfile), "XF", XF, "S", S)
	end
	
	return S, XF, W
end