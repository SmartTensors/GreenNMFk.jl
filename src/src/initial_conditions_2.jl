#                        ---80---    ---80---  ---80---   --80--     --80---
#                           1           2         3         4          5
# The structure of S is [first mix   0,0...0,   0,0...0,  0,0...0    0,0...0]
#                       [0, 0...0   second mix  0,0...0,  0,0...0    0,0...0]
#                       [0, 0...0    0,0...0   third mix  0,0...0    0,0...0]
#                       [0, 0...0    0,0...0    0,0...0, fourth mix  0,0...0]
#                       [0, 0...0    0,0...0    0,0...0,  0,0...0  fifth mix]

function initial_conditions_2(As, Xs, xD, D, t0, u, numT, noise, time)
	# Calculation of the observation matrix
	nd = size(xD, 1)
	S = zeros(nd, nd*numT)
	W = zeros(nd, length(As))
	
	# Below, we are ordering the signals consecutively in a single vector
	#========= Construction of the mixes at the detectors =========#
	for i = 1:nd # Loop over number of detectors (nd)
		
		for d = 1:length(As) # Loop over the count of sources
			if (d == 1) # If first source...
				Mix = source(time, Xs[d,:], xD[i,:], D, t0, u)
				W[i,d] = sum(Mix)
			else # If second or higher source...
				Mix = Mix + source(time, Xs[d,:], xD[i,:], D, t0, u)
				W[i,d] = sum(Mix)
			end
		end
		
		S[i,:] = [ zeros(1, (i-1)*numT) Mix zeros(1, nd-i)*numT ]

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
	XF = reshape(S', 1, size(S,2) * nd)
	
	# The original source constants, [Ai Xi Yi]
	xtrue = reshape(Xs', 1, size(Xs, 2) * size(Xs, 1))
	
	return S, XF, W
end