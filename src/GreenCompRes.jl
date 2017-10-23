function comp_res(cent, xD, solution, t0, numT, noise, S)
	
	time = collect(linspace(0, 20, numT))
	
	if (size(solution,1) == 1)
		As = solution[4]
		Xs = solution[4:6]
		D = solution[1:2]
		u = solution[3]
	else
		As = cent[:,1]
		Xs = cent
		D = solution[1,end-2:end-1]
		u = solution[1,end]
	end
	
	Sf, XFf, Wf = initial_conditions_2(As, Xs, xD, D, t0, u, numT, noise, time)
	comp = sum((Sf.^2 - S.^2), 2)
	
	# Define empty S & Sf matrices
	Det = Array{Float64}(size(xD, 1), size(S[1,1:80],1))
	Dr = Array{Float64}(size(xD, 1), size(Sf[1:80],1))
	
	for i = 1:size(xD, 1)
		a = (i-1) * 80 + 1
		b = 80 * i
		Det[i,:] = S[i,a:b]
		Dr[i,:] = Sf[a:b]
	end
	
	if show_plots == true
		# Generate a data frame to produce a stacked bar plot
		# Dependent on Wf being square matrix
		data = DataFrames.DataFrame()

		s = size(Wf,1)
		data[:Mixes_Normalized] = ["Wf_$(Int64(floor(i/s)))" for i=s:s*s+s-1] # 'Labels' for columns
		data[:Index] = [(i%s+1) for i=0:s*s-1] # Ranges from 1:N; N number of times (i.e. 1,2,3,4, 1,2,3,4, 1,2...)
		data[:Value] = vec(Wf) # Condense matrix into vector

		plt = Gadfly.plot(data, x="Index", y="Value", color="Mixes_Normalized", Gadfly.Geom.bar(position=:stack))
		display(plt)
	
		# Build second plot
		# This is a stacked plot, showing an overlay of a pulse and signal for each source 
		s = size(xD,1)
		local plots = Array{Any}(s, 1)
		X = collect(linspace(0,s*s-1,s*s-1))
	
		# Generate function for returning some plot, which will then be stacked over prior generated plots
		plt_func(Y1, Y2) = Gadfly.plot(Gadfly.layer(x=X, y=Y1, Gadfly.Geom.point), Gadfly.layer(x=X, y=Y2, Gadfly.Geom.line))
		[plots[i] = plt_func(Det[i,:], Dr[i,:]) for i=1:s] # Add plots to array
		
		stack = Gadfly.vstack([plots[i] for i=1:s]) # Decompose plots into a form vstack can read
		display(stack)
	end
	
	return Sf, comp, Dr, Det, Wf
end