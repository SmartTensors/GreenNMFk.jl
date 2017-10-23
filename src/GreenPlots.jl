function plot(Wf, xD)
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

function plot(RECON,SILL_AVG,max_number_of_sources)
    X = Array{Float64}(1:1:max_number_of_sources)
	plt = Gadfly.plot(Gadfly.layer(x=X, y=RECON[:], Gadfly.Geom.line), Gadfly.layer(x=X, y=SILL_AVG[:], Gadfly.Geom.line))
	
    display(plt)
end