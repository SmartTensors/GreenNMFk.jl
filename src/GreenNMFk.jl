# Module file for Green-NMFk
module GreenNMFk

	import DocumentFunction
	import NMFk # Parent library
	import Gadfly # Plotting library
	import JLD # Used for file IO
	import Base # Used for system paths
	import DataFrames # Helpful in Gadfly plotting
	import Clustering # Silhouette clustering
	import Mads # Has Levenberg–Marquardt solver
	import Distributions

	# IO level: 0 = no IO, 1 = some IO, 2 = verbose
	io_level = 2

	# Save output from simulation? true/false
	save_output = true

	# Show Gadfly plots?
	show_plots = false

	# Set to true to test against Matlab output
	matlab_tests = true

	const nmfkdir = splitdir(splitdir(Base.source_path())[1])[1]

	# Edit to change desired output directory
	working_dir = joinpath(nmfkdir, "Results")

	# If directory doesn't exist, create it
	if !isdir(working_dir)
		Base.mkdir(working_dir)
	end

	if matlab_tests
		import MAT
		matlab_dir = joinpath(nmfkdir, "matlab", "Results")
	end

	include("GreenSource.jl")
	include("GreenInitial.jl")
	include("GreenCalculationsNMF.jl")
	include("GreenCompRes.jl")
	include("GreenLudmilCluster.jl")
	include("GreenClustering.jl")
	#include("src/tests.jl")
	include("GreenExecute.jl")

end