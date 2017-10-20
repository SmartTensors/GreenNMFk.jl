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

# Set up IO filepaths and directories
# Procedure inspired by Mads.jl/src/MadsIO.jl

# Get directory of where the source is running
source_path = Base.source_dir()
println(source_path)
if typeof(source_path) == Void
	root_dir = "."
else
	root_dir = dirname(source_path)
end

# Edit to change desired output directory
working_dir = joinpath(root_dir, "Results")

# If directory doesn't exist, create it
if !isdir(working_dir)
	Base.mkdir(working_dir)
end

if matlab_tests
	import MAT
	matlab_dir = joinpath(root_dir, "matlab", "Results")
end

include("src/source.jl")
include("src/initial_conditions.jl")
include("src/calculations_nmf.jl")
include("src/comp_res.jl")
include("src/ludmil_cluster.jl")
include("src/clustering_the_solutions.jl")
include("src/tests.jl")
include("src/NMFkExecute.jl")

end