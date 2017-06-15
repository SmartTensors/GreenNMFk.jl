# Module file for Green-NMFk
module GreenNMFk

import JLD # Used for file IO
import Base # Used for system paths
using JuMP # Used for nonlinear modeling

# Save output from simulation? true/false
save_output = true

# Logging function for debugging
# Currently print, can be moved to file IO
function log(instream)
	println(instream)
end

# Set up IO filepaths and directories
# Procedure inspired by Mads.jl/src/MadsIO.jl

# Get directory of where the source is running
source_path = Base.source_path()
if typeof(source_path) == Void
	root_dir = "."
else
	root_dir = dirname(source_path)
end

# Edit to change desired output directory
working_dir = joinpath(root_dir, "..", "Results")

# If directory doesn't exist, create it
if !isdir(working_dir)
	Base.mkdir(working_dir)
end

include("src/source.jl")
include("src/initial_conditions.jl")
include("src/calculations_nmf.jl")
#include("src/clustering_the_solutions.jl")

end