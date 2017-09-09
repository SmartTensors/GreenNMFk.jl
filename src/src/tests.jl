# Using the module MAT, which can read in Matlab files,
# test generated Julia files against known good Matlab files
function test_results_mat(infile::String; ml_dir::String=matlab_dir,wk_dir::String=working_dir)
	prefix = split(infile, ".")[end-1]

	# Open corresponding Matlab and Julia files
	matlab_file = MAT.matopen(joinpath(ml_dir, prefix * ".mat"))
	julia_file = JLD.load(joinpath(wk_dir, infile))

	# Get variable names as dict
	mat_keys = MAT.names(matlab_file)

	@Base.Test.testset "JLD/MAT Test for: $(prefix)" begin
		# Iterate over keys and test matching ones
		for key in mat_keys
			mat_result = MAT.read(matlab_file, key)
			jld_result = julia_file[key]

			rtype = string(typeof(mat_result))
			if contains(rtype, "Array")

				# If dimensions don't match, throw an error
				# May be better just to let the test fail
				if size(mat_result) != size(jld_result)
					println("ERROR: Dimension mismatch for $(key) in $(prefix) ($(size(mat_result)) != $(size(jld_result)))")
				else
					@show mean(sum(mat_result .- jld_result))
					@Base.Test.test isapprox(mean(sum(mat_result .- jld_result)), 0.0, atol=1e-2)
				end

			elseif contains(rtype, "Float")
				@Base.Test.test isapprox(mat_result, jld_result, atol=1e-6)
			elseif contains(rtype, "Int")
				@Base.Test.test mat_results == jld_result
			end
		end
	end
end