include("../GreenNMFk.jl")
import GreenNMFk
import Base.Test

function get_mat_value(infile::String, key::String)
	matlab_file = MAT.matopen(joinpath(GreenNMFk.matlab_dir, infile))
	mat_keys = MAT.names(matlab_file)
	return MAT.read(matlab_file, key)
end

# Initial_example9d4sdiff
Nsim = 100
t0 = -10
As = [0.5; 0.5; 0.5; 0.5]
D = [0.005 0.00125]
u = 0.05
numT = 80
noise = 0E-3
xD = [0.0 0.0; -0.5 -0.5; -0.5 0.5; 0.5 0.5;
		0.5 -0.5;0.0 0.5; 0.0 -0.5; -0.5 0.0; 0.5 0.0]
Xn = [[-0.3; -0.4] [0.4; -0.3] [-0.1; 0.25] [-0.3; 0.65]];

S, sol, normF, sol_real, normF_real, lb, ub, Qyes = GreenNMFk.execute(Nsim, t0, As, D, u, numT, noise, xD, Xn)

#----- Error in solution --------------#
#   This code block is only present for debugging
ML_sol = get_mat_value("Results_9det_4sources.mat","sol")

avg_ML = Array{Float64}(size(sol)[2])
avg_sol = Array{Float64}(size(sol)[2])

for i=1:size(sol)[2]
	avg_ML[i] = round(mean(ML_sol[:,i]),5)
	avg_sol[i] = round(mean(sol[:,i]),5)
end

println("Debugging:\n","-"^40)
println(" Column-based average of:")
println("    Julia solution:\n    $(avg_sol)\n")
println("    Matlab solution:\n    $(avg_ML)\n")
println(" Delta between solutions:")
println("    ML - Julia = $(avg_ML-avg_sol)")
#--------------------------------------

@Base.Test.testset "GreenNMFk" begin
	@Base.Test.test isapprox(sum(S),231.992505, atol=1e-5)
	@Base.Test.test Qyes == 0
	@Base.Test.test isapprox(mean(normF),16.15177,atol=1e-5)
	@Base.Test.test isapprox(mean(sol_real),0.132844,atol=1e-4)

	@Base.Test.test lb == [1e-6 1e-6 1e-6 1e-6 -1 -1 1e-6 -1 -1 1e-6 -1 -1 1e-6 -1 -1]
	@Base.Test.test ub == [1.0 1.0 1.0 1.5 1.0 1.0 1.5 1.0 1.0 1.5 1.0 1.0 1.5 1.0 1.0]
	
	ML_average_sol = [0.02508,0.06989,0.10817,0.93109,-0.13471,0.14874,0.96308,-0.14126,0.15853,0.95115,-0.10073,0.10416,0.98761,-0.10117,0.04129]
	for i=1:size(sol)[2]
		@Base.Test.test isapprox(mean(sol[:,i]),ML_average_sol[i],atol=1e-3)
	end
	
	if GreenNMFk.matlab_tests == true
		ML_S = get_mat_value("Results_9det_4sources.mat","S")
		ML_sol = get_mat_value("Results_9det_4sources.mat","sol")
		ML_normF = get_mat_value("Results_9det_4sources.mat","normF")
		ML_lb = get_mat_value("Results_9det_4sources.mat","lb")
		ML_ub = get_mat_value("Results_9det_4sources.mat","ub")
		ML_sol_real = get_mat_value("Results_9det_4sources.mat","sol_real")
		ML_Qyes = get_mat_value("Results_9det_4sources.mat","Qyes")

		@Base.Test.test isapprox(S,ML_S,atol=1e-3)
		@Base.Test.test isapprox(sol,ML_sol,atol=1e-3)
		@Base.Test.test isapprox(normF,ML_normF,atol=1e-2)
		@Base.Test.test isapprox(lb,ML_lb,atol=1e-6)
		@Base.Test.test isapprox(ub,ML_ub,atol=1e-6)
		@Base.Test.test isapprox(sol_real,ML_sol_real,atol=1e-3)
		@Base.Test.test isapprox(Qyes,ML_Qyes,atol=1e-3)
	end
end