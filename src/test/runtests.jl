include("../GreenNMFk.jl")
import GreenNMFk
import Base.Test

# Setting this false runs a larger (but much longer!) test suite
simple_tests = false

# Epsilon for calculated vs. actual solution values
sol_tolerance = 1e-3

# columnwise_approx(matrix, solution vector, tolerance)
# Takes column-wise average of matrix, then compares to
# that index in real solution
function columnwise_approx(sol,real,tol)
	mean_sol = Array{Float64}(size(sol)[2])
	N = size(sol)[2]

	# Iterate over each column
	for i=1:N
		mean_sol[i] = mean(sol[:,i])

		# Return false if difference > tolerance
		if abs(mean_sol[i] - real[i]) >= tol
			return false
		end
	end

	return true
end

# Values sourced from Initial_example9d4sdiff.m
Nsim = 100
t0 = -10
As = [0.5; 0.5; 0.5; 0.5]
D = [0.005 0.00125]
u = 0.05
numT = 80
noise = 0E-3
xD = [0.0 0.0; -0.5 -0.5; -0.5 0.5; 0.5 0.5;
		0.5 -0.5;0.0 0.5; 0.0 -0.5; -0.5 0.0; 0.5 0.0]
Xn = [[-0.3; -0.4] [0.4; -0.3] [-0.1; 0.25] [-0.3; 0.65]]

# Initial conditions for simple test
nsources = 2
x_init = [0.00990707 0.0819074 0.109797 0.410209 -0.200658 -0.429437 0.532161 0.470071 -0.323975]

# Run Green-NMFk solver
S, sol, normF, sol_real, normF_real, lb, ub, Qyes = GreenNMFk.execute(Nsim, t0, As, D, u, numT, noise, xD, Xn, number_of_sources=nsources, x_init=x_init)

@Base.Test.testset "Green-NMFk" begin
	#@Base.Test.test isapprox(sum(S),231.992505, atol=1e-5)
	#@Base.Test.test isapprox(mean(normF),16.15177,atol=1e-5)
	#@Base.Test.test isapprox(mean(sol_real),0.132844,atol=1e-4)

	@Base.Test.test lb == [1e-6 1e-6 1e-6 1e-6 -1.0 -1.0 1e-6 -1.0 -1.0]
	@Base.Test.test ub == [1.0 1.0 1.0 1.5 1.0 1.0 1.5 1.0 1.0]

	true_sol = [0.00548 0.01408 0.04827 1.22777 -0.34733 -0.44911 0.79602 0.50358 -0.42129]
	@Base.Test.test columnwise_approx(sol,true_sol,sol_tolerance) == true

	# This portion of the test runs more complicated tests, but they take much longer to converge
	if simple_tests == false
		nsources = 3
		x_init = [0.0194544 0.0879254 0.0435591 0.506107 -0.279446 -0.371166 0.555605 0.429977 -0.398611 0.519732 -0.191911 0.270702]
		S,sol,_ = GreenNMFk.execute(Nsim, t0, As, D, u, numT, noise, xD, Xn, number_of_sources=nsources, x_init=x_init)
		true_sol = [0.00413 0.03965 0.04242 0.77465 -0.29941 -0.36911 0.61061 0.43244 -0.41551 0.78064 -0.24139 0.27156]
		@Base.Test.test columnwise_approx(sol,true_sol,sol_tolerance) == true

		nsources = 4
		x_init = [0.0758843 0.0568675 0.048459 0.561233 -0.387407 -0.468881 0.460469 0.363911 -0.258558 0.416207 -0.147628 0.304434 0.543241 -0.383654 0.720799]
		S,sol,_ = GreenNMFk.execute(Nsim, t0, As, D, u, numT, noise, xD, Xn, number_of_sources=nsources, x_init=x_init)
		true_sol = [0.00629 0.02891 0.05117 0.75646 -0.44025 -0.50155 0.32589 0.38881 -0.33513 0.40995 -0.23239 0.36519 0.75455 -0.39463 0.7222]
		@Base.Test.test columnwise_approx(sol,true_sol,sol_tolerance) == true
	end
end














#=



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
=#