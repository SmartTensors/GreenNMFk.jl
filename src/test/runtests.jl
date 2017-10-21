include("../GreenNMFk.jl")
import GreenNMFk
import Base.Test

# Setting this false runs a larger (but much longer!) test suite
simple_tests = true

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
	@Base.Test.test isapprox(sum(S),231.9925, atol=1e-3)
	#@Base.Test.test isapprox(sum(normF),551.0656,atol=1e-5)
	#@Base.Test.test isapprox(sum(sol_real),2.0759e3,atol=100)

	@Base.Test.test lb == [1e-6 1e-6 1e-6 1e-6 -1.0 -1.0 1e-6 -1.0 -1.0]
	@Base.Test.test ub == [1.0 1.0 1.0 1.5 1.0 1.0 1.5 1.0 1.0]

	# Julia solution = [0.00548 0.01408 0.04827 1.22777 -0.34733 -0.44911 0.79602 0.50358 -0.42129]
	true_sol = [0.0054 0.0019 0.0522 0.8902 -0.2413 -0.3211 1.2067 -0.2580 0.7398] # From Matlab
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