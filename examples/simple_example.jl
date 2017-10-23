import GreenNMFk

# Initialize variables
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

# Optional initial conditions
nsources = 2
x_init = [0.00990707 0.0819074 0.109797 0.410209 -0.200658 -0.429437 0.532161 0.470071 -0.323975]

# Run Green-NMFk solver
S, sol, normF, sol_real, normF_real, lb, ub, Qyes = GreenNMFk.execute(Nsim, t0, As, D, u, numT, noise, xD, Xn, number_of_sources=nsources, x_init=x_init)

println("GreenNMFk complete.")
println("-"^20)
println("Solution is: $(sol)")