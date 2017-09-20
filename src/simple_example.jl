include("GreenNMFk.jl")
import GreenNMFk

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

GreenNMFk.execute(Nsim, t0, As, D, u, numT, noise, xD, Xn)

#=
# -- Set up variables -------------------------
number_of_sources = 1
Nsim = 100
aa = 1
t0 = -10
numT = 80	
As = [0.5; 0.5; 0.5; 0.5]
D = [0.005 0.00125]
u = 0.05
numT = 80
noise = 0E-3
time = collect(linspace(0, 20, numT))
xD = [0.0 0.0; -0.5 -0.5; -0.5 0.5; 0.5 0.5;
	  0.5 -0.5;0.0 0.5; 0.0 -0.5; -0.5 0.0; 0.5 0.0]
nd = size(xD,1)
ns = length(As) # 'Real' number of sources

Xs = Array{Float64}(length(As),3)
Xn = [[-0.3; -0.4] [0.4; -0.3] [-0.1; 0.25] [-0.3; 0.65]];
for k = 1:size(Xs,1)
	Xs[k,:] = [As[k] Xn[1,k] Xn[2,k]]
end

x_true = [D[1], D[2], u]
for k = 1:size(Xs,1)
	x_true = [x_true..., As[k], Xn[1,k], Xn[2,k]]
end

number_of_sources = ns
Nsim = 100

srand(2015)








# -- Run core functions -----------------------
tic()
S, XF = GreenNMFk.initial_conditions(As,Xs,xD,D,t0,u,numT,noise,time)

print("S = $(mean(S))")
S = GreenNMFk.create_problem(x_true, nd, numT, ns, xD, t0, time)
print("\n\n\n")
print("S = $(mean(S))")

quit()
tc1 = toc()
tic()
#sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all, normF_abs, Qyes = GreenNMFk.calculations_nmf_v02(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT, x_true)
tc2 = toc()
tic()
sol, normF, lb, ub, AA, sol_real, normF_real, normF1, sol_all, normF_abs, Qyes = GreenNMFk.calculations_nmf_v02(number_of_sources, nd, Nsim, aa, xD, t0, time, S, numT, x_true)
#solution, vect_index, cent, reconstr, mean_savg, number_of_clust_sim = GreenNMFk.clustering_the_solutions(number_of_sources+1, nd, sol_real, normF_real, Qyes)
tc3 = toc()

println(sol)

# -- Verify that functions are correct --------
println("Testing function initial conditions: ")
#outfile = "test_results/initial_conditions.jld"
JLD.save("test_results/initial_conditions.jld", "S", S, "XF", XF)
JLD.save("test_results/calculations_nmf.jld","sol",sol,"normF",normF,"lb",lb,"ub",ub,
		"AA",AA,"sol_real",sol_real,"normF_real",normF_real,"normF1",normF1,"sol_all",
		sol_all,"normF_abs",normF_abs,"Qyes",Qyes)

println(normF)


GreenNMFk.test_results_mat("initial_conditions.jld",ml_dir="test_results/",wk_dir="test_results/")
GreenNMFk.test_results_mat("calculations_nmf.jld",ml_dir="test_results/",wk_dir="test_results/")

println("Time to compute initial conditions: $(tc1) s")
println("Time to compute calculations_nmf:   $(tc2) s")
println("Time to compute clustering:         $(tc3) s")
=#