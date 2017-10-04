include("../GreenNMFk.jl")
import GreenNMFk
import Base.Test

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

S, sol, normF, sol_real = GreenNMFk.execute(Nsim, t0, As, D, u, numT, noise, xD, Xn)

@Base.Test.testset "simulation" begin
	@Base.Test.test isapprox(sum(S), 231.992505, atol=1e-5)
    @Base.Test.test isapprox(mean(S),0.161846,atol=1e-3)
    @Base.Test.test isapprox(mean(sol),0.13041667,atol=1e-4)
    #@Base.Test.test isapprox(mean(normF),0.7337,atol=1e-2)
end