# import Gadfly
import JuMP
import Ipopt

# srand(1)
d = collect(linspace(0,3,100))
# y = exp(-1.3*d) + 0.05*randn(size(d))
y = readdlm("y.dat")'
# Gadfly.plot(x=d,y=y)

function fun(r...)
	x = collect(r)
	sum((exp(-d*x[1])-y).^2)
end

nvar = 1
m = JuMP.Model(solver=Ipopt.IpoptSolver())
JuMP.register(m, :fun, nvar, fun, autodiff=true)
@JuMP.variable(m, x[i=1:nvar], start=4)
# @eval @JuMP.NLobjective(m, Min, $(Expr(:call, :fun, [Expr(:ref, :x, i) for i=1:nvar]...)))
JuMP.setNLobjective(m, :Min, Expr(:call, :fun, [x[i] for i=1:nvar]...))
JuMP.solve(m)
JuMP.getvalue(x)