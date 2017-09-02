# import Gadfly
import JuMP
import Ipopt

# srand(1)
d = collect(linspace(0,3,100))
# y = exp(-1.3*d) + 0.05*randn(size(d))
y = readdlm("y.dat")'
# Gadfly.plot(x=d,y=y)

function source(t::Vector, fs::Number, xs::Vector, xd::Vector, Dx::Number, Dy::Number, t0::Number, u::Number)
	# Cheap solution - if Dx || Dy < 0 then Dx || Dy = 0
	(Dx < 0) && (Dx = 0)
	(Dy < 0) && (Dy = 0)

	# The solver will occasionally pass negative values for Dx/Dy. Try/catch them.
	# ERROR: DomainError:
	# sqrt will only return a complex result if called with a complex argument. Try sqrt(complex(x)).
	comp_coefficients() = try
		coeff1 = (fs./(4*pi*sqrt(Dx*Dy)*(t-t0)))
		coeff2 = exp(-(((xd[1]-(xs[1] + u*(t-t0))).^2)./(4*Dx*(t-t0))))
		coeff3 = exp(-((xd[2] - xs[2]).^2)./(4*Dy*(t-t0)))
		#@show coeff1
		#@show coeff2
		#@show coeff3

		#@show coeff1.*coeff2.*coeff3
		return coeff1.*coeff2.*coeff3
	catch errmsg
		Base.showerror(Base.STDERR, errmsg)
		println("Error in calculation:")
		println("  Dx = $(Dx)\n  Dy = $(Dy)")
		return nothing
	end

	comp_coefficients()
end

function nl_func(a...)
	x = collect(a)
	min_sum = 0
	fun_sum = 0
	i = 1
	# @show x[4], x[5:6], xD[i,:], x[1], x[2], t0, x[3]
	# @show time
	# @show source(time, x[4], x[5:6], xD[i,:], x[1], x[2], t0, x[3])
	for i=1:nd
		if number_of_sources == 1
			min_sum += source(time, x[4], x[5:6], xD[i,:], x[1], x[2], t0, x[3]) # replace with true variable names
		else
			for d=1:number_of_sources
				if (d == 1)
					min_sum += source(time, x[4], x[5:6], xD[i,:], x[1], x[2], t0, x[3])
				else
					min_sum += source(time, x[d*3+1], x[d*3+2:d*3+3], xD[i,:], x[1], x[2], t0, x[3])
				end
			end
		end

		if (i == 1)
			fun_sum += [min_sum; zeros((nd-1)*numT)] - S[i,:]
		else
			fun_sum += [zeros((i-1)*numT); min_sum; zeros((nd-i)*numT)] - S[i,:]
		end
	end
	# @show fun_sum
	return sum(fun_sum.^2)
end


nvar = 1
m = JuMP.Model(solver=Ipopt.IpoptSolver())
JuMP.register(m, :nl_func, nvar, nl_func, autodiff=true)
@JuMP.variable(m, x[i=1:nvar], start=4)
# @eval @JuMP.NLobjective(m, Min, $(Expr(:call, :nl_func, [Expr(:ref, :x, i) for i=1:nvar]...)))
JuMP.setNLobjective(m, :Min, Expr(:call, :nl_func, [x[i] for i=1:nvar]...))
JuMP.solve(m)
JuMP.getvalue(x)