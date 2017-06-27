# Returns source
# Inputs:
#		t	- time vector from t_initial to t_final
#       fs  - source strength
#		xs	- sources positions
#		xd	- detector positions
#		Dx	- diffusion coefficient x
#       Dy  - diffusion coefficient y
#		t0	- initial time of sources
#		u	- flow speed [km/year]
function source(t::Vector, fs::Number, xs::Vector, xd::Vector, Dx::Number, Dy::Number, t0::Number, u::Number)
	
	# Cheap solution - if Dx, Dy < 0 then Dx, Dy = 0
	(Dx < 0) && (Dx = 0)
	(Dy < 0) && (Dy = 0)
	
	# The solver will occasionally pass negative values for Dx/Dy. Try/catch them.
	# ERROR: DomainError:
	# sqrt will only return a complex result if called with a complex argument. Try sqrt(complex(x)).
	comp_coefficients() = try
		coeff1 = (fs./(4*pi*sqrt(Dx*Dy)*(t-t0)))
		coeff2 = exp(-(((xd[1]-(xs[1] + u*(t-t0))).^2)./(4*Dx*(t-t0))))
		coeff3 = exp(-((xd[2] - xs[2]).^2)./(4*Dy*(t-t0)))
		
		return coeff1.*coeff2.*coeff3
	catch y
		println("Error in calculation:")
		println("  Dx = $(Dx)\n  Dy = $(Dy)")
		return nothing
	end
	
	comp_coefficients()
end