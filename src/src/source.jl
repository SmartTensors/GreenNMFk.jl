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
function source(t::Vector, fs::Vector, xs::Matrix, xd::Matrix, Dx::Number, Dy::Number, t0::Number, u::Number)
	coeff1 = (fs./(4*pi*sqrt(D[1]*D[2])*(t-t0)))
	coeff2 = exp(-(((xd[1]-(xs[1] + u*(t-t0))).^2)./(4*Dx*(t-t0))))
	coeff3 = exp(-((xd[2] - xs[2]).^2)./(4*Dy*(t-t0)))
	return coeff1.*coeff2.*coeff3
end