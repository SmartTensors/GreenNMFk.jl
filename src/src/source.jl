# Returns source
# Inputs:
#		t	- time vector from t_initial to t_final
#		xs	- sources matrix
#		xd	- detector positions
#		D	- diffusion coefficient
#		t0	- initial time of sources
#		u	- flow speed [km/year]
function source(t, xs, xd, D, t0, u)
	coeff1 = (xs[1]./(4*pi*sqrt(D[1]*D[2])*(t-t0)))
	coeff2 = exp(-(((xd[1]-(xs[2] + u*(t-t0))).^2)./(4*D[1]*(t-t0))))
	coeff3 = exp(-((xd[2] - xs[3]).^2)./(4*D[2]*(t-t0)))
	
	return coeff1.*coeff2.*coeff3
end