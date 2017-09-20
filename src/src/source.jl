"""
Returns observation matrix S. The observation matrix is a sparse matrix
with a diagonal populated by mixes.
$(DocumentFunction.documentfunction(source;
argtext=Dict("t"=>"time vector from t_initial to t_final",
			"fs"=>"source strength",
			"xs"=>"sources positions",
			"xd"=>"detector positions",
			"Dx"=>"diffusion coefficient x",
			"Dy"=>"diffusion coefficient y",
			"t0"=>"initial time of sources",
            "u"=>"flow speed in km/year")))
Returns:
- mix
"""

function source(t::Vector, fs::Number, xs::Vector, xd::Vector, Dx::Number, Dy::Number, t0::Number, u::Number)

	coeff1 = (fs./(4*pi*sqrt(Dx*Dy)*(t-t0)))
	coeff2 = exp(-(((xd[1] - (xs[1] + u*(t-t0))).^2)./(4*Dx*(t-t0))))
	coeff3 = exp(-((xd[2] - xs[2]).^2)./(4*Dy*(t-t0)))

	return coeff1.*coeff2.*coeff3
end