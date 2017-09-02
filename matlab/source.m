function s = source(t, fs, xs, xd, D, t0, u)
s = (fs./(4*pi*sqrt(D(1)*D(2))*(t-t0))).*exp(-(((xd(1)-(xs(1) + u*(t-t0))).^2)./(4*D(1)*(t-t0)))).*exp(-((xd(2)- xs(2)).^2)./(4*D(2)*(t-t0)));
end

%s = (xs(1)./(4*pi*D*(t-t0))).*exp(-(((xd(1)-(xs(2) + u(1)*(t-t0))).^2 + (xd(2)-(xs(3) + u(2)*(t-t0))).^2)./(4*D*(t-t0))));
%exp(-((xd(2)- xs(3) - u(2)*(t-t0)).^2)./(4*D*(t-t0)));