% Codes are written by: Agniv Bandyopadhyay
% Pls make sure xarr is of the form -L:(L/N):L
function[theta] = exact_theta(xarr,T,sig)
tau = sig^2*T;
zarr = 0.5*sqrt(tau)-log(1+abs(xarr))/sqrt(tau);
theta = -0.25*sig^2*normcdf(zarr)-(0.25*sig^2/sqrt(tau))*(log(1+abs(xarr))/tau+1.5).*normpdf(zarr)+(0.25*sig^2/sqrt(tau))*(1+abs(xarr)).*(log(1+abs(xarr))/tau-0.5).*normpdf(zarr-sqrt(tau));
end
