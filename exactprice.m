% Codes are written by: Agniv Bandyopadhyay
%% When r=alpha=0, this function gives the deflated price
function[price] = exactprice(xarr,T,sig)
tau = sig^2*T;
z=0.5*sqrt(tau)-log(1+abs(xarr))/sqrt(tau);
price = max(xarr,0)+0.5*(z*sqrt(tau)+1).*normcdf(z)-0.5*(1+abs(xarr)).*normcdf(z-sqrt(tau))+0.5*sqrt(tau)*normpdf(z);
end
