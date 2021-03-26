% Codes are written by: Agniv Bandyopadhyay
% Pls make sure xarr is of the form -L:(L/N):L
function[delta] = exact_delta(xarr,T,sig)
delta = zeros(size(xarr));
tau = sig^2*T;
N = size(xarr,2);
n = (N-1)/2;
xn = xarr(1,1:n);
xp = xarr(1,n+2:N);
zarrp = 0.5*sqrt(tau)-log(1+xp)/sqrt(tau);
zarrn = 0.5*sqrt(tau)-log(1-xn)/sqrt(tau);
dzarrp = -1./(sqrt(tau)*(1+xp));
dzarrn = 1./(sqrt(tau)*(1-xn));
delta(1,n+2:N) = 1-0.5*normcdf(zarrp-sqrt(tau))+dzarrp.*(0.5*(zarrp*sqrt(tau)+1).*normpdf(zarrp)-0.5*(1+xp).*normpdf(zarrp-sqrt(tau))+0.5*sqrt(tau)*(-zarrp.*normpdf(zarrp))+0.5*sqrt(tau)*normcdf(zarrp));
delta(1,n+1) = 0.5;
delta(1,1:n) = 0.5*normcdf(zarrn-sqrt(tau))+dzarrn.*(0.5*(zarrn*sqrt(tau)+1).*normpdf(zarrn)+0.5*sqrt(tau)*normcdf(zarrn)-0.5*(1-xn).*normpdf(zarrn-sqrt(tau))+0.5*sqrt(tau)*(-zarrn.*normpdf(zarrn)));
end
