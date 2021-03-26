% Pls make sure xarr is of the form -L:(L/N):L
function[gamma] = exact_gamma(xarr,T,sig)
gamma = zeros(size(xarr));
tau = sig^2*T;
N = size(xarr,2);
n = (N-1)/2;
xn = xarr(1,1:n);
xp = xarr(1,n+2:N);
zarrp = 0.5*sqrt(tau)-log(1+xp)/sqrt(tau);
zarrn = 0.5*sqrt(tau)-log(1-xn)/sqrt(tau);
dzarrp = -1./(sqrt(tau)*(1+xp));
dzarrn = 1./(sqrt(tau)*(1-xn));
d2zarrp = 1./(sqrt(tau)*(1+xp).^2);
d2zarrn = 1./(sqrt(tau)*(1-xn).^2);
gamma(1,n+1) = 0.5*normcdf(0.5*sqrt(tau))+normpdf(0.5*sqrt(tau))/sqrt(tau);
gamma(1,n+2:N) = 0.5*sqrt(tau)*(dzarrp.^2).*normpdf(zarrp)+0.5*(zarrp*sqrt(tau)+1).*(dzarrp.^2).*(-zarrp.*normpdf(zarrp))+0.5*(zarrp*sqrt(tau)+1).*normpdf(zarrp).*d2zarrp+0.5*sqrt(tau)*d2zarrp.*normcdf(zarrp)+0.5*sqrt(tau)*(dzarrp.^2).*normpdf(zarrp)-0.5*normpdf(zarrp-sqrt(tau)).*dzarrp-0.5*(1+xp).*(-(zarrp-sqrt(tau)).*normpdf(zarrp-sqrt(tau))).*(dzarrp.^2)-0.5*(1+xp).*normpdf(zarrp-sqrt(tau)).*d2zarrp-0.5*normpdf(zarrp-sqrt(tau)).*dzarrp+0.5*sqrt(tau)*(zarrp.^2-1).*normpdf(zarrp).*(dzarrp.^2)+0.5*sqrt(tau)*(-zarrp.*normpdf(zarrp)).*d2zarrp;
gamma(1,1:n) = 0.5*sqrt(tau)*normpdf(zarrn).*(dzarrn.^2)+0.5*(zarrn*sqrt(tau)+1).*(-zarrn.*normpdf(zarrn)).*(dzarrn.^2)+0.5*(zarrn*sqrt(tau)+1).*normpdf(zarrn).*d2zarrn+0.5*sqrt(tau)*d2zarrn.*normcdf(zarrn)+0.5*sqrt(tau)*(dzarrn.^2).*normpdf(zarrn)+0.5*normpdf(zarrn-sqrt(tau)).*dzarrn-0.5*(1-xn).*(-(zarrn-sqrt(tau)).*normpdf(zarrn-sqrt(tau))).*(dzarrn.^2)-0.5*(1-xn).*normpdf(zarrn-sqrt(tau)).*d2zarrn+0.5*normpdf(zarrn-sqrt(tau)).*dzarrn+0.5*sqrt(tau)*(zarrn.^2-1).*normpdf(zarrn).*(dzarrn.^2)+0.5*sqrt(tau)*(-zarrn.*normpdf(zarrn)).*d2zarrn;
end