function[price] = fdm_kspc(xmin,xmax,tmin,tmax,A,alpha,Nsp,Nt)
%% Grids and initializations
spgrids = linspace(xmin,xmax,Nsp+1);
dx = (xmax-xmin)/Nsp; dt = (tmax-tmin)/Nt; lambda = dt/dx^2;
tgrids = linspace(tmin,tmax,Nt+1);
price = zeros(size(tgrids,2),size(spgrids,2));
Axarr = A(spgrids(2:end-1));

%% Initial Conditions(level j=1)
price(1,:) = max(spgrids,0);

%% Crank Nicholson for(level j=2)
M=diag([1+0.5*alpha*dt+lambda*Axarr,3])+diag(-0.5*lambda*Axarr,1)+diag([-0.5*lambda*Axarr(2:end),-4],-1);
M(end,end-2)=1;
b=0.5*lambda*(price(1,1:end-2)+price(1,3:end)).*Axarr+(1-0.5*alpha*dt-lambda*Axarr).*price(1,2:end-1);
b=[b,2*dx*exp(-alpha*tgrids(1,2))];
b=b';
price(2,2:end) = (M\b)';

%% Three time level scheme(j>=3)
for j=3:Nt+1
    M=diag([-(1.5+alpha*dt+2*lambda*Axarr),3])+diag(lambda*Axarr,1)+diag([lambda*Axarr(2:end),-4],-1);
    M(end,end-2)=1;
    b = -2*price(j-1,2:end-1)+0.5*price(j-2,2:end-1);
    b= [b,2*dx*exp(-alpha*tgrids(1,j))];
    b=b';
    price(j,2:end) = (M\b)';
end

end