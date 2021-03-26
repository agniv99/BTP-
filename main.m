%% Initialization 
% Codes are written by: Agniv Bandyopadhyay
close all; clear all; clc;
S_init = 100; sig=0.3;T=1;Nsp=1000;Nt=50; 
xmin = -1; xmax = +1;
dx = (xmax-xmin)/Nsp; dt = T/Nt;
A = @(x) 0.5*sig^2*(1+abs(x)).^2;

%% Computing Prices Numerically
deflated_price_kspc = fdm_kspc(xmin,xmax,0,1,A,0,Nsp,Nt);
deflated_price_cn = fdm_cn(xmin,xmax,0,1,A,0,Nsp,Nt);
price_kspc = S_init*deflated_price_kspc(end,:);
price_cn = S_init*deflated_price_cn(end,:);
price_exact = S_init*exactprice(linspace(xmin,xmax,1001),1,sig);

%% Error while evaluating price
err_kspc = abs(price_kspc-price_exact);
err_cn = abs(price_cn-price_exact);
gain = linspace(-30,30,301);
gain2 = linspace(-1,1,11);
%% Exact Value of Greeks
delta_exact = S_init*exact_delta(gain/100,T,sig);
gamma_exact = S_init*exact_gamma(gain/100,T,sig);
theta_exact = S_init*exact_theta(gain/100,T,sig);

%% Numerically Computed Value of Greeks
delta_kspc = (price_kspc(1,352:652)-price_kspc(1,350:650))/(2*dx); 
delta_cn = (price_cn(1,352:652)-price_cn(1,350:650))/(2*dx);
gamma_kspc = (price_kspc(1,352:652)-2*price_kspc(1,351:651)+price_kspc(1,350:650))/(dx^2);
gamma_cn = (price_cn(1,352:652)-2*price_cn(1,351:651)+price_cn(1,350:650))/(dx^2);
theta_kspc = -S_init*(deflated_price_kspc(end,351:651)-deflated_price_kspc(end-1,351:651))/dt;
theta_cn = -S_init*(deflated_price_cn(end,351:651)-deflated_price_cn(end-1,351:651))/dt;

%% Error while computing Greeks
err_delta_kspc = abs(delta_kspc-delta_exact);
err_delta_cn = abs(delta_cn-delta_exact);
err_gamma_kspc = abs(gamma_kspc-gamma_exact);
err_gamma_cn = abs(gamma_cn-gamma_exact);
err_theta_kspc = abs(theta_kspc-theta_exact);
err_theta_cn = abs(theta_cn-theta_exact);

%% Plotting the results
% Price
pic=figure(1);
subplot(1,2,1);
plot(gain2,price_kspc(1,496:506));
hold on;
plot(gain2,price_exact(1,496:506),'r--');
hold off;
legend(["numerical","exact"]);
title("Three time level Scheme");
subplot(1,2,2);
plot(gain2,price_cn(1,496:506));
hold on;
plot(gain2,price_exact(1,496:506),'r--');
hold off;
legend(["numerical","exact"]);
title("Crank-Nicholson Scheme");
%saveas(pic,"plots/price","epsc");
%saveas(pic,"plots/price.png");
pic=figure(2);
subplot(1,2,1);
plot(gain,err_kspc(1,351:651));
title("Three time level");
subplot(1,2,2);
plot(gain,err_cn(1,351:651));
title("Crank-Nicholson");
%saveas(pic,"plots/err_price","epsc");
%saveas(pic,"plots/err_price.png");
%% Error while evaluating greeks
% Delta
pic=figure(3);
subplot(1,2,1);
plot(gain,delta_kspc);
hold on;
plot(gain,delta_exact,'r--');
hold off;
legend(["numerical","exact"]);
title("Three time level Scheme");
subplot(1,2,2);
plot(gain,delta_cn);
hold on;
plot(gain,delta_exact,'r--');
hold off;
legend(["numerical","exact"]);
title("Crank-Nicholson Scheme");
%saveas(pic,"plots/delta","epsc");
%saveas(pic,"plots/delta.png");
% Gamma
pic=figure(4);
subplot(1,2,1);
plot(gain,gamma_kspc);
hold on;
plot(gain,gamma_exact,'r--');
hold off;
legend(["numerical","exact"]);
title("Three time level Scheme");
subplot(1,2,2);
plot(gain,gamma_cn);
hold on;
plot(gain,gamma_exact,'r--');
hold off;
legend(["numerical","exact"]);
title("Crank-Nicholson Scheme");
%saveas(pic,"plots/gamma","epsc");
%saveas(pic,"plots/gamma.png");
% Theta
pic=figure(5);
subplot(1,2,1);
plot(gain,theta_kspc);
hold on;
plot(gain,theta_exact,'r--');
hold off;
legend(["numerical","exact"]);
title("Three time level Scheme")
subplot(1,2,2);
plot(gain,theta_cn);
hold on;
plot(gain,theta_exact,'r--');
hold off;
legend(["numerical","exact"]);
title("Crank-Nicholson Scheme");
%saveas(pic,"plots/theta","epsc");
%saveas(pic,"plots/theta.png");
