% res = LFPE_Rx(5, 1, 3, 1);
% x = res.quadpt;
% y = res.Rx;
% yfit = res.Rpoly;
% plot(x, y);
% hold on;
% plot(x, yfit);


rm = 1;
Es = [2,5,8,10];
vs = [0.158, 0.0551, 0.0132, 0.00508];

for i = 1:length(Es)
    E = Es(i);
    res = NEQ3(rm, E, 1, 1, 2);
    % disp(res.A);
    mtx = res.mtx2;
    eigvalues = eig(mtx);
    lam0 = min(abs(eigvalues));
    % disp(eigvalues);
    disp(lam0);
    % v = vs(i);
    % k = lam0 / (1 - v);
    % disp(k);
    % fac = k / (exp(-E));

end



% xm1 = rm/(1+rm);
% xm2 = 1- xm1;
% rtm1 = sqrt(xm1);
% rtm2 = sqrt(xm2);
% res = LFPE_Rx(30,1,E,rm);
% gp = res.Gp;
% gm = res.Gm;
% pt = res.quadpt;
% 
% 
% rtx=sqrt(pt);
% t0=(sqrt(eps)+rtm2*rtx)/rtm1;
% i0=sqrt(pi)*erfc(t0)/2;
% i1=0.5*exp(-t0.^2);
% i2=0.5*(exp(-t0.^2).*t0 + i0);
% g1=xm1*i2;
% g2=2*rtm1*rtm2*rtx.*i1;
% g3=(xm2*pt-eps).*i0;
% Gp=rtm1*(g1-g2+g3);
