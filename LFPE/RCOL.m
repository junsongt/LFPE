function RCOL(estar,gam)
% =========================================
% The reactive collision frequency R(x) from
% Can. J. Phys. paper to the matrix rep of the LFPE - D'.D
% =========================================
%
%Get the I_n integrals
% gam is the mass ratio m1/m2 - mass fractions M1 and M2
%
pt=[0:.1:10];
nmax=length(pt);
xm2=1/(1+gam); xm1=gam/(1+gam);
% estar = e*/kT
% Get G(-sqrt(x)) and G(sqrt(x))
for nn=1:nmax
    nn;
    pt(nn);
t0=(sqrt(estar)+sqrt(xm2*pt(nn)))/sqrt(xm1);
Int0=sqrt(pi)*erfc(t0)/2;
Int1=0.5*(exp(-t0^2));
Int2=0.5*(exp(-t0^2)*t0+Int0);
%pause
Gp(nn)=sqrt(xm1)*(xm1*Int2-2*sqrt(xm1*xm2*pt(nn))*Int1+(xm2*pt(nn)-estar)*Int0);
Gm(nn)=sqrt(xm1)*(xm1*Int2+2*sqrt(xm1*xm2*pt(nn))*Int1+(xm2*pt(nn)-estar)*Int0);
end
% ==========
% Now get R(x)
% ==========
for nn=1:nmax
%    R(nn)= exp(-pt(nn))*(Gm(nn)-Gp(nn))/(sqrt(xm1*xm2*pt(nn)));
    R(nn)= (Gm(nn)-Gp(nn))/(sqrt(xm1*xm2));
    fprintf('%12.4e,%15.5e,%15.5e,%15.5e\n',pt(nn),Gm(nn),Gp(nn),R(nn))
    pause
end
Fac=1000*pt.*exp(-pt);
plot(pt,R.*Fac,'-k')
disp('Reactive Collision Frequency')
pause