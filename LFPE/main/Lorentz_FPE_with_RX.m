function result = Lorentz_FPE_with_RX(nmax,npoly,estar,gam)
%
% ----------- Continued with this March 6, 2022 --------------------
%
% estar is epsilon* in reactive cross section
% gam is the mass ratio m1/m2 - mass fractions M1 and M2
%
%xm1=1/(1+gam); xm2=gam/(1+gam);
% This version adds the reactive collision frequency to the
% program Lorentz-FPE.m: This code works and gives the correct
% eigenvalues
format long e
% -- nmax is the number of quad pts; npoly is to plot P_npoly(x)
output=fopen('eig_fpe.dat','wt');
% -- Load the alpha_n and beta_n recurrence coefficints
% -- Weight function w(x)=x^2exp(-x^2)
ab = ab_maxwell_p2(nmax,100,100,0,10);
a = ab(:,1);
b = ab(:,2);
% load ab2-70.dat;
% a=ab2_70(1:nmax,1); b=ab2_70(1:nmax+2,2);
rtb=sqrt(b); rtb(1)=[];
t=diag(rtb(1:nmax-1),-1)+diag(a)+diag(rtb(1:nmax-1),1);
[f,lambda]=eig(t); pt=diag(lambda);
% -- Quadrature points and weights from diagonlization of Jacobi matrix
mu_0=sqrt(pi)/4; mu_1=0.5;
wt=mu_0*f(1,:).^2;
% -- First three moments of the weight function
rtpi=sqrt(pi); xn=[0:1:nmax]; 
v=ones(1,nmax);w = zeros(1,nmax); poly=[]; polyp=[]; 
% Construct the SPEED POLYNOMIAL Basis Set
poly=[poly v./sqrt(mu_0)']';
poly=[poly (pt-a(1))/sqrt(b(2)*mu_0)];
polyp=[polyp w'];
polyp=[polyp (v./sqrt(b(2)*mu_0))'];
% Recurrence relation for the Speed Polynomials
for n=2:nmax-1
    xp=(pt-a(n)).*poly(:,n)/rtb(n)-poly(:,n-1)*rtb(n-1)/rtb(n);
    poly=[poly xp];
    yp=(pt-a(n)).*polyp(:,n)/rtb(n)-polyp(:,n-1)*rtb(n-1)/rtb(n)+...
        poly(:,n)/rtb(n);
    polyp=[polyp yp];
end
% Construct the DERIVATIVE MATRIX D:
for i=1:nmax
    for j=1:nmax
        s=sum(poly(j,:).*polyp(i,:)); %sum over P'_n(x_i)P_n(x_j)
d(i,j)=wt(j)*s;
%d(i,j)=s;
    end
end


for i=1:nmax
    for j=1:nmax
        s=0;
        for l=1:nmax
            %s=s+wt(l)*d(l,i)*d(l,j); % this will give Rayleigh
            s=s+pt(l)*wt(l)*d(l,i)*d(l,j); 
        end
        d2(i,j)=s/sqrt(wt(i)*wt(j));
    end
end


% This gives the well-known eigenvalues of the LFPE
[eigfcns,eigen]=eig(d2);
eig1=diag(eigen)';
eig2=sort(eig1,'ascend')';
for i=1:nmax
    if eig1(i)==eig2(npoly)
        indxeigf=i;
    else
    end
end
imax=nmax;

%{
% if nmax >= 20
%     imax=20;
% else
% end
for i=1:imax
   fprintf('%2i   %12.6f \n' , i, eig2(i));
end
disp('Eigenvalues - No Rx')
pause
%}


% =========================================
% Now add the reactive collision frequency R(x) from
% Can. J. Phys. paper to the matrix rep of the LFPE - D'.D
% =========================================
%
%Get the I_n integrals
% gam is the mass ratio m1/m2 - mass fractions M1 and M2
%
m2=1/(1+gam); m1=gam/(1+gam);
rtm1=sqrt(m1); rtm2=sqrt(m2);
eps=estar;


xm1=m1;xm2=m2;
% estar = e*/kT
% Get G(-sqrt(x)) and G(sqrt(x))


%{
% Transfer from roc_loc_ch3.m 
%
%Calculate gp, gm and rcol
%
h=0.01;
x=[0.00001:h:16];
nx=length(x); n=nx-1;
% Calculate gp
rtx=sqrt(x);
t0=(sqrt(eps)+rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*exp(-t0.^2).*t0+i0;
g1=m1*i2;
g2=2*rtm1*rtm2*rtx.*i1;
g3=(m2*x-eps).*i0;
gp=rtm1*(g1-g2+g3);
%Calculate gm
t0=(sqrt(eps)-rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*exp(-t0.^2).*t0+i0;
g1=m1*i2;
g2=-2*rtm1*rtm2*rtx.*i1;
g3=(m2*x-eps).*i0;
gm=rtm1*(g1-g2+g3);
rcol=(gm-gp)./(rtm1*rtm2*x);
rcol=2*x.*exp(-x).*rcol/sqrt(pi);
rate=h/3*(rcol(1)+4*sum(rcol(2:2:n))+2*sum(rcol(3:2:n-1))+rcol(n+1));
ratio=rate/(2*sqrt(2)*exp(-eps));
disp('rate-and-ratio')
pause
plot(x,1000*rcol/sqrt(2),'-k','linewidth',1.6)
disp('Line 126 - Plot of R(x)-and Eq. Rate')
axis([0 10 0 6])
set(gca,'FontSize',24)
set(gca,'Ytick',[0:2:6],'linewidth',1.6)
set(gca,'Xtick',[0:2:10],'linewidth',1.6)
%str2=num2str(npoly);
%str1={'$n = $'};
%str=strcat(str1,str2);
%text(3,.3,0,str,'Interpreter','latex','fontsize',26)
xlabel('$x$','Interpreter','latex','fontsize',30)
ylabel('$R(x)$','Interpreter','Latex','fontsize',30)
pause
%
% End of transfer
%
%plot(x,1000*rcol/sqrt(2),'-k','linewidth',1.6)
%for nn=1:nmax
%
%}


%{
% Now for the reactive collision Freq R(x) at quadrature points
%
% Calculate gp
rtx=sqrt(pt);
t0=(sqrt(eps)+rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*exp(-t0.^2).*t0+i0;
g1=m1*i2;
g2=2*rtm1*rtm2*rtx.*i1;
g3=(m2*pt-eps).*i0;
gp=rtm1*(g1-g2+g3);
%Calculate gm
t0=(sqrt(eps)-rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*exp(-t0.^2).*t0+i0;
g1=m1*i2;
g2=-2*rtm1*rtm2*rtx.*i1;
g3=(m2*pt-eps).*i0;
gm=rtm1*(g1-g2+g3);
rcol=(gm-gp)./(rtm1*rtm2*pt);
rcol=2*pt.*exp(-pt).*rcol/sqrt(pi);
%rate=h/3*(rcol(1)+4*sum(rcol(2:2:n))+2*sum(rcol(3:2:n-1))+rcol(n+1))
%rate=rate*2*sqrt(2)*exp(-eps)
%ratio=rate/(2*sqrt(2)*exp(-eps));
plot(pt,1000*rcol/sqrt(2),'-ok','linewidth',1.6)
disp('Line 173 - Plot of R(p_i) - at quad pt')
axis([0 10 0 6])
set(gca,'FontSize',24)
set(gca,'Ytick',[0:2:6],'linewidth',1.6)
set(gca,'Xtick',[0:2:10],'linewidth',1.6)
%str2=num2str(npoly);
%str1={'$n = $'};
%str=strcat(str1,str2);
%text(3,.3,0,str,'Interpreter','latex','fontsize',26)
xlabel('$x$','Interpreter','latex','fontsize',30)
ylabel('$R(x)$','Interpreter','Latex','fontsize',30)
pause
%}



%{
% These are the quadrature points
x=pt;
%
% COPY PREVIOUS LINES TO HERE - PTS are the QUAD PTD
%
% Calculate gp
rtx=sqrt(x);
t0=(sqrt(eps)+rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*exp(-t0.^2).*t0+i0;
g1=m1*i2;
g2=2*rtm1*rtm2*rtx.*i1;
g3=(m2*x-eps).*i0;
gp=rtm1*(g1-g2+g3);
%Calculate gm
t0=(sqrt(eps)-rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*exp(-t0.^2).*t0+i0;
g1=m1*i2;
g2=-2*rtm1*rtm2*rtx.*i1;
g3=(m2*x-eps).*i0;
gm=rtm1*(g1-g2+g3);
rcol=(gm-gp)./(rtm1*rtm2*x);
%rcol=2*x.*exp(-x).*rcol1/sqrt(pi);
n=length(x);
for k=1:n
rate=2*wt(k)*x(k)*exp(-x(k))*rcol(k)/sqrt(pi);
end
% Calculate k(T) with quadrature
ratio=rate/(2*sqrt(2)*exp(-estar));
disp('Line 216-Check rate cft')
pause
plot(x,1000*rcol/sqrt(2),'-k','linewidth',1.6)
disp('Line 126 - Plot of R(x)')
axis([0 10 0 6])
set(gca,'FontSize',24)
set(gca,'Ytick',[0:2:6],'linewidth',1.6)
set(gca,'Xtick',[0:2:10],'linewidth',1.6)
%str2=num2str(npoly);
%str1={'$n = $'};
%str=strcat(str1,str2);
%text(3,.3,0,str,'Interpreter','latex','fontsize',26)
xlabel('$x$','Interpreter','latex','fontsize',30)
ylabel('$R(x)$','Interpreter','Latex','fontsize',30)
pause
%}

%{
% OLD CODE BELOW
t0=(sqrt(estar)+sqrt(xm2*pt))/sqrt(xm1);
Int0=sqrt(pi)*erfc(t0)/2;
Int1=0.5*exp(-t0.^2);
Int2=0.5*exp(-t0.^2).*t0+Int0;
Gp=sqrt(xm1)*(xm1*Int2-2*sqrt(xm1*xm2*pt).*Int1+(xm2*pt-estar).*Int0);
Gm=sqrt(xm1)*(xm1*Int2-2*sqrt(xm1*xm2*pt).*Int1-(xm2*pt-estar).*Int0);
%}

% Calculate gp
rtx=sqrt(pt);
t0=(sqrt(eps)+rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*exp(-t0.^2).*t0+i0;
g1=m1*i2;
g2=2*rtm1*rtm2*rtx.*i1;
g3=(m2*pt-eps).*i0;
Gp=rtm1*(g1-g2+g3);
%Calculate gm
t0=(sqrt(eps)-rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*exp(-t0.^2).*t0+i0;
g1=m1*i2;
g2=-2*rtm1*rtm2*rtx.*i1;
g3=(m2*pt-eps).*i0;
Gm=rtm1*(g1-g2+g3);
% rcol=(gm-gp)./(rtm1*rtm2*pt);
% rcol=2*pt.*exp(-pt).*rcol/sqrt(pi);



% ==========
% Now get R(x)
% ==========

for nn=1:nmax
    FacR=2*pt(nn)*exp(-pt(nn))/sqrt(pi);
    R(nn)= FacR*(Gm(nn)-Gp(nn))/(sqrt(xm1*xm2*pt(nn)));
end


%{
% plot R(x)
disp('R(x)')
pause
plot(pt,R,'-ok')
disp('Line 166 - Reactive Collision Frequency')
pause
%}

%
% Add this to the diagonal elements of the Fokker-Planck operator with 
% the factor (de/dr)^2
%


% subtracting Average of R(x) when discretiing the operator
avgRx = 0;
for i = 1:nmax
    avgRx = avgRx + wt(i)*R(i);
end
avgRx = avgRx / (pt(nmax) - pt(1));

% for i=1:nmax
%     d2(i,i)=d2(i,i)-R(i);
% end

for i=1:nmax
    d2(i,i)=d2(i,i)- (R(i) - avgRx);
end






[eigfcns,eigen]=eig(d2);
eig1=diag(eigen)';
eig2=sort(eig1,'ascend')';

%{
% print out modified eigenvalues
i=1;
fprintf('%2i   %12.6e \n' ,i, eig2(i));
for i=2:imax
    fprintf('%2i   %12.6f \n' ,i, eig2(i));
end
disp('Modified Eigenvalues')
%}


%{
% check if the index of wanted eigenfunction stay unchanged even if with
% R(x)
pause
disp(indxeigf)
pause
for i=1:nmax
    if eig1(i)==eig2(npoly)
        mod_indxeigf=i;
    else
    end
end
disp(mod_indxeigf)
pause
%}


for i=1:nmax
    fac2=sqrt(pt(i)^2*exp(-pt(i)^2)/sqrt(pi));
    fac1=sqrt(wt(i));
    eigf2(i)=eigfcns(i,npoly);
    eigf1(i)=eigfcns(i,indxeigf);
%   fprintf('%9.4f %12.5e %12.4e %12.5e %12.4e %12.4e\n' ,pt(i),...
%       eigf2(i),eigf1(i),eigf1(i)/eigf2(i),fac1,fac2);
%   pause
end

% KEEP THESE???
% eig2(1)
% exp(-eps)


% THIS ETA IS WRONG!!!
eta=(exp(-eps)+eig2(1))/exp(-eps);


% disp('Correction-to-Rx-rate')
% pause
% %fac=eigf2(14)/abs(eigf2(14));

%{
% plot the npoly-th eigenfunctions
plot(pt,eigf1,'-ok','linewidth',1.2,'markersize',7,'markerfacecolor','k')
axis([0 5 -.5 .5])
set(gca,'FontSize',24)
set(gca,'Ytick',[-.4:.2:.4],'linewidth',1.6)
set(gca,'Xtick',[0:1:5],'linewidth',1.6)
str2=num2str(npoly);
str1={'$n = $'};
str=strcat(str1,str2);
text(3,.3,0,str,'Interpreter','latex','fontsize',26)
xlabel('$x$','Interpreter','latex','fontsize',30)
ylabel('$\psi_n(x)$','Interpreter','Latex','fontsize',30)
%}


% forming result
% result.eta = eta;
result.eigenvalues = eig2';
result.eigenfunction = eigf1;
result.quadpt = pt';
result.weight = wt;
result.Rx = R;
result.avgRx = avgRx;


end