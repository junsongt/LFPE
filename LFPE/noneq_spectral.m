function noneq_spectral(nmax, eps)
format long e 
% =========================================================================
%          Code run three times for three different eps values
%          PLOTTING SYMBOL (line 91) CHANGED EACH TIME
%       The exact value of eta computed to 16 sig figs
% =========================================================================
% nmax=16;
% eps=8;
exact=3.6571372511842871e-002;
% For the next run - uncomment the next two lines 
%eps=16;
%exact=1.0178842461950339e-003;
% and for the third run  - uncomment the next two.
%eps=32;
%exact= 1.4278326257466724e-006;
myfile=fopen('noneq.dat','wt');
% The B_j vector in Eq. (13)
b=(0:nmax+1);
b=[0,b];
% 
% MATRIX ELEMENTS of the LINEARIZED Boltzmann collision operator given 
% by Eq. (13) in Gust and Reichl % Phys Rev E79, 031202 (2009)
% with ell = 0, that is the isotropic operator.
%
for n=2:nmax+1
    for m=2:n
        sum=0;
   for jp1=1:(m+1)
        j=jp1-1;
        arg=n+m-2*j-0.5;
% Replace gamma(arg) with
% newgamma=gamma(z)*gamma(z+0.5)/(sqrt(pi)*2^(1-2*z));
        sum=sum+b(jp1)*gamma(arg)/(factorial(n-j)*factorial(m-j)*2^(n+m-2*j));
   end
%   c(n-1,m-1)=-sum*sqrt(factorial(n)*factorial(m)/(8*gamma(n+1.5)*gamma(m+1.5)));
c(n-1,m-1)=-2*sum/sqrt(pi); 
c(m-1,n-1)=c(n-1,m-1);
    end
end
% Polynomial coefficients of the Sonine-Laguerre polynomials
for np1=1:nmax+3
    n=np1-1;
    for ip1=1:np1
        i=ip1-1;
        son(np1,ip1)=(-1)^i*gamma(n+1.5)/(gamma(i+1.5)*factorial(n-i)*factorial(i));
    end
end
es=exp(-eps);
%myfile=fopen('c:\eig12.dat','wt')
jk(1)=es*(eps+1)/2;
for kp1=2:nmax+3
    term=es*(eps^kp1)/2;
    jk(kp1)= es*(eps^kp1)/2+kp1*jk(kp1-1);
end
% The K_k integrals with the reactive cross section; Eq. (5.83)
xk(1)=es/2;
for kp1=2:nmax+3
    k=kp1-1;
    xk(kp1)=jk(kp1)-jk(kp1-1)*eps;
end
% The A_k integrals from the K_k integrals; Eq. (5.82) 
for np1=3:nmax+2
    sa=0;
    for kp1=1:np1
        sa=sa+son(np1,kp1)*xk(kp1);
    end
    a(np1-2)=8*sa/2^(np1-1);
end
nvec=[];
etavec=[];
acc=[];

for nn=1:nmax
    nvec=[nvec nn];
    ax=a(1:nn);
    cx=c(1:nn,1:nn);
    an=linsolve(cx,ax');
    eta=0;
    for kk=1:nn
        eta=eta+an(kk)*a(kk);
    end
    eta=eta/(2*es);
    accuracy=log10(abs(1-eta/exact));
    acc=[acc accuracy];
    etavec=[etavec eta];
end
nvec=[1:1:nmax];

% Symbol on the graph changed for each eps value; also markerfacecolor

plot(nvec,acc,'-^k','linewidth',2.2,'markersize',12,'markerfacecolor','k')
axis([0 16 -10 -2])
set(gca,'FontSize',28)
set(gca,'Ytick',[-10:2:-2],'linewidth',2.2)
set(gca,'Xtick',[0:4:16],'linewidth',2.2)
ylabel('Accuracy','Interpreter','LaTex','FontSize',32)
xlabel('$N$','Interpreter','LaTex','FontSize',32)
leg=legend('$E^{*}/k_BT=8$','$E^{*}/k_BT=16$','$E^{*}/k_BT=32$','Location','SouthWest');
set(leg,'Interpreter','latex','fontsize',26);
str={'(C)'};
text(13.5,-2.5,str,'Interpreter','latex','fontsize',36)
set(gcf, 'Units','centimeters','Papersize',[36,36])
set(gcf, 'Units','centimeters','Position',[3 3 24 20]) 
set(gca,'OuterPosition',[0.1 0.1 .9 0.9])
        