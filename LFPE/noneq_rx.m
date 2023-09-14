function result = noneq_rx(nmax)
%The eigenvalues from Gust-Reichl-2009 
format short e
% The B_j vector in Eq. (13)
b=(0:nmax+1);
b=[0,b];
% MATRIX ELEMENTS given by Eq. (13) in Gust and Reichl
% Phys Rev E79, 031202 (2009) with ell = 0, that is the isotropic operator.
for n=2:nmax+1
    for m=2:n
        sum=0;
    for jp1=1:(m+1)
        j=jp1-1;
        arg=n+m-2*j-0.5;
        sum=sum+b(jp1)*gamma(arg)/(factorial(n-j)*factorial(m-j)*2^(n+m-2*j));
    end
%   c(n-1,m-1)=-sum*sqrt(factorial(n)*factorial(m)/(8*gamma(n+1.5)*gamma(m+1.5)));
    c(n-1,m-1)=-2*sum/sqrt(pi); 
    c(m-1,n-1)=c(n-1,m-1);
    end
end
for np1=1:nmax+3
    n=np1-1;
    for ip1=1:np1
        i=ip1-1;
        son(np1,ip1)=(-1)^i*gamma(n+1.5)/(gamma(i+1.5)*factorial(n-i)*factorial(i));
    end
end
% Store eps values in a vector
epsvec=[0:0.1:20];
% Store eta values in a vector
etavec=[];
neps=length(epsvec);
% Loop through the chosen eps values
for nneps=1:neps
    eps=epsvec(nneps);
    es=exp(-eps);
    jk(1)=es*(eps+1)/2;
    for kp1=2:nmax+3
        term=es*(eps^kp1)/2;
        jk(kp1)= es*(eps^kp1)/2+kp1*jk(kp1-1);
    end

    xk(1)=es/2;
    for kp1=2:nmax+3
        k=kp1-1;
        xk(kp1)=jk(kp1)-jk(kp1-1)*eps;
    end

    for np1=3:nmax+2
        sa=0;
        for kp1=1:np1
            sa=sa+son(np1,kp1)*xk(kp1);
        end
        a(np1-2)=8*sa/2^(np1-1);
    end

    an=linsolve(c,a');
    % an
    eta=0;
    for kk=1:nmax
        eta=eta+an(kk)*a(kk);
    end
    eta=eta/(2*es);
    etavec=[etavec eta];
end


result.etas = etavec;

plot(epsvec,etavec,'-k','linewidth',1.6)
hold on


%{
%
% Attach variational solution
%
eps2=[];
eta2=[];
for eps=0:0.5:20
    s = [-1:0.001:0.5];
    fac1 = 2*(1-2*s).^2*exp(-eps)./(s.^4.*sqrt(1-s));
    fac2 = sqrt((2-s)/2).*exp(eps*s./(2-s));
    fac3 = (1-0.75*s+0.5*eps*s)./sqrt(1-s);
    eta = fac1.*(fac2-fac3).^2;
    eta2 = [eta2 max(eta)];
    eps2 = [eps2 eps];
end


plot(eps2,eta2,'ok','markersize',10,'markerfacecolor','k')
axis([0 20 0 0.1])
set(gca,'FontSize',36)
xlabel('$E^*/k_BT$','Interpreter','LaTex','FontSize',36)
ylabel('$\eta$','Interpreter','LaTex','FontSize',36)
set(gca,'Ytick',[0:.02:.1],'linewidth',1.6)
set(gca,'Xtick',[0:5:20],'linewidth',1.6)
%axis([-3 3 -1 1.2])
set(gcf, 'Units','centimeters','Papersize',[36,36])
set(gcf, 'Units','centimeters','Position',[3 3 24 20]) 
set(gca,'OuterPosition',[0.1 0.1 .9 0.9])
%}
end
        