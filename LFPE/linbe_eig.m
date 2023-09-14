function linbe_eig
%
% This code can be significantly vectorized
%
clc
disp('Insert npts = 10 - 200')
npts = input('Please enter a value for npts:');
format long e
% Load alpha_n and beta_n for Maxwell quadrature p = 2
load ab_maxp2_200.dat; n=200;
alf=ab_maxp2_200(:,1);
b=ab_maxp2_200(:,2); 
b([1])=[];
rtb=sqrt(b);
rtbx=rtb; 
rtbx(npts:n-1,:)=[];  %delete the bottom npoly to n rows
alfx=alf; 
alfx(npts+1:n,:)=[];
% Delete last element of the off-diagonal elements
% Quadrature points and weights for Maxwell polynomials
t=diag(rtbx,-1)+diag(alfx)+diag(rtbx,1);
[f,lambda]=eig(t); 
pt=diag(lambda); 
wt=sqrt(pi)*f(1,:).^2/4;
wtfcn=pt.^2.*exp(-pt.^2);
rtpi = sqrt(pi);
fac=4/(3*rtpi);
for i=1:npts
    bigwt(i)=wt(i)/wtfcn(i); 
end

% Construct the kernel
for i=1:npts
    x1=pt(i);
    for j=1:npts
        x2=pt(j); 
        bk = wwkern(x1,x2); 
        a(i,j)=bk;
    end
end
% --- This is the kernel in Siewert; Calculate the collision frequency
for j=1:npts
    x2=pt(j); 
    s=0;
    for i=1:npts
        s=s+wt(i)*a(i,j)*pt(i)^2;
    end
    znum(j)=(1/x2^2)*s;
    zexact(j)=(2*x2^2+1)*rtpi*erf(x2)/(2*x2)+exp(-x2^2);
end
% out=[pt,znum,zexact];

% Uncomment the next three lines to see comparison with numerical and
% exact collision frequency. The numerical collision frequency ensures that
% one eigenvalue is exactly zero.
% for i=1:npts
% fprintf('%10.7f %10.5f %10.5f %10.5f\n' ,pt(i),znum(i),zexact(i),znum(i)/zexact(i));
% end
for i=1:npts
    x1=pt(i);
    for j=1:npts
        x2=pt(j);
        bk = wwkern(x1,x2); 
        a(i,j)=sqrt(wt(i)*wt(j))*bk;
        if i == j
            a(i,j)=a(i,j)-znum(j);
        else
        end
    end
end
[eigfcn,lambda]=eig(a);
%pt
%eigfcn(:,1)
pause
eigen=-diag(lambda)/2;
eigen2=sort(eigen);
eigen2(2)=[];
for k=2:npts
    if eigen2(4)==eigen(k)
        kk=k;
    else
    end
end
%fprintf(' %0.6f & %0.6f & %0.6f & %0.6f & %0.6f & %0.6f & %0.6f\\\n',eigen2(1),...
%    eigen2(2),eigen2(3),eigen2(4),eigen2(5),eigen2(6),eigen2(7))
disp('The first 7 nonzero eigenvalues of the')
disp('linearized Boltzmann collision operator')
for i=1:9
    fprintf('%i   %0.6f\n',i,eigen2(i)); 
end
disp('Eigenvalues < 1 are in the discrete spectrum whereas')
disp('eigenvalues > 1 are in the continuum')
subplot(3,1,1)
plot(pt,eigfcn(:,kk),'-k','linewidth',1.6)
hold on
axis([0 2 -.4 .4])
set(gca,'FontSize',26)
set(gca,'Ytick',[-.4:.4:.4],'linewidth',1.6)
set(gca,'Xtick',[0:1:2],'linewidth',1.6)
for k=2:npts
    if eigen2(5)==eigen(k)
        kk=k;
    else
    end
end
subplot(3,1,2)
plot(pt,eigfcn(:,kk),'-k','linewidth',1.6)
axis([0 4 -.4 .4])
set(gca,'FontSize',26)
set(gca,'Ytick',[-.4:.4:.4],'linewidth',1.6)
set(gca,'Xtick',[0:1:2],'linewidth',1.6)
for k=2:npts
    if eigen2(6)==eigen(k)
        kk=k;
    else
    end
end
subplot(3,1,3)
plot(pt,eigfcn(:,kk),'-k','linewidth',1.6)
axis([0 1 -.5 .5])
set(gca,'FontSize',26)
set(gca,'Ytick',[-.5:.5:.5],'linewidth',1.6)
set(gca,'Xtick',[0:.5:1],'linewidth',1.6)
    function bk = wwkern(u,v)
        % --- FROM THE PAPER BY SIEWERT - JQSRT 74, 789 (2002)
        % --- Eqs. (7c), (8a) and 8(b)
        if(u <= v)
            xkern=4*exp(u^2)*erf(u)/(u*v)-fac*(u*u+3*v^2)/v;
            xk=rtpi*xkern;
        elseif (u > v)
            xkern=4*exp(v^2)*erf(v)/(u*v)-fac*(3*u*u+v^2)/u;
            xk=rtpi*xkern;
        end
        bk=xk;
    end
end
