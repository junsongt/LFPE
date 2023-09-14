function SE_electron_potl
% Convergence of the eigenvalues for the Schroedinger equation 
% for the potential in Eq. (6.180) with the pseudospectral solution
% based on quadratures for the weight function Eq. (6.181)
% w(x) = x**5*exp(-x**4/16) - Order, N, of the calculation is in maxvector
nmaxvec=[4 5 6 8 10 15 20 30 40 45 50 60];
kkmax=length(nmaxvec);
% Scan the particular N values above
for kk=1:kkmax
nmax=nmaxvec(kk);
fprintf('nmax = %2i\n',nmax)
%nmax is the number of quad pts; npoly is to plot P_npoly(x)
format long e
% Recurrence coefficients for the polynomials orthogonal wrt to w(x).
load ab_se_electron.dat
a=ab_se_electron(1:nmax,1);
b=ab_se_electron(1:nmax,2);
rtb=sqrt(b);
t=diag(rtb(2:nmax),-1)+diag(a)+diag(rtb(2:nmax),1);
[f,lambda]=eig(t);
% Quadrature points and weights
pt=diag(lambda);
mu_0=8*sqrt(pi); wt=mu_0*f(1,:).^2; mu_1=32*gamma(7/2); mu_2=64;
rtpi=sqrt(pi); omptsq=diag(1-pt.*pt); xn=[0:1:nmax];
% First two P_0 and P_1 polynomials and derivatives P_0' and P_1' 
v=ones(1,nmax); poly=[]; polyp=[]; 
poly=[poly v'./sqrt(mu_0)];poly=[poly (pt-a(1))/sqrt(b(2)*mu_0)];
w = zeros(1,nmax); 
polyp=[polyp w']; polyp=[polyp (v./sqrt(b(2)*mu_0))'];
% The polynomials calculated by recurrence
for n=2:nmax-1
    xp=(pt-a(n)).*poly(:,n)/rtb(n+1)-poly(:,n-1)*rtb(n)/rtb(n+1);
    poly=[poly xp];
    yp=(pt-a(n)).*polyp(:,n)/rtb(n+1)-polyp(:,n-1)*rtb(n)/rtb(n+1)+...
        poly(:,n)/rtb(n+1);
    polyp=[polyp yp];
end
% Calculation of the derivative operator, Eq. (3.138)
for i=1:nmax
    xppi=polyp(i,:);
    for j=1:nmax
        ypj=poly(j,:);
        s=sqrt(wt(i)*wt(j))*sum(ypj.*xppi);
d(i,j)=s;
    end
end
% The representation of the Hamiltonian with the potential Eq. (6.190)
% WITHOUT reference to the potential
mat2=d'*d;
[eigfcns,eigen]=eig(mat2);
eig1=diag(eigen)';
eig2=sort(eig1,'ascend');
imax=nmax;
% Print eigenvalues as reported in Table 6.7
if nmax >=31
    imax=31;
end
for i=2:imax
   fprintf('%2i %12.5f\n' ,i-1,eig2(i));
end
pause
end