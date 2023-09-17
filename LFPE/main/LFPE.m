% Author: Junsong Tang
% Description: This function gives an object which contains:
% the eigenvalues and eigenfunctions of Lorentz Fokker-Planck operator

function result = LFPE(nmax,npoly)
% -- nmax is the number of quad pts; npoly is to plot P_npoly(x)

% -- Load the alpha_n and beta_n recurrence coefficints
% -- Weight function w(x)=x^2exp(-x^2)

% load ab2-70.dat;
% a=ab2_70(1:nmax,1); b=ab2_70(1:nmax+2,2);


% If users dont have "ab2-70.dat", then call ab_maxwell_p2 to produce
% quadrature points
ab = ab_maxwell_p2(nmax,100,100,0,10);
a = ab(:,1);
b = ab(:,2);

rtb=sqrt(b);
rtb(1)=[];
t=diag(rtb(1:nmax-1),-1)+diag(a)+diag(rtb(1:nmax-1),1);
[f,lambda]=eig(t); pt=diag(lambda);
% Quadrature points and weights from diagonlization of Jacobi matrix
mu_0=sqrt(pi)/4; 
wt=mu_0*f(1,:).^2;
v=ones(1,nmax);
w = zeros(1,nmax); 
poly=[]; 
polyp=[]; 
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
    end
end


for i=1:nmax
    for j=1:nmax
        s=0;
        for l=1:nmax
            %s=s+wt(l)*d(l,i)*d(l,j); % this will give Rayleigh eigenvalues!
            s=s+pt(l)*wt(l)*d(l,i)*d(l,j); 
        end
        d2(i,j)=s/sqrt(wt(i)*wt(j));
    end
end

% return a matrix of col vector as eigenvector, and diag matrix as eigenvalues
[eigfcns,eigen]=eig(d2);
% eig1 is a col vector of eigenvalues
eig1=diag(eigen)';
eig2=sort(eig1);
% match the npoly-th eigenvalue to its pos in the unsorted eigenvalues,
% since the index in the unsorted matched the eigenfunction
for i=1:nmax
    if eig1(i)==eig2(npoly)
            indxeigf=i;
    else
    end
end



% pick the npoly-th eigenfunction using index: indxeigf
for i=1:nmax
    fac2=sqrt(pt(i)^2*exp(-pt(i)^2)/sqrt(pi));
    fac1=sqrt(wt(i));
    eigf2(i)=eigfcns(i,npoly);
    eigf1(i)=eigfcns(i,indxeigf);
  % fprintf('%9.4f %12.5e %12.4e %12.5e %12.4e %12.4e\n' ,pt(i),...
  %     eigf2(i),eigf1(i),eigf1(i)/eigf2(i),fac1,fac2);
  % pause
end


% Forming results
result.eigenvalues = eig2;
result.eigenfunction = eigf1;
% quadrature pts
result.quadpt = pt';
% Maxwell basis (column j --> [Mj(x1),...Mj(xn)]^T)
result.poly = poly;

end