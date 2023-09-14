function linbe_eig_gamma1
disp('Insert npts = 10 - 200');
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
pause
%delete last element of the off-diagonal elements
t=diag(rtbx,-1)+diag(alfx)+diag(rtbx,1);
[f,lambda]=eig(t);
pt=diag(lambda);
wt=sqrt(pi)*f(1,:).^2/4.d0;
wtfcn=pt.^2.*exp(-pt.^2);
rtpi = sqrt(pi);
fac=4/(3*rtpi);
for i=1:npts
    bigwt(i)=wt(i)/wtfcn(i); 
end

rtpi=sqrt(pi);
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
    s=0.;
    for i=1:npts
        x1=pt(i); 
        s=s+bigwt(i)*a(j,i)*pt(i);
    end
    znum(j)=2*s;
    zexact(j)=(2*x2^2+1)*rtpi*erf(x2)/(2*x2)+exp(-x2^2);
end
% out=[pt,znum,zexact];

% Uncomment the next three lines to see comparison with numerical and
% exact collision frequency. The numerical collision frequency ensures that
% one eigenvalue is exactly zero.
for i=1:npts
    fprintf('%10.7f %10.5f %10.5f %10.5f\n' ,pt(i),znum(i),zexact(i),znum(i)/zexact(i));
end
for i=1:npts
    x1=pt(i);
    for j=1:npts
        x2=pt(j);
        bk = 2*wwkern(x1,x2)*sqrt(x1*x2)*exp(-x1^2/2+x2^2/2);
        a(i,j)=sqrt(bigwt(i)*bigwt(j))*bk;
        if i == j
            a(i,j)=a(i,j)-znum(j);
        else
        end
    end
    %    pause
end
[eigfcn,lambda]=eig(a);
eigen=-diag(lambda)/2;
eigen2=sort(eigen);
%eigen2(2)=[];
%fprintf(' %0.7f & %0.7f & %0.7f & %0.7f & %0.7f & %0.7f & %0.7f & %0.7f\\\n',eigen2(1),...
%    eigen2(2),eigen2(3),eigen2(4),eigen2(5),eigen2(6),eigen2(7),eigen2(8))
disp('The first 7 nonzero eigenvalues of the')
disp('linear Boltzmann collision operator for mass ratio 1')
fprintf('N = %i\n',npts)
for i=1:8
    fprintf('%i   %0.7f\n',i-1,eigen2(i)); 
end
disp('Eigenvalues < 1 are in the discrete spectrum whereas')
disp('eigenvalues > 1 are in the continuum')
    function bk = wwkern(x,y)
        % --- FROM THE PAPER BY SIEWERT - JQSRT 74, 789 (2002)
        % --- Eqs. (7c), (8a) and 8(b)
        if(y <= x)
            xkern=erf(y)/x;
            xk=rtpi*xkern;
        elseif (y > x)
            xkern=exp(x^2-y^2)*erf(x)/x;
            xk=rtpi*xkern;
        end
        bk=xk;
    end
end
