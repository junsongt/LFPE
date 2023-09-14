function linbe_gr_eig(nmax)
% Calculate the eigenvalues of the LINEARIZED Boltzmann collision operator
% MATRIX ELEMENTS given by Eq. (13) in Gust and Reichl
% Phys Rev E79, 031202 (2009) with ell = 0, that is the isotropic operator.
format long e
tic
% Value of 1st eigenalvues with nmax = 96
exact= 6.712273681229485e-01;
% First row of the eigenvalue matrix to identify the first 5 eigenvalues
header = [1 2 3 4 5];
eigmat=[];
% The B_j vector in Eq. (13)
b=(0:nmax+1);
b=[0,b];
xn=[];
acc=[];
% Construct the matrix representation of the collision operator
% increasing the dimension of the matrix by 5 
for n=2:nmax+1
    for m=2:n
        sum=0;
        for jp1=1:(m+1)
            j=jp1-1;
            arg=n+m-2*j-0.5;
            z=arg/2;
            % Replace gamma(arg) with
            % newgamma=gamma(z)*gamma(z+0.5)/(sqrt(pi)*2^(1-2*z));
            sum=sum+b(jp1)*(gamma(z)/(factorial(n-j))*gamma(z+0.5)/(factorial(m-j)*2^(1-2*z-2*j)));
        end
        c(n-1,m-1)=-(0.5)^(n+m)*sum*sqrt(factorial(n)*factorial(m)/(8*gamma(n+1.5)*gamma(m+1.5)))/sqrt(pi);
        c(m-1,n-1)=c(n-1,m-1);
    end
end
% Caculate the eigenvalues
for nmx=5:5:nmax
    xn=[xn nmx];
    e=eig(c(1:nmx,1:nmx));
    % Store the 1st five eigenvalues as the rows of the matrix eigmat
    eigmat=[eigmat; e(1:5)'];
    accnew=log10(abs(1-e(1)/exact));
    acc=[acc accnew];
    % Store the accuracy of the first eigenvalue to plot spectral convergence
end
%xline=eigmat(1,:);
% PRINT THE HEADER LINE OF THE MATRIX
fprintf('%10i %13i %13i %13i %13i \n',header)
% nbf is the number of basis functions
nbf=0;
for i=1:nmax/5
    xline=eigmat(i,:);
    nbf=nbf+5;
    % Print the first 5 eigenvalues
    fprintf('%2i %13.7f %13.7f %13.7f %13.7f %13.7f \n',nbf, xline);
end
toc
% A graphical representation of the spectral convergence
% of the lowest nonzero eigenvalue
plot(xn,acc,'-ok','linewidth',1.6,'markersize',10,'markerfacecolor','k')
xlabel('${\rm N}$','Interpreter','latex','fontsize',26)
ylabel('$\log_{10}[|1 - \lambda_1^{(N)}/\lambda_1^{(Exact)}|]$','Interpreter','latex','fontsize',26)
axis([0 80 -14 -2])
set(gca,'FontSize',26)
set(gca,'Ytick',[-14:2:-2],'linewidth',1.6)
set(gca,'Xtick',[0:10:80],'linewidth',1.6)
set(gcf, 'Units','centimeters','Papersize',[36,36])
set(gcf, 'Units','centimeters','Position',[3 3 24 20])
%fprintf(myfile,'%i & %16.12f & %16.12f & %9.5f & %9.5f & %9.5f  & %9.5f\n',nmax,e(2),e(3),e(4),e(5),e(6),e(7));
