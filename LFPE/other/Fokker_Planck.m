%Author: Lucas Philipp
%Description: this code numerically diagonalizes the hermitian
%Fokker-Planck operator L (eigenvalues and eigenfunctions)
% =====================================================================

%change B(x) in the Lij expression
%choose w(x) = P_0(x)

function Fokker_Planck(xmin, xmax, numbasis, nint, npts)
%xmin/xmax defines the domain of integration. the function can be truncated
%at a xmin if beyond which the function can be considered negligable, similarly with xmax
%numbasis is the dimension number representing the total number of basis polynomials, as well as the number of quadrature points/weights 
%numeig specifies the eigenfunction # that you wish to plot
%nint is the number of subintervals used exclusively in the Fejer quadrature subroutine
%npts is the number of points per subinterval used exclusively in the Fejer quadrature subroutine

%Gautschi Stieltjes procedure for caluation of recurrence coefficients alpha and beta
tic %start of first timer (can guage performance)
format long e %format data type (e means scientific notation)

pwcat = multidomquad(nint,npts,xmin,xmax); %multi-dimensional matrix of quad pts p and wts w
%this calls functions multidomquad. see multidomquad.

%Fejer quadrature used to numerically evaluate integrals and to initialize
%the three term recurrence relation
ntot=nint*npts; %total number of Fejer quadrature points in [xmin, xmax]
p=pwcat(:,1); %quadrature points are the first column (MATLAB has no 0 index)
Nzeros=pwcat(:,2); %quadrature weights are the second column (MATLAB has no 0 index)

%SPECIFY THE WEIGHT FUNCTION HERE
thermobeta = 1.080;
a = 460;
b = -48;
c = -115;
wtfcn=(1./(thermobeta)).*exp(-(1./thermobeta).*((a.*p.^4)+(b.*p.^3)+(c.*p.^2)));

%orthogonal polynomial Q_0(x)
q0=ones(ntot,1); %returns an ntot by 1 array of ones (Q_0(x)=1, Q_-1=0)
%evaluate using Fejer quadrature
s1=sum(Nzeros.*wtfcn); %gamma_1 i.e. the first norm
s2=sum(Nzeros.*(p.*wtfcn)); %<Q_0(x)|x|Q_0(x)>
mu_0=s1; %the 0th moment with x^n, n=0 is s1
gamma(1)=s1; %gamma stores the norm values %(MATLAB has no 0 index) gamma(1) is really more like gamma(0)
alpha(1)=s2/s1; %(MATLAB has no 0 index) alpha(1) is really more like alpha(0)
beta(1)=0; %(MATLAB has no 0 index) beta(1) is really more like beta(0). Since there is no beta(0) set matrix element to 0
k=1; %index variable

%print alpha's and beta's as a list
disp(' ')
disp('alpha''s and beta''s')
fprintf('%20.12f %20.12f\n',alpha(k),beta(k)); %\n starts a new line
%20.12f means first value in each line is floating-point with field width of 20 digits, 12 past the decimal

%orthogonal polynomial Q_1(x)
q1=p-alpha(1); %using the three term recurrence relation
%Next two integrals 
s1=sum(Nzeros.*(wtfcn.*(q1.^2)));
s2=sum(Nzeros.*(p.*(wtfcn.*(q1.^2))));
gamma(2)=s1; %(MATLAB has no 0 index) alpha(2) is really more like alpha(1)
alpha(2)=s2/s1; %(MATLAB has no 0 index) alpha(2) is really more like alpha(1)
beta(2)=gamma(2)/gamma(1); %(MATLAB has no 0 index) beta(2) is really more like beta(1)
k=2; %increment index variable

%write alpha(1), beta(1) to text file
fprintf('%20.12f %20.12f\n',alpha(k),beta(k));

%iterate using for loop
for k=3:numbasis+1
    pma=p-alpha(k-1); %p minus alpha
    %Recurrence for the remaining polynomials (up to numbasis)
    q2=pma.*q1-beta(k-1)*q0;
    s1=sum(Nzeros.*(wtfcn.*(q2.^2)));
    s2=sum(Nzeros.*(p.*(wtfcn.*(q2.^2))));
    alpha(k)=s2/s1; gamma(k)=s1; beta(k)=gamma(k)/gamma(k-1);
    %alpha_k, norm and beta_k
    fprintf('%20.16f %20.16f\n',alpha(k),beta(k));
    q0=q1; q1=q2; %increment polynomials
end
disp(' ')
disp('recurrence relation timer')
toc %end of first timer
tic %start of second timer

% Get the quadrature points and weights
rtbeta=sqrt(beta); %square root of all betas
rtbeta(1)=[]; %deletes the first column of the matrix i.e. beta(0) which is just a place holder. rtbeta runs from rtbeta(2), rtbeta(3), ... rtbeta(N)
%which is really more like rtbeta(1), rtbeta(2), ... rtbeta(N-1)
J=diag(rtbeta(1:numbasis-1),-1)+diag(alpha(1:numbasis))+diag(rtbeta(1:numbasis-1),1); %construct the Jacobi matrix %second argument > 0 above main diagonal, < 0 below main diagonal
size(J)
[reigvector,lambda]=eig(J); %lambda is the diagonal matrix containing eigenvalues, reigvector is a matrix whose columns are the right eigenvectors s.t. Jdiag = reigvector J reigvector^-1
pt=diag(lambda); %creates a row vector of the eigenvalues, the quadrature points are the eigenvalues of the Jacobi matrix
disp(' ')
disp('quad pts and wts')
for i=1:numbasis %i is an index
wt(i)=mu_0*reigvector(1,i)^2; %first element of ith eigenvector^2 the first moment gives the weights %wt is a row vector
fprintf('%2i %13.5e %13.5e\n', i, pt(i), wt(i)) %prints quad points and weights in list format
end %end of for loop
disp(' ')
disp('quad pts and wts timer')
toc %end of second timer

%T_ni transformation matrix
tic %start of third timer
%The orthogonal polynomials are evaluated at the quadrature points (eigenvalues of the Jacobi matrix)
Nones=ones(1,numbasis); %row vector with N ones
poly=[]; %initialize polynomial matrix
polyprime=[]; %initialize polynomial derivative matrix
poly=[poly Nones./sqrt(mu_0)]'; %poly is a column vector, mu_0 = =gamma_1 which is more like gamma_0, P_0(x) is caculated in poly 
T(:,1)=poly(:,1).*sqrt(wt'); %colon represents changing the entire column, T has 1 column
poly=[poly (pt-alpha(1))/sqrt(beta(2)*mu_0)]; %P_1(x) is calculated in poly, alpha(1) is more like alpha(0), beta(2) is more like beta(1), mu_0 = = gamma_1 which is more like gamma_0
Nzeros = zeros(1,numbasis);
polyprime=[polyprime Nzeros']; %polyprime is a column vector %Nzeros in the first column which is taking the derivative of P_0(x)
polyprime=[polyprime (Nones./sqrt(beta(2)*mu_0))']; %taking the derivative of P_1(x)
for n=2:numbasis
    xp=(pt-alpha(n)).*poly(:,n)/rtbeta(n)-poly(:,n-1)*rtbeta(n-1)/rtbeta(n); %calculate next P_n(x) %the rtbeta(n) index is shifted 1 relative to alpha(n) and beta(n)
    poly=[poly xp]; %append next P_n(x)
    yp=(pt-alpha(n)).*polyprime(:,n)/rtbeta(n)-polyprime(:,n-1)*rtbeta(n-1)/rtbeta(n)+...
        poly(:,n)/rtbeta(n); %calculate next P'_n(x)
    polyprime=[polyprime yp]; %append next P'_n(x)
    T(:,n)=poly(:,n)'.*sqrt(wt); %append additional columns in T
end
disp(' ')
disp('orthogonal basis polynomials and transformation matrix timer')
toc %end of third timer

%these are the D's without the hat
tic %start of fourth timer
%first order derivative matrix
for i=1:numbasis
    for j=1:numbasis
        s=sum(poly(j,:).*polyprime(i,:)); %sum over P'_n(x_i)P_n(x_j)
d(i,j)=sqrt(wt(i)*wt(j))*s;
    end
end

%second order derivative matrix %I HAVE A PROLEM WITH THIS LINE
d2=d'*d;
disp(' ')
disp('derivative matrix timer')
toc %end of fourth timer

tic % start of fifth timer
[Lfcns,Leig]=eig(thermobeta.*d2); %B(x_k) %Leig are the eigenvalues, Lfcns are the right eigenfunctions
disp(' ')
disp('diagonalization timer')
toc % end of fifth timer

%Evaluate the equilibrium distribution at the quadrature points, pt.
%w(x)=P_0(x)
Pnaught=(1./thermobeta).*exp(-(1./thermobeta).*((a.*pt.^4)+(b.*pt.^3)+(c.*pt.^2)));
plot(pt, Pnaught./885.656,'k','linewidth',1.6)
axis([-0.6 0.6 -0.1 9]) %adjust axis here, try xmin and xmax
set(gca,'FontSize',20)
set(gca,'Xtick',[-0.6:0.2 :0.6],'linewidth',2.0) %set tick marks
set(gca,'Ytick',[0:1:9],'linewidth',2.0)
xlabel('$x$','Interpreter','latex','fontsize',28) %LaTeX independent variable
ylabel('$P_{eq}(x)$','Interpreter','Latex','fontsize',28) %LaTeX dependent variable
pause

%display eigenvalues
roweig=diag(Leig);
roweig=sortrows(roweig);
disp(' ')
disp('eigenvalues')
for i=1:numbasis
fprintf('%2i %20.16e\n', i-1, roweig(i)) %i in the print command means integer, e means scientific notation
end

%extract the numeig'th eigenfunction
for i=1:numbasis
    indexedeigf2(i)=Lfcns(i,1+2); %MATLAB has no zero index
    indexedeigf4(i)=Lfcns(i,1+4);
    indexedeigf6(i)=Lfcns(i,1+6);
    indexedeigf8(i)=Lfcns(i,1+8);
end

figure(1);
subplot(2,2,1);
%plot the numeig'th eigenfunction
plot(pt,indexedeigf2,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([-0.6 0.6 -0.4 0.2]) %adjust axis here, try xmin and xmax
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.1:1],'linewidth',1.6)
set(gca,'Xtick',[-1:0.25:1],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_2(x)$','Interpreter','Latex','fontsize',20)


subplot(2,2,2);
plot(pt,indexedeigf4,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([-0.6 0.6 -0.3 0.3]) %adjust axis here, try xmin and xmax
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.1:1],'linewidth',1.6)
set(gca,'Xtick',[-1:0.25:1],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_4(x)$','Interpreter','Latex','fontsize',20)

subplot(2,2,3);
plot(pt,indexedeigf6,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([-0.6 0.6 -0.3 0.3]) %adjust axis here, try xmin and xmax
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.1:1],'linewidth',1.6)
set(gca,'Xtick',[-1:0.25:1],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_6(x)$','Interpreter','Latex','fontsize',20)

subplot(2,2,4);
plot(pt,indexedeigf8,'-ok','linewidth',1.2,'markersize',3,'markerfacecolor','k')
axis([-0.6 0.6 -0.3 0.3]) %adjust axis here, try xmin and xmax
set(gca,'FontSize', 16)
set(gca,'Ytick',[-1:0.1:1],'linewidth',1.6)
set(gca,'Xtick',[-1:0.25:1],'linewidth',1.6)
xlabel('$x$','Interpreter','latex','fontsize',24) %LaTeX independent variable
ylabel('$\psi_8(x)$','Interpreter','Latex','fontsize',20)
end

% ====================================================================
function pwcat = multidomquad(nint,npts,xmin,xmax)
format long e        
ntot=nint*npts; %calculates the total number of quadrature points
%there are nint intervals and npts per interval
dx=(xmax-xmin)/ntot; %calculate distance between quadrature points
for i=1:nint
    a=xmin+(i-1)*npts*dx; %left end of subinterval 
    b=a+npts*dx; %right end of subinterval
    pw=fejer(a,b,npts); %Fejer quadrature used for each interval
    if i==1
        pwcat=pw;
    else
        pwcat=cat(1,pwcat,pw); %concatenates along the first dimension i.e. under the double column vector
    end
end
end
% =====================================================================

%Fejer quadrature subroutine (according to Fejer's 1st rule)
%pw is matrix with two coloumn vectors, col 1 quad points col 2 quad weights
function pw=fejer(a,b,npts) %generates the n-point quadrature in [a,b]
format long e
i=npts:-1:1; %i is the quadrature index
%backwards construction of the array (total N-1 entries)
j=1:floor(npts/2); %j is the sum index

theta=(2*i-1)*pi./(2*npts); %theta is a row vector containing all indexed thetas

x=cos(theta'); %x are the quadrature points %column vector of quadrature points (' denotes the transpose)
for i_loop=npts:-1:1 %same as i just defined within the loop
  s=sum(cos(2*j*theta(i_loop))./(4*(j.^2)-1)); %handle the sum separately
  w(i_loop)=2*(1-2*s)/npts; %calculate the weights
end
%Mapping [-1,1] onto [a,b]. 
%a is the left bound of the subinterval, b is the right bound of the subinterval
%the subinterval begins at a on the domain and ends at b, b-a is the length
r1=(b-a)/2.0; %size scaling %multiply by half the length because the interval grows in the negative direction %want r1 to be a double.
r2=((b-a)/2.0)+a; %location shift %shift by half the length to shift out of the negative region, then shift by
sclqd=r1*x+r2; %sclqd are the location of the scaled quadrature points on the entire domain
sclwts=r1*w; %sclwts are the scaled weights (multiplied only by the scaling factor)
pw=[sclqd,sclwts']; %sclqd is a column vector, sclwts is a row vector
end