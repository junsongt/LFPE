function maxwell_p0_with_pts(npoly)
xmax=10; ngrid=1501; dx=xmax/(ngrid-1); 
x=[0:dx:xmax]; wtfcn=exp(-x.*x); srwtfcn=sqrt(wtfcn);
%grid to plot polynomial and root of the weight function
load abmaxp0.dat; n=100;
disp(abmaxp0) ;                  %delete the bottom m+1 to n rows
alf=abmaxp0(:,1);
pause
b=abmaxp0(:,2);
rtb=sqrt(b); 
rtbx=rtb; rtbx(1)=[]; rtbx(npoly:n-1,:)=[]
alfx=alf; alfx(npoly+1:n,:)=[]
pause
%delete last element of the off-diagonal elements
t=diag(rtbx,-1)+diag(alfx)+diag(rtbx,1);
[f,lambda]=eig(t);
disp('quad pts');
pt=diag(lambda);
pause
for i=1:npoly
wt(i)=sqrt(pi)*f(1,i)^2/2.d0; end
ptwt=[pt,wt'];
zero2=0*ones(1,ngrid); %To draw dashed line at y = 0 
a=[];
gam=sqrt(pi/4); srgam=sqrt(gam); rtpi=sqrt(pi);
b1=-sqrt(4/(rtpi*(pi-2)));
a1=-rtpi*b1;
pause
%Create the first polynomial as a vector with unit components:
p1=sqrt((2/rtpi))*ones(1,ngrid);
%Create the matrix A with the first column x and the second p(1):
a=[a p1'];
%Create the second polynomial as a vector:
p2=a1*x+b1*ones(1,ngrid);
%Add this vector to the 3rd column of A:
%display('B_2(x)')
%plot(x,srwtfcn.*p2,'-k')
%pause
a=[a p2'];
%display('B_2(x)')
%plot(sqrt(x),srwtfcn'.*a(:,2),'-k')
for n=2:npoly+1
%Use the recursion relation to get the next polynomial:
    p3=(x.*p2-alf(n)*p2)/rtb(n+1) - p1*rtb(n)/rtb(n+1);
    %Add this polynomial to the next column of A:    
    a=[a p3'];
%Store the previous two polynomials so as to use in the recursion
%the next time around in the loop
    p1=p2; p2=p3;
end
npt=length(pt); zero=0*ones(1,npt);
maxpolyn=srwtfcn'.*a(:,npoly+1);
plot(x,maxpolyn,'-k','linewidth',1.6)
hold on
plot(pt,zero,'ok','markersize',9,'markerfacecolor','k')
hold on
plot(x,zero2,'--k','linewidth',1.2)
hold off
set(gca,'FontSize',28)
axis([0 8 -0.8 0.8])
set(gca,'Ytick',[-.8:.4:.8],'linewidth',1.6)
set(gca,'Xtick',[0:2:8],'linewidth',1.6)
axis([0 8 -0.8 0.8])
xlabel('$x$','Interpreter','LaTex','FontSize',28)
ylabel('$e^{-x^2/2}M_{10}^{(0)}(x)$','Interpreter','LaTex','FontSize',28)
end

     